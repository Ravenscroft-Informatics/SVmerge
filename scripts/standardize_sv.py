#!/usr/bin/env python3
"""
Stage 0.5 post-processing: harmonize REF/ALT/SVLEN/END representation
across callers (Manta, Dysgu, ...) *before* svdb merges them, so that
(a) svdb's bnd_distance/overlap matching compares like-for-like SVLEN/END
    instead of caller-specific artifacts, and
(b) the final cohort VCF has a uniform representation that plays nicely
    with both symbolic-ALT-expecting annotation tools (AnnotSV, VEP,
    SVAFotate) and truvari benchmarking.

Why not just harmonize everything to one style?
    Truvari can only recover sequence for symbolic <DEL>/<INV> (from
    --reference) and <DUP> (via --dup-to-ins); it cannot recover a novel
    insertion's sequence from the reference. So collapsing INS to a
    symbolic <INS> would silently disable truvari's sequence-similarity
    scoring (PctSeqSimilarity) for every insertion. This script therefore:

      - DEL / DUP / INV, size >= --symbolic-threshold (default 50bp):
            rewritten to symbolic ALT (<DEL>, <DUP>, <INV>), REF trimmed
            to the single anchor base, SVLEN/END recomputed and made
            consistent.
      - INS (any size):
            kept sequence-resolved (REF=anchor base, ALT=anchor+insert),
            but rewritten to the *minimal* VCF-spec form. Callers like
            Manta sometimes emit a CIGAR (e.g. 1M376I4D) describing a
            multi-base match/microdeletion alongside the insertion,
            which bloats REF/ALT beyond the clean single-base-anchor
            form Dysgu already uses. When CIGAR is present, the true
            inserted sequence is extracted by walking the CIGAR ops
            over ALT (M/I consume ALT, D does not); any co-located
            deletion length is preserved in INFO/CIGAR_DEL_BP rather
            than silently dropped, since it's a real (if small)
            reference-length change the caller's assembly implied.
            If the CIGAR's M+I bases don't add up to the actual ALT
            length (seen in real Manta output -- the CIGAR and REF/ALT
            can be internally inconsistent), the record is left
            untouched and flagged INFO/SV_STD_SKIPPED=CIGAR_ALT_MISMATCH
            rather than guessing at a possibly-truncated sequence.
      - Everything else (small DEL/DUP/INV below threshold, BND, symbolic
        input already): passed through, only SVLEN/END are verified /
        recomputed for consistency.

If --primary-contigs is given, records are also dropped when their *mate*
falls outside that set: BND ALT bracket notation (e.g.
A[chrUn_JTFH01000016v1_decoy:1[) is parsed for the mate contig, and
INFO/CHR2 is checked for TRA/other records that carry it. The upstream
normalize_caller_vcf rule already restricts each record's own CHROM to
primary contigs, but that doesn't catch a BND/TRA whose own CHROM is
primary while its mate points into a decoy/unplaced contig -- exactly
the kind of alignment artifact decoys exist to absorb, seen in real
Manta output (a chrX BND with a mate on chrUn_JTFH01000016v1_decoy).

Usage:
    python3 standardize_sv.py input.vcf[.gz] output.vcf \
        [--symbolic-threshold 50] [--primary-contigs chr1,chr2,...]
"""
import argparse
import re
import sys

import pysam

CIGAR_OP_RE = re.compile(r"(\d+)([MIDNSHP=X])")
# Matches the contig in BND ALT bracket notation, any orientation:
# A[chr2:123[, ]chr2:123]A, A]chr2:123]A, [chr2:123[A, etc.
BND_MATE_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")
SYMBOLIC_TYPES = {"DEL", "DUP", "INV"}


def parse_cigar(cigar):
    """Return list of (op, length) tuples, or [] if unparseable.

    pysam returns INFO fields declared Number=. or Number=A as a tuple
    even when there's a single value (e.g. CIGAR=('1M376I4D',) for
    Manta, since CIGAR is nominally per-ALT-allele). Unwrap that before
    regex-matching; other multi-value cases (true multi-allelic CIGAR)
    are conservatively treated as unparseable rather than guessed at.
    """
    if not cigar:
        return []
    if isinstance(cigar, tuple):
        if len(cigar) != 1 or not cigar[0]:
            return []
        cigar = cigar[0]
    if not isinstance(cigar, str):
        return []
    ops = CIGAR_OP_RE.findall(cigar)
    if not ops:
        return []
    return [(op, int(length)) for length, op in ops]


def extract_insertion_from_cigar(alt: str, cigar_ops):
    """
    Walk CIGAR ops over the ALT string to recover the true inserted
    sequence(s) (I ops consume ALT; M ops consume ALT but are matches,
    not insert; D ops consume neither, i.e. reference-only bases).

    Returns (inserted_seq, deleted_bp, consistent) where deleted_bp is the
    total length of D ops (reference bases the caller's alignment implied
    were deleted alongside the insertion), and consistent is False if the
    CIGAR's M+I bases don't exactly account for the full ALT string -- in
    that case the caller's CIGAR and its own reported REF/ALT have gone
    out of sync (seen in real Manta output) and the slice is not
    trustworthy; callers of this function should not use inserted_seq
    when consistent is False.
    """
    alt_pos = 0
    inserted_chunks = []
    deleted_bp = 0
    for op, length in cigar_ops:
        if op in ("M", "=", "X"):
            alt_pos += length
        elif op == "I":
            inserted_chunks.append(alt[alt_pos : alt_pos + length])
            alt_pos += length
        elif op == "D":
            deleted_bp += length
        # N, S, H, P: not expected from these callers' SV CIGARs; ignore.
    consistent = alt_pos == len(alt)
    return "".join(inserted_chunks), deleted_bp, consistent


def get_svtype(rec):
    svtype = rec.info.get("SVTYPE")
    return svtype


def is_symbolic_alt(alt: str) -> bool:
    return alt.startswith("<") and alt.endswith(">")


def get_mate_contig(rec, alt, svtype):
    """
    Best-effort extraction of the "other side" contig for records that
    describe a connection between two loci, so it can be checked against
    the primary contig set. Returns None when not applicable/determinable
    (most records, including plain DEL/DUP/INV/INS, have no mate).
    """
    if svtype == "BND" and alt:
        m = BND_MATE_RE.search(alt)
        if m:
            return m.group(1)
        return None
    # Guard with "in rec.info" first: pysam needs a header declaration to
    # decode a value, so calling .get() on a key the header never declared
    # raises ValueError("Invalid header") rather than returning None --
    # e.g. Manta's own VCFs never declare INFO/CHR2 at all (only Dysgu
    # emits it), so a bare rec.info.get("CHR2") blows up on Manta input.
    if "CHR2" not in rec.info:
        return None
    chr2 = rec.info.get("CHR2")
    if chr2 is None:
        return None
    if isinstance(chr2, tuple):
        chr2 = chr2[0] if chr2 else None
    return chr2 if isinstance(chr2, str) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("output")
    ap.add_argument(
        "--symbolic-threshold",
        type=int,
        default=50,
        help="Minimum SVLEN (bp) at which DEL/DUP/INV are rewritten to "
        "symbolic ALT. Below this, sequence-resolved representation is "
        "kept as-is (default: 50).",
    )
    ap.add_argument(
        "--primary-contigs",
        default=None,
        help="Comma-separated primary contig list. If given, records "
        "whose mate (BND ALT bracket target, or INFO/CHR2) falls outside "
        "this set are dropped as likely decoy/unplaced-mate artifacts. "
        "If omitted, no mate-contig filtering is applied.",
    )
    args = ap.parse_args()

    primary_contigs = None
    if args.primary_contigs:
        primary_contigs = {c for c in args.primary_contigs.split(",") if c}

    vin = pysam.VariantFile(args.input)

    header = pysam.VariantHeader()
    for line in str(vin.header).rstrip("\n").split("\n"):
        if not line or line.startswith("#CHROM"):
            continue
        header.add_line(line)
    for s in vin.header.samples:
        header.add_sample(s)

    if "END" not in header.info:
        header.info.add("END", 1, "Integer", "End position of the variant")
    if "SVLEN" not in header.info:
        header.info.add("SVLEN", ".", "Integer", "Length of structural variant")
    if "CIGAR_DEL_BP" not in header.info:
        header.info.add(
            "CIGAR_DEL_BP",
            1,
            "Integer",
            "Reference bases deleted alongside an insertion, per the "
            "caller's CIGAR, that were folded into a clean INS record "
            "by standardize_sv.py rather than emitted as a separate "
            "deletion.",
        )
    if "SV_STD_SKIPPED" not in header.info:
        header.info.add(
            "SV_STD_SKIPPED",
            1,
            "String",
            "standardize_sv.py left REF/ALT untouched for this record "
            "because it could not safely rewrite it. CIGAR_ALT_MISMATCH: "
            "the record's own CIGAR implies a different ALT length than "
            "the ALT actually present (caller's CIGAR and REF/ALT are "
            "internally inconsistent), so the CIGAR-derived insertion "
            "sequence was not trusted.",
        )

    vout = pysam.VariantFile(args.output, "w", header=header)

    n_ins_rewritten = n_symbolic_rewritten = n_passthrough = n_cigar_mismatch = 0
    n_dropped_nonprimary_mate = 0

    for rec in vin:
        svtype = get_svtype(rec)
        alt = rec.alts[0] if rec.alts else None
        ref = rec.ref

        if primary_contigs is not None:
            mate_contig = get_mate_contig(rec, alt, svtype)
            if mate_contig is not None and mate_contig not in primary_contigs:
                n_dropped_nonprimary_mate += 1
                continue

        cigar_ops = parse_cigar(rec.info.get("CIGAR")) if "CIGAR" in rec.info else []

        new_ref = ref
        new_alt = alt
        new_svlen = rec.info.get("SVLEN")
        new_end = rec.info.get("END")
        cigar_del_bp = None
        std_skipped = None

        if alt is None or is_symbolic_alt(alt) or svtype == "BND":
            # Already symbolic or a breakend: leave REF/ALT alone, just
            # make sure END is populated when missing -- required by
            # truvari and most annotation tools to size a symbolic SV.
            if new_end is None:
                if svtype in SYMBOLIC_TYPES and new_svlen is not None:
                    # DEL/DUP/INV span [POS, POS+SVLEN-1].
                    svlen_val = new_svlen[0] if isinstance(new_svlen, tuple) else new_svlen
                    new_end = rec.pos + abs(svlen_val) - 1
                elif svtype in ("INS", "TRA"):
                    # Single-point events: no reference span, END == POS
                    # (matches how both Manta and Dysgu report these when
                    # they do include END, e.g. the Dysgu INS/TRA examples
                    # above both had END equal to POS).
                    new_end = rec.pos
                # BND: no fill -- a breakend is a single point and doesn't
                # carry a meaningful END per the VCF spec.
            n_passthrough += 1

        elif svtype == "INS":
            if cigar_ops:
                inserted_seq, deleted_bp, consistent = extract_insertion_from_cigar(
                    alt, cigar_ops
                )
                if not consistent:
                    # The CIGAR's M+I bases don't add up to the actual ALT
                    # length -- this caller's CIGAR and REF/ALT have gone
                    # out of sync (observed in real Manta INS records).
                    # Don't guess: leave REF/ALT as reported, keep the
                    # caller's own SVLEN/END, and flag it so these records
                    # are easy to find and inspect rather than silently
                    # carrying a truncated/wrong insertion sequence.
                    std_skipped = "CIGAR_ALT_MISMATCH"
                    n_cigar_mismatch += 1
                elif inserted_seq:
                    new_ref = ref[0]
                    new_alt = new_ref + inserted_seq
                    new_svlen = len(inserted_seq)
                    new_end = rec.pos
                    if deleted_bp:
                        cigar_del_bp = deleted_bp
            if new_alt == alt and not cigar_ops:
                # No CIGAR (e.g. Dysgu): already minimal spec form in
                # practice, but recompute SVLEN/END for consistency.
                new_svlen = len(alt) - len(ref)
                new_end = rec.pos
            n_ins_rewritten += 1

        elif svtype in SYMBOLIC_TYPES:
            # Recompute SVLEN from CIGAR (if present) or REF/ALT length
            # diff, then decide symbolic vs sequence-resolved by size.
            if cigar_ops:
                span = sum(length for op, length in cigar_ops if op in ("M", "D", "=", "X"))
                svlen_val = span if span else abs(len(alt) - len(ref))
            else:
                svlen_val = abs(len(ref) - len(alt)) if alt else abs(
                    (new_svlen[0] if isinstance(new_svlen, tuple) else new_svlen) or 0
                )
            if svlen_val >= args.symbolic_threshold:
                new_ref = ref[0]
                new_alt = f"<{svtype}>"
                new_svlen = svlen_val
                new_end = rec.pos + svlen_val - 1
                n_symbolic_rewritten += 1
            else:
                new_svlen = svlen_val
                new_end = rec.pos + svlen_val - 1
                n_passthrough += 1
        else:
            n_passthrough += 1

        new_rec = vout.new_record(
            contig=rec.contig,
            start=rec.start,
            stop=rec.start + len(new_ref),
            alleles=(new_ref, new_alt) if new_alt else (new_ref,),
            id=rec.id,
            qual=rec.qual,
            filter=list(rec.filter.keys()) if rec.filter else None,
        )
        for key, value in rec.info.items():
            if key in ("SVLEN", "END"):
                continue
            try:
                new_rec.info[key] = value
            except (TypeError, OverflowError):
                pass
        if new_svlen is not None:
            new_rec.info["SVLEN"] = new_svlen
        if new_end is not None:
            # pysam reserves INFO/END: it's not settable via info[...],
            # only via the .stop property (same underlying value/integer
            # semantics as VCF END).
            new_rec.stop = new_end
        if cigar_del_bp:
            new_rec.info["CIGAR_DEL_BP"] = cigar_del_bp
        if std_skipped:
            new_rec.info["SV_STD_SKIPPED"] = std_skipped

        for s in rec.header.samples:
            for key in rec.samples[s].keys():
                try:
                    new_rec.samples[s][key] = rec.samples[s][key]
                except TypeError:
                    pass

        vout.write(new_rec)

    vin.close()
    vout.close()

    print(
        f"[standardize_sv] INS rewritten={n_ins_rewritten}, "
        f"symbolic-rewritten={n_symbolic_rewritten}, passthrough={n_passthrough}, "
        f"CIGAR/ALT mismatch (left as-is)={n_cigar_mismatch}, "
        f"dropped (non-primary mate contig)={n_dropped_nonprimary_mate}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()