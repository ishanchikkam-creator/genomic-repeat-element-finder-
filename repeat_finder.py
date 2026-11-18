"""Core algorithms for the Genomic Repeat Finder.

Dispersed repeat finder (Tanishka Kumar) and tandem repeat finder
(Ishan Chikkam), combined into one module with GFF3 export.
"""
import io

try:  # Biopython is optional; a built-in parser is used if it is missing.
    from Bio import SeqIO
    HAVE_BIOPYTHON = True
except ImportError:  # pragma: no cover
    HAVE_BIOPYTHON = False


# --------------------------------------------------------------------------
# FASTA input
# --------------------------------------------------------------------------
def _parse_fasta_builtin(text):
    records, rid, chunks = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if rid is not None:
                records.append((rid, "".join(chunks)))
            rid, chunks = line[1:].split()[0] if line[1:].split() else "unnamed", []
        else:
            chunks.append(line)
    if rid is not None:
        records.append((rid, "".join(chunks)))
    return records


def parse_sequences(text):
    """Return [(id, sequence)] from FASTA text or a bare pasted sequence."""
    text = text.strip()
    if not text:
        return []
    if not text.startswith(">"):  # plain pasted sequence
        seq = "".join(text.split())
        return [("pasted_sequence", seq)]
    if HAVE_BIOPYTHON:
        return [(r.id, str(r.seq)) for r in SeqIO.parse(io.StringIO(text), "fasta")]
    return _parse_fasta_builtin(text)


# --------------------------------------------------------------------------
# Dispersed repeats
# --------------------------------------------------------------------------
def find_dispersed_repeats(seq, k_len, min_dist):
    """Find k-mers that occur more than once, at least min_dist apart."""
    if k_len < 1:
        raise ValueError("k-mer length must be at least 1")
    if min_dist < 1:
        raise ValueError("Minimum distance must be at least 1")
    seq = seq.upper()
    repeats = {}
    for i in range(len(seq) - k_len + 1):
        kmer = seq[i:i + k_len]
        if "N" in kmer:
            continue
        repeats.setdefault(kmer, []).append(i)

    dispersed = {}
    for kmer, positions in repeats.items():
        if len(positions) > 1:
            valid = [positions[0]]
            for pos in positions[1:]:
                if pos - valid[-1] >= min_dist:
                    valid.append(pos)
            if len(valid) > 1:
                dispersed[kmer] = valid
    return dispersed


def summarize(results):
    """Summary statistics for a dispersed-repeat result dict."""
    if not results:
        return {"unique_repeats": 0, "total_occurrences": 0,
                "most_frequent": None, "most_frequent_count": 0}
    best = max(results, key=lambda k: len(results[k]))
    return {
        "unique_repeats": len(results),
        "total_occurrences": sum(len(v) for v in results.values()),
        "most_frequent": best,
        "most_frequent_count": len(results[best]),
    }


# --------------------------------------------------------------------------
# Tandem repeats
# --------------------------------------------------------------------------
def _is_periodic(unit):
    """True if unit is itself a repeat of a shorter unit."""
    n = len(unit)
    return any(n % d == 0 and unit == unit[:d] * (n // d) for d in range(1, n))


def find_tandem_repeats(seq, min_unit_len=1, max_unit_len=6, min_repeats=2):
    """Find consecutive copies of a motif."""
    if min_unit_len < 1 or max_unit_len < min_unit_len:
        raise ValueError("Invalid repeat unit length range")
    if min_repeats < 2:
        raise ValueError("Minimum repeats must be at least 2")
    seq = seq.upper().strip()
    repeats = []
    for unit_len in range(min_unit_len, max_unit_len + 1):
        for i in range(len(seq) - unit_len + 1):
            unit = seq[i:i + unit_len]
            if "N" in unit or _is_periodic(unit):
                continue
            if i >= unit_len and seq[i - unit_len:i] == unit:
                continue
            count, j = 1, i + unit_len
            while seq[j:j + unit_len] == unit:
                count += 1
                j += unit_len
            if count >= min_repeats:
                repeats.append({
                    "start": i + 1,
                    "end": j,
                    "motif": unit,
                    "repeat_count": count,
                    "repeat_length": count * unit_len,
                    "sequence": seq[i:j],
                })
    return sorted(repeats, key=lambda r: (r["start"], r["end"], r["motif"]))


# --------------------------------------------------------------------------
# GFF3 export
# --------------------------------------------------------------------------
def to_gff3(analyses, k_len):
    lines = ["##gff-version 3"]
    for a in analyses:
        lines.append(f"##sequence-region {a['id']} 1 {a['length']}")
    n = 0
    for a in analyses:
        for kmer, positions in a["dispersed"].items():
            for p in positions:
                n += 1
                lines.append("\t".join([
                    a["id"], "GenomicRepeatFinder", "dispersed_repeat",
                    str(p + 1), str(p + k_len), ".", ".", ".",
                    f"ID=disp{n};Name={kmer};copies={len(positions)}"]))
        for t in a["tandem"]:
            n += 1
            lines.append("\t".join([
                a["id"], "GenomicRepeatFinder", "tandem_repeat",
                str(t["start"]), str(t["end"]), ".", ".", ".",
                f"ID=tan{n};Name={t['motif']};copies={t['repeat_count']}"]))
    return "\n".join(lines) + "\n"