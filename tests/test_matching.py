"""
Check find_bead (STEP 3 of slidetags_flex.py) against the Julia approach (a table of every
1-mismatch neighbour of every bead, spatial-count.jl create_SBtoindex), and the UMI encoding
against Julia's UMItoindex. Beads with N are removed in STEP 2, so the whitelist has none.

Usage: python tests/test_matching.py
"""
import random
import sys

sys.path.insert(0, ".")
from slidetags_flex import find_bead, umi_to_number  # noqa: E402


def julia_table(whitelist):
    """Literal port of create_SBtoindex: key -> (index, mismatch position); index -2 = ambiguous."""
    table = {}
    def put(key, i, pos):
        table[key] = (-2, 0) if key in table else (i, pos)
    for i, sb in enumerate(whitelist):
        n = sb.count("N")
        if n == 0:
            for j in range(len(sb)):
                for c in "ACGTN":
                    if c != sb[j]:
                        put(sb[:j] + c + sb[j + 1:], i, j + 1)
        elif n == 1:
            k = sb.index("N")
            put(sb, i, k + 1)
            for c in "ACGT":
                put(sb[:k] + c + sb[k + 1:], i, k + 1)
    for i, sb in enumerate(whitelist):
        if "N" not in sb:
            table[sb] = (i, 0)
    return table


def main():
    rnd = random.Random(1)
    L = 8  # short barcodes so that collisions / ambiguities actually happen
    wl = sorted({"".join(rnd.choice("ACGT") for _ in range(L)) for _ in range(3000)})
    table = julia_table(wl)
    index = {b: i for i, b in enumerate(wl)}

    # queries: whitelist, all their neighbours, and random strings (some with N)
    queries = set(wl) | set(table) | {"".join(rnd.choice("ACGTN") for _ in range(L)) for _ in range(20000)}
    bad = 0
    for q in queries:
        expected = table.get(q, (-1, 0))[0]
        expected = None if expected < 0 else expected  # -1 no bead, -2 ambiguous: both unused
        got = find_bead(q, index)
        if got != expected:
            bad += 1
            if bad <= 5:
                print("MISMATCH", q, "julia", expected, "ours", got)
    print(f"bead matching: {len(queries) - bad}/{len(queries)} queries identical to the Julia table")

    # UMI encoding: Julia dot(4^i, (codeunits>>1)&3)
    ok_umi = True
    for _ in range(2000):
        u = "".join(rnd.choice("ACGT") for _ in range(12))
        julia = sum((4 ** i) * ((ord(c) >> 1) & 3) for i, c in enumerate(u))
        ok_umi &= umi_to_number(u) == julia
    print(f"UMI encoding identical to Julia: {ok_umi}")
    return bad == 0 and ok_umi


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
