"""
primer_blast_check_spyder.py

Checks a list of primer pairs against multiple species genomes using NCBI BLAST.

I suggest you run this script to flag primer pairs that show high-identity hits (pct_identity > 85%) in non-target species, 
then manually run those flagged pairs through the Primer-BLAST web interface to confirm whether they'd actually form an amplicon.


Input:  CSV with one row for each primer pair and 3 columns: pair_name, forward, reverse

Output: summary CSV + per-primer full BLAST XML results

Output columns:
    pair_name — the name of your primer pair from your input CSV
    species — which species this BLAST query was run against
    forward — your forward primer sequence
    reverse — your reverse primer sequence
    hit_id — the NCBI accession number of the database sequence that matched
    hit_def — a description of that sequence
    hit_length — the total length of that database sequence in base pairs
    hsp_score — a raw alignment score; higher = better match (but not directly interpretable on its own)
    hsp_evalue — statistical significance of the match; lower = more significant (<0.01 is worth flagging, <0.001 is strong)
    hsp_identities — the raw count of matching bases in the alignment
    hsp_align_len — the total length of the alignment in base pairs
    pct_identity — hsp_identities / hsp_align_len × 100; the most intuitive quality measure (>85% concerning)
    query_start / query_end — where the alignment starts and ends on your concatenated primer query sequence
    sbjct_start / sbjct_end — where the alignment starts and ends on the database sequence (the genome)
    spans_junction — True if the alignment overlaps both the forward and reverse primer portions of your concatenated query
    

Configured for running directly in Spyder — edit the USER SETTINGS section below.

Dependencies:
    pip install biopython (note: biopython is not in Spyder's bundled python)
"""

import csv
import io
import os
import time
from pathlib import Path

from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML


# ============================================================
#  USER SETTINGS — edit these before running
# ============================================================

# (note: default working directory will be wherever this script is saved)
INPUT_CSV = r"seaotter_primers_biopy_4.28.26.csv"          # path to your primer CSV file
EMAIL     = r"ellen.dymit@oregonstate.edu"       # your email (required by NCBI)
OUTDIR    = r"blast_results"        # output folder (will be created if needed)
DELAY     = 10.0                    # seconds between BLAST calls (keep >= 3)

# If your CSV columns have different names, update these:
COL_NAME    = "pair_name"
COL_FORWARD = "forward"
COL_REVERSE = "reverse"

# ============================================================
#  SPECIES TO CHECK
#  Add or remove species here as needed.
#  Format: ("label", "NCBI organism search term")
# ============================================================

SPECIES = [
    ("sea_otter",  "Enhydra lutris[Organism]"),
    ("human",      "Homo sapiens[Organism]"),
    ("wolf",       "Canis lupus[Organism]"),
    ("coyote",     "Canis latrans[Organism]"),
    ("brown_bear", "Ursus arctos[Organism]"),
]

# ============================================================
#  BLAST PARAMETERS — no need to change these normally
# ============================================================

BLAST_PROGRAM  = "blastn"
BLAST_DATABASE = "nt"       # nucleotide collection; covers all above genomes
WORD_SIZE      = 7          # short word size needed for primers (~20 bp)
EXPECT         = 10       # permissive e-value to catch near-misses
HITLIST_SIZE   = 10         # hits returned per query

# ============================================================


def load_primers(csv_path):
    """Load primer pairs from CSV. Returns list of (name, fwd, rev) tuples."""
    primers = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Normalise column names
            row = {k.strip().lower(): v.strip() for k, v in row.items()}
            name = row.get(COL_NAME.lower())
            fwd  = row.get(COL_FORWARD.lower())
            rev  = row.get(COL_REVERSE.lower())
            if not all([name, fwd, rev]):
                print(f"  [WARNING] Skipping malformed row: {row}")
                continue
            primers.append((name, fwd.upper(), rev.upper()))
    return primers


def blast_primer_pair(name, fwd, rev, species_name, entrez_query):
    """
    Concatenate fwd + short linker + rev and BLAST as a single query.
    Returns raw XML string.
    """
    query_seq = fwd + "NNNNNN" + rev

    print(f"    BLASTing {name} vs {species_name} ... ", end="", flush=True)
    result_handle = NCBIWWW.qblast(
        program      = BLAST_PROGRAM,
        database     = BLAST_DATABASE,
        sequence     = query_seq,
        entrez_query = entrez_query,
        word_size    = WORD_SIZE,
        expect       = EXPECT,
        hitlist_size = HITLIST_SIZE,
        format_type  = "XML",
    )
    xml_content = result_handle.read()
    result_handle.close()
    print("done.")
    return xml_content


def summarise_record(blast_record, fwd_len):
    """
    Extract key info from a BLAST record.
    Returns a list of dicts (one per HSP).

    spans_junction = True means a hit overlaps both the forward and reverse
    primer regions of the concatenated query — suggesting both primers could
    bind nearby on that genome (i.e. a real amplification risk).
    """
    rows = []
    junction_start = fwd_len + 1
    junction_end   = fwd_len + 6   # the NNNNNN linker

    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            rows.append({
                "hit_id":         alignment.hit_id,
                "hit_def":        alignment.hit_def[:120],
                "hit_length":     alignment.length,
                "hsp_score":      hsp.score,
                "hsp_evalue":     hsp.expect,
                "hsp_identities": hsp.identities,
                "hsp_align_len":  hsp.align_length,
                "pct_identity":   round(hsp.identities / hsp.align_length * 100, 1),
                "query_start":    hsp.query_start,
                "query_end":      hsp.query_end,
                "sbjct_start":    hsp.sbjct_start,
                "sbjct_end":      hsp.sbjct_end,
                "spans_junction": (hsp.query_start <= junction_start
                                   and hsp.query_end >= junction_end),
            })
    return rows


def main():
    Entrez.email = EMAIL

    # Set working directory to the folder containing this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    print(f"Working directory: {script_dir}")

    outdir  = Path(OUTDIR)
    xml_dir = outdir / "xml"
    xml_dir.mkdir(parents=True, exist_ok=True)

    primers = load_primers(INPUT_CSV)
    print(f"\nLoaded {len(primers)} primer pairs from '{INPUT_CSV}'.")
    print(f"Species to check: {[s[0] for s in SPECIES]}")
    print(f"Total BLAST calls to make: {len(primers) * len(SPECIES)}\n")

    summary_path = outdir / "summary.csv"
    summary_fields = [
        "pair_name", "species", "forward", "reverse",
        "hit_id", "hit_def", "hit_length",
        "hsp_score", "hsp_evalue", "hsp_identities",
        "hsp_align_len", "pct_identity",
        "query_start", "query_end", "sbjct_start", "sbjct_end",
        "spans_junction",
    ]

    with open(summary_path, "w", newline="") as summary_f:
        writer = csv.DictWriter(summary_f, fieldnames=summary_fields)
        writer.writeheader()

        total = len(primers) * len(SPECIES)
        count = 0

        for pair_name, fwd, rev in primers:
            for species_name, entrez_query in SPECIES:
                count += 1
                print(f"[{count}/{total}] {pair_name} × {species_name}")

                xml_path = xml_dir / f"{pair_name}__{species_name}.xml"

                # ── Skip if already done (safe to re-run after interruption) ──
                if xml_path.exists():
                    print(f"    Cached result found, skipping BLAST.")
                    with open(xml_path) as xf:
                        blast_records = list(NCBIXML.parse(xf))
                else:
                    try:
                        xml_content = blast_primer_pair(
                            pair_name, fwd, rev, species_name, entrez_query
                        )
                        with open(xml_path, "w") as xf:
                            xf.write(xml_content)
                        blast_records = list(NCBIXML.parse(io.StringIO(xml_content)))
                        time.sleep(DELAY)
                    except Exception as e:
                        print(f"    [ERROR] {e}")
                        continue

                # ── Parse and write to summary ──────────────────────────────
                hits_written = 0
                for blast_record in blast_records:
                    rows = summarise_record(blast_record, len(fwd))
                    for row in rows:
                        row.update({
                            "pair_name": pair_name,
                            "species":   species_name,
                            "forward":   fwd,
                            "reverse":   rev,
                        })
                        writer.writerow(row)
                        hits_written += 1
                summary_f.flush()

                if hits_written == 0:
                    writer.writerow({
                        "pair_name": pair_name, "species": species_name,
                        "forward": fwd, "reverse": rev,
                        "hit_id": "NO_HIT", "hit_def": "No significant hits found",
                        "hit_length": "", "hsp_score": "", "hsp_evalue": "",
                        "hsp_identities": "", "hsp_align_len": "", "pct_identity": "",
                        "query_start": "", "query_end": "", "sbjct_start": "",
                        "sbjct_end": "", "spans_junction": "",
                    })
                    summary_f.flush()
                    print(f"    No hits.")
                else:
                    print(f"    {hits_written} HSP(s) recorded.")

    print(f"\n✓ All done!")
    print(f"  Summary CSV : {summary_path.resolve()}")
    print(f"  Full XML    : {xml_dir.resolve()}")


# ── Run ──────────────────────────────────────────────────────────────────────
main()
