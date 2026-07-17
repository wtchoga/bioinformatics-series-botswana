#!/usr/bin/env python3
"""
Practical 2 | Dataset C - prepare the 2009 pandemic H1N1 data.

The supplied file is called "H1N1pdm_2009.fasta" but it is NOT a FASTA file.
It is a NEXUS file, written with classic-Mac CR (\\r) line endings, which is
why grep and most tools report "0 sequences" in it. This is a real-world data
problem, so we fix it explicitly rather than hiding it.

  Input : New Folder With Items 3/H1N1pdm_2009.fasta   (NEXUS, CR endings)
  Output: data/h1n1/H1N1_alignment.fasta   50 whole genomes, 13109 bp
          data/h1n1/H1N1_dates.tsv         decimal sampling dates
          data/h1n1/H1N1_metadata.csv      location + region

Taxon names look like:   'A/Beijing/01/2009_135_2009.37'
                          <-- strain -->  ^dayofyear ^decimal date

Run:  python3 scripts/00b_prepare_h1n1.py
"""
import os, re, csv, sys
from collections import Counter

SRC = "/Users/wchoga/Downloads/New Folder With Items 3/H1N1pdm_2009.fasta"
OUT = "data/h1n1"
os.makedirs(OUT, exist_ok=True)

if not os.path.exists(SRC):
    sys.exit(f"Source not found: {SRC}")

# --- 1. read and normalise line endings -------------------------------------
with open(SRC, "rb") as fh:
    raw = fh.read().decode("ascii", errors="replace")
text = raw.replace("\r\n", "\n").replace("\r", "\n")

ntax = int(re.search(r"ntax\s*=\s*(\d+)", text, re.I).group(1))
nchar = int(re.search(r"nchar\s*=\s*(\d+)", text, re.I).group(1))
print(f"NEXUS header claims: ntax={ntax}  nchar={nchar}")

# --- 2. pull the matrix -----------------------------------------------------
body = text.split("Matrix", 1)[1].split("\n;", 1)[0]

seqs = {}
for line in body.split("\n"):
    line = line.strip()
    if not line or line.startswith("["):        # [ ... ] are ruler comments
        continue
    parts = line.split(None, 1)
    if len(parts) != 2:
        continue
    name, seq = parts
    name = name.strip("'\"")
    seq = re.sub(r"\s+", "", seq).upper()
    if not re.fullmatch(r"[ACGTUNRYSWKMBDHV\-\?]+", seq):
        continue
    seqs[name] = seqs.get(name, "") + seq        # interleaved-safe

print(f"parsed sequences   : {len(seqs)}")
widths = {len(s) for s in seqs.values()}
print(f"sequence widths    : {widths}")
assert len(seqs) == ntax, f"expected {ntax} taxa, parsed {len(seqs)}"
assert widths == {nchar}, f"expected width {nchar}, got {widths}"

# --- 3. parse dates out of the taxon names ----------------------------------
rows = []
for name, seq in seqs.items():
    m = re.match(r"^(.*)_(\d+)_(\d{4}\.\d+)$", name)
    if not m:
        raise ValueError(f"cannot parse date from taxon name: {name!r}")
    strain, doy, date = m.group(1), int(m.group(2)), float(m.group(3))
    loc = strain.split("/")[1].replace("_", " ")   # A/<location>/<n>/<year>
    rows.append({"strain": strain, "doy": doy, "date": date,
                 "location": loc, "seq": seq})
rows.sort(key=lambda r: r["date"])

# --- 4. tidy tip labels: short AND date-bearing ------------------------------
# Keep the date IN the label - TempEst and TreeTime can then read it directly.
clean = lambda s: re.sub(r"[^A-Za-z0-9]+", "-", s).strip("-")
used = set()
for r in rows:
    label, i = f"{clean(r['strain'])}|{r['date']:.3f}", 2
    while label in used:
        label = f"{clean(r['strain'])}-{i}|{r['date']:.3f}"; i += 1
    used.add(label)
    r["id"] = label

# --- 5. write outputs -------------------------------------------------------
with open(f"{OUT}/H1N1_alignment.fasta", "w") as fh:
    for r in rows:
        fh.write(f">{r['id']}\n")
        for i in range(0, len(r["seq"]), 60):
            fh.write(r["seq"][i:i + 60] + "\n")

with open(f"{OUT}/H1N1_dates.tsv", "w") as fh:
    fh.write("name\tdate\n")
    for r in rows:
        fh.write(f"{r['id']}\t{r['date']:.3f}\n")

# continent grouping from the strain location, for colouring only
REGION = {
    "Beijing": "Asia", "Fujian": "Asia", "Guangdong": "Asia", "Hong Kong": "Asia",
    "Hyogo": "Asia", "Japan": "Asia", "Osaka": "Asia", "Shanghai": "Asia",
    "Zhejiang": "Asia", "Nagano": "Asia", "Nagasaki": "Asia", "Kobe": "Asia",
    "Thailand": "Asia", "Singapore": "Asia", "Taiwan": "Asia", "Korea": "Asia",
    "Israel": "Asia", "Narita": "Asia", "Osaka-C": "Asia", "Shandong": "Asia",
    "Sichuan": "Asia", "Nonthaburi": "Asia",
    "California": "North America", "Texas": "North America",
    "New York": "North America", "Canada-NS": "North America",
    "Canada-ON": "North America", "Mexico": "North America",
    "Wisconsin": "North America", "Illinois": "North America",
    "Ohio": "North America", "Massachusetts": "North America",
    "Indiana": "North America", "Michigan": "North America",
    "England": "Europe", "Italy": "Europe", "Spain": "Europe",
    "Germany": "Europe", "Netherlands": "Europe", "Denmark": "Europe",
    "Finland": "Europe", "Russia": "Europe", "Moscow": "Europe",
    "Norway": "Europe", "Paris": "Europe",
    "Auckland": "Oceania", "Australia": "Oceania", "New Zealand": "Oceania",
}
with open(f"{OUT}/H1N1_metadata.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["sample_id", "strain", "location", "region", "date",
                "day_of_year", "year"])
    for r in rows:
        w.writerow([r["id"], r["strain"], r["location"],
                    REGION.get(r["location"], "Other"),
                    f"{r['date']:.3f}", r["doy"], int(r["date"])])

span_days = (rows[-1]["date"] - rows[0]["date"]) * 365
print(f"\ndate range         : {rows[0]['date']:.3f} - {rows[-1]['date']:.3f}"
      f"  ({span_days:.0f} days)")
print(f"locations          : {len(set(r['location'] for r in rows))}")
print(f"regions            : {dict(Counter(REGION.get(r['location'],'Other') for r in rows))}")
print("\nNOTE: the sampling window is only ~6 weeks. Keep that in mind when you")
print("      interpret the root-to-tip R^2 in Exercise 6 - a short window")
print("      flattens R^2 even when the molecular clock itself is fine.")
print(f"\nWrote: {OUT}/H1N1_alignment.fasta")
print(f"Wrote: {OUT}/H1N1_dates.tsv")
print(f"Wrote: {OUT}/H1N1_metadata.csv")
