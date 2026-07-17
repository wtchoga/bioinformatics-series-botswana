#!/usr/bin/env python3
"""
Practical 2 | Dataset A - build the real HIV-1 teaching alignment.  [OPTIONAL]

Not run by default. `run_all.sh` uses Datasets B (simulated) and C (H1N1).
Build this one only if you want the NEGATIVE CONTROL - and it is worth it.

WHY IT EXISTS
-------------
Dataset A is real, unlinked, cross-sectional HIV-1 gag spanning 20 years. It
looks like exactly the sort of data you would date. It has almost no temporal
signal (R^2 ~ 0.06), and TreeTime cheerfully dates its common ancestor to
around 1793 - centuries before HIV-1 entered humans - with a ~200-year
confidence interval and no warning whatsoever.

The reason is structural, not a bug: two random subtype-C patients differ by
~0.145 substitutions/site, while 20 years of clock evolution adds only ~0.02.
Ancient divergence swamps the date signal. Subsetting to a single subtype does
not fix it. Nothing in this workshop teaches scepticism faster.

Source : GAG_sdb_n417.fasta - 417 REAL HIV-1 gag sequences (1503 bp, aligned),
         6 African countries, 2005-2025, plus the HXB2 reference.
Output : data/HIV_alignment.fasta  80 sequences, ANONYMISED
         data/reference.fasta      HXB2 (K03455) gag
         data/dates.tsv            TreeTime format, YEAR precision only
         data/metadata.csv         country + year are REAL; province,
                                   risk_group and subtype are SIMULATED
                                   (see data/README_DATA.md)

Sample IDs are re-issued as <ISO3>_WS###_<year> so no original patient or study
identifier is carried into workshop materials.

THREE teaching artefacts are planted deliberately - listed at the end of this
script's output and in notes/INSTRUCTOR_KEY.md. Don't show students the key.

Deterministic: fixed seed, same dataset every run.
Run:  python3 scripts/00a_prepare_hiv.py
Then: bash run_all.sh sim h1n1 real
"""
import random, re, csv, os, sys
from collections import defaultdict, Counter

random.seed(20260716)

SRC = "/Users/wchoga/Downloads/GAG_sdb_n417.fasta"
OUT = "data"

ISO3 = {"BW": "BWA", "KE": "KEN", "RW": "RWA", "UG": "UGA", "ZA": "ZAF", "ZM": "ZMB"}
COUNTRY_NAME = {"BWA": "Botswana", "KEN": "Kenya", "RWA": "Rwanda",
                "UGA": "Uganda", "ZAF": "South Africa", "ZMB": "Zambia"}
# ---- SIMULATED metadata (see data/README_DATA.md) --------------------------
# These columns are FABRICATED so Exercise 13 has something to colour by.
# They are NOT attributes of these samples. Never present them as real.
PROVINCES = {
    "BWA": ["South-East", "Central", "North-East"],
    "KEN": ["Nairobi", "Nyanza", "Coast"],
    "RWA": ["Kigali", "Southern", "Eastern"],
    "UGA": ["Central", "Western", "Northern"],
    "ZAF": ["KwaZulu-Natal", "Gauteng", "Western Cape"],
    "ZMB": ["Lusaka", "Copperbelt", "Southern"],
}
RISK = ["Heterosexual", "MSM", "PWID", "Unknown"]
RISK_W = [0.72, 0.14, 0.06, 0.08]
SUBTYPE = {  # regionally plausible frequencies - NOT called from the sequence
    "BWA": (["C"], [1.0]),
    "ZAF": (["C"], [1.0]),
    "ZMB": (["C", "G"], [0.92, 0.08]),
    "UGA": (["A1", "D", "C"], [0.55, 0.35, 0.10]),
    "KEN": (["A1", "D", "C"], [0.65, 0.25, 0.10]),
    "RWA": (["A1", "C"], [0.80, 0.20]),
}


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:], []
        else:
            buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def write_fasta(path, records, width=60):
    with open(path, "w") as fh:
        for k, v in records:
            fh.write(f">{k}\n")
            for i in range(0, len(v), width):
                fh.write(v[i:i + width] + "\n")


if not os.path.exists(SRC):
    sys.exit(f"Source not found: {SRC}\n"
             "Dataset A is optional - run_all.sh works without it.")

raw = read_fasta(SRC)

# ---- reference: HXB2 -------------------------------------------------------
hxb2 = next(v for k, v in raw.items() if "K03455" in k)

# ---- candidate pool: drop short/empty sequences -----------------------------
pool = []
for k, v in raw.items():
    if "K03455" in k:
        continue
    V = v.upper()
    ung = len(V) - V.count("-")
    m = re.search(r"_(\d{4})$", k)
    if not m:
        continue
    yr, cc = int(m.group(1)), k.split("_")[0]
    if ung >= 1400 and 2005 <= yr <= 2025 and cc in ISO3:
        pool.append({"orig": k, "seq": V, "cc": ISO3[cc], "year": yr})

# ---- stratified subsample: spread across years, balanced across countries ---
by_country = defaultdict(list)
for r in pool:
    by_country[r["cc"]].append(r)

TARGET = 77          # + 3 planted teaching cases = 80
per = TARGET // len(by_country)
picked = []
for cc, rows in sorted(by_country.items()):
    byyr = defaultdict(list)
    for r in rows:
        byyr[r["year"]].append(r)
    years = sorted(byyr)
    take, yi = [], 0
    while len(take) < per and years:
        y = years[yi % len(years)]
        if byyr[y]:
            take.append(byyr[y].pop(random.randrange(len(byyr[y]))))
        else:
            years.remove(y)
            if not years:
                break
            yi -= 1
        yi += 1
    picked.extend(take)

# top up with the oldest available seqs (widens the temporal span)
rest = sorted([r for r in pool if r not in picked], key=lambda r: r["year"])
for r in rest:
    if len(picked) >= TARGET:
        break
    picked.append(r)

picked.sort(key=lambda r: (r["cc"], r["year"]))

# ---- anonymise -------------------------------------------------------------
records = []
for i, r in enumerate(picked, start=1):
    r["id"] = f"{r['cc']}_WS{i:03d}_{r['year']}"
    records.append((r["id"], r["seq"]))

# ---- PLANTED TEACHING CASES (see notes/INSTRUCTOR_KEY.md) ------------------
planted = []

# 1) High-N sequence: take an UNUSED sequence and mask a 300 bp window.
#    It must come from the unused pool, NOT be a copy of a sequence already in
#    the dataset - otherwise it is also a duplicate, and the "bad quality"
#    lesson gets tangled up with the "duplicate" lesson (it would cluster with
#    its own twin at zero distance).
donor = min([r for r in pool if "id" not in r], key=lambda r: r["year"])
s = list(donor["seq"])
s[400:700] = ["N"] * 300
qc_id = f"{donor['cc']}_WS078_{donor['year']}"
records.append((qc_id, "".join(s)))
donor["id"] = qc_id                      # mark as used
picked.append({"orig": "PLANTED-highN", "seq": "".join(s), "cc": donor["cc"],
               "year": donor["year"], "id": qc_id})
planted.append((qc_id, "high-N (300 Ns masked in) -> Exercise 1 / QC removal"))

# 2) Date outlier: a genuinely OLD sequence relabelled with a RECENT date.
#    Also from the unused pool, so it is not an exact duplicate of anything.
#    NOTE: this one is deliberately (and instructively) HARD to detect - see
#    the key. Dataset A has no temporal signal, so root-to-tip outlier
#    detection quietly stops working. That is the lesson of Exercise 6.
unused = [r for r in pool if "id" not in r]
old = min(unused, key=lambda r: r["year"])
out_id = f"{old['cc']}_WS079_2024"          # true year is old['year']
records.append((out_id, old["seq"]))
old["id"] = out_id
picked.append({"orig": "PLANTED-dateoutlier", "seq": old["seq"], "cc": old["cc"],
               "year": 2024, "id": out_id})
planted.append((out_id, f"date outlier: sequence truly from {old['year']}, "
                        f"metadata says 2024 -> Exercise 6 (hard to detect ON PURPOSE)"))

# 3) Exact duplicate of an existing sample (same patient, resequenced)
dup_src = picked[30]
dup_id = f"{dup_src['cc']}_WS080_{dup_src['year']}"
records.append((dup_id, dup_src["seq"]))
picked.append({"orig": "PLANTED-duplicate", "seq": dup_src["seq"],
               "cc": dup_src["cc"], "year": dup_src["year"], "id": dup_id})
planted.append((dup_id, f"exact duplicate of {dup_src['id']} -> Exercise 15 removal"))

# ---- metadata + dates ------------------------------------------------------
meta = []
for r in picked:
    cc = r["cc"]
    sts, stw = SUBTYPE[cc]
    meta.append({
        "sample_id": r["id"],
        "country": COUNTRY_NAME[cc],       # REAL
        "country_code": cc,                # REAL
        "year": r["year"],                 # REAL, year-precision only
        "province": random.choice(PROVINCES[cc]),          # SIMULATED
        "risk_group": random.choices(RISK, RISK_W)[0],     # SIMULATED
        "subtype": random.choices(sts, stw)[0],            # SIMULATED
        "gene": "gag",
        "length_bp": 1503,
    })

os.makedirs(OUT, exist_ok=True)
write_fasta(f"{OUT}/HIV_alignment.fasta", records)
write_fasta(f"{OUT}/reference.fasta", [("HXB2_K03455_gag", hxb2)])

with open(f"{OUT}/metadata.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(meta[0].keys()))
    w.writeheader()
    w.writerows(meta)

# dates.tsv - YEAR precision. The source headers carry no month/day, so we do
# not invent one. Year-only dates mean "January 1st", injecting up to a year of
# error per tip. Worth remembering when you look at the temporal signal.
with open(f"{OUT}/dates.tsv", "w") as fh:
    fh.write("name\tdate\n")
    for m in meta:
        fh.write(f"{m['sample_id']}\t{m['year']}\n")

print(f"sequences written : {len(records)}")
print(f"alignment width   : {len(records[0][1])}")
print(f"year span         : {min(m['year'] for m in meta)} - {max(m['year'] for m in meta)}")
print(f"by country        : {dict(Counter(m['country_code'] for m in meta))}")
print("\nPLANTED teaching cases (INSTRUCTOR ONLY):")
for pid, why in planted:
    print(f"  {pid:22s} {why}")
print("\nReminder: province / risk_group / subtype in metadata.csv are SIMULATED.")
print("See data/README_DATA.md. Do not analyse them as real data.")
print("\nNext:  bash run_all.sh sim h1n1 real")
