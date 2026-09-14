#!/usr/bin/env python3
"""Build the CWL job file (YAML) listing haplotype assemblies from several cohorts.

Naming rules:
  APR      aprNNN.{1,2}.polished.fa.gz            -> sample aprNNN, hap 1/2
  HPRC r2  SAMPLE_{hap1,hap2}_hprc_r2_*.fa         -> hap 1/2
           SAMPLE_{pat,mat}_hprc_r2_*.fa           -> pat=1, mat=2 (HPRC convention)
  JaSaPaGe ksaNNN.hapN.asm.clean.fasta            -> hap N
           NAxxxxx.hifiasm.*.hapN.clean.fasta      -> hap N
  REF      GRCh38 / CHM13 PanSN fasta               -> cohort REF, hap 0
"""
import argparse
import os
import re
import sys

import yaml

RULES = [
    ("APR", re.compile(r"^(apr[-\w]*?)\.(\d)\.polished\.fa(\.gz)?$")),
    ("HPRC_r2", re.compile(r"^([A-Za-z0-9]+)_hap(\d)_hprc_r2.*\.fa(\.gz)?$")),
    ("HPRC_r2", re.compile(r"^([A-Za-z0-9]+)_(pat|mat)_(?:hprc_r2|v1\.0).*\.fa(\.gz)?$")),
    ("HPRC_r2", re.compile(r"^(hg002)v1\.1\.(pat|mat).*\.PanSN\.fa(\.gz)?$")),
    ("JaSaPaGe", re.compile(r"^(ksa\d+)\.hap(\d)\.asm\.clean\.fasta(\.gz)?$")),
    ("JaSaPaGe", re.compile(r"^([A-Za-z0-9]+)\.hifiasm\..*hap(\d)\.clean\.fasta(\.gz)?$")),
    ("REF", re.compile(r"^(GCA_000001405\.15_GRCh38)_no_alt_analysis_set\.PanSN\.fa(?:\.gz)?$")),
    ("REF", re.compile(r"^(chm13)v2\.0_maskedY_rCRS\.fa\.PanSN\.fa(?:\.gz)?$")),
]
HAP = {"pat": "1", "mat": "2"}


def classify(fn):
    for cohort, rx in RULES:
        m = rx.match(fn)
        if m:
            hap = m.group(2) if m.re.groups > 1 and m.group(2) else "0"
            sample = {"GCA_000001405.15_GRCh38": "GRCh38"}.get(m.group(1), m.group(1))
            return cohort, sample, HAP.get(hap, hap)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", action="append", required=True, help="cohort=path (cohort label only used for reporting)")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--reference-chrom", default="GRCh38#0#chr6")
    ap.add_argument("--region", default="GRCh38#0#chr6:28510120-33480577")
    ap.add_argument("--immuannot-dir", required=True)
    ap.add_argument("--immuannot-ref", required=True)
    ap.add_argument("--mhc-dir", help="directory with precomputed <sample>_<hap>.mhc.fa/.mhc.tsv; skips extraction")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--per-cohort-limit", type=int, default=0, help="keep at most N haplotypes per cohort (references always kept)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    files, samples, haps, cohorts = [], [], [], []
    skipped = []
    for spec in a.dir:
        label, path = spec.split("=", 1)
        for fn in sorted(os.listdir(path)):
            if fn.endswith((".fai", ".gzi")):
                continue
            if fn.endswith((".fa.gz", ".fasta.gz")) and os.path.exists(os.path.join(path, fn[:-3])):
                continue  # prefer the uncompressed copy while both exist
            c = classify(fn)
            if not c:
                skipped.append(fn)
                continue
            rule_cohort, sample, hap = c
            if rule_cohort == "REF":
                label_use = "REF"
            else:
                label_use = label
            if a.per_cohort_limit and label_use != "REF" and cohorts.count(label_use) >= a.per_cohort_limit:
                continue
            files.append({"class": "File", "path": os.path.abspath(os.path.join(path, fn))})
            samples.append(sample)
            haps.append(hap)
            cohorts.append(label_use)
    # the same individual can be assembled by two projects (e.g. 1000G JPT samples in
    # HPRC r2 and JaSaPaGe); suffix the sample with the cohort for later duplicates so
    # sample#haplotype labels stay unique
    seen = {}
    for i, (smp, hap, co) in enumerate(zip(samples, haps, cohorts)):
        key = (smp, hap)
        if key in seen and cohorts[seen[key]] != co:
            samples[i] = f"{smp}.{co}"
            print(f"duplicate sample {smp}#{hap}: {co} relabelled {samples[i]}", file=sys.stderr)
        else:
            seen[key] = i
    if a.limit:
        files, samples, haps, cohorts = files[:a.limit], samples[:a.limit], haps[:a.limit], cohorts[:a.limit]
    mhc_fa, mhc_tsv = [], []
    if a.mhc_dir:
        for smp, hap in zip(samples, haps):
            fa = os.path.join(a.mhc_dir, f"{smp}_{hap}.mhc.fa")
            tsv = os.path.join(a.mhc_dir, f"{smp}_{hap}.mhc.tsv")
            if not (os.path.exists(fa) and os.path.exists(tsv)):
                sys.exit(f"missing precomputed MHC files for {smp}#{hap}: {fa}")
            mhc_fa.append({"class": "File", "path": os.path.abspath(fa)})
            mhc_tsv.append({"class": "File", "path": os.path.abspath(tsv)})
    job = {
        "assemblies": files,
        "samples": samples,
        "haplotypes": haps,
        "cohorts": cohorts,
        "reference": {"class": "File", "path": os.path.abspath(a.reference)},
        "reference_chrom": a.reference_chrom,
        "mhc_region": a.region,
        "immuannot_dir": {"class": "Directory", "path": os.path.abspath(a.immuannot_dir)},
        "immuannot_ref": {"class": "Directory", "path": os.path.abspath(a.immuannot_ref)},
    }
    if a.mhc_dir:
        job["mhc_fastas"] = mhc_fa
        job["mhc_tsvs"] = mhc_tsv
    with open(a.out, "w") as fh:
        yaml.safe_dump(job, fh, sort_keys=False)
    print(f"{len(files)} haplotypes written to {a.out}", file=sys.stderr)
    for s in skipped:
        print(f"skipped: {s}", file=sys.stderr)


if __name__ == "__main__":
    main()
