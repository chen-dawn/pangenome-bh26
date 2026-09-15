#!/usr/bin/env python3
"""Tabulate the coding differences Immuannot reports for novel alleles.

For every transcript with a cds_mut attribute, the CDS-relative mutation
offsets are converted into contig coordinates using the CDS rows of the same
transcript (strand-aware), giving one row per substitution:
sample haplotype gene contig pos strand ref_allele closest_allele cds_offset ref_base alt_base aa_change
"""
import gzip, re, sys, csv

def parse_attrs(s):
    return dict(re.findall(r"(\w+) \"?([^\";]*)\"?;", s))

def main(paths, out):
    w = csv.writer(out, delimiter="\t")
    w.writerow(["sample", "haplotype", "gene", "contig", "pos", "strand", "consensus", "closest_allele", "cds_offset", "known_base", "assembly_base", "aa_change"])
    for p in paths:
        label = p.split("/")[-1].replace(".gtf.gz", "")
        sample, hap = label.rsplit("_", 1)
        tx = {}
        for line in gzip.open(p, "rt"):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            a = parse_attrs(f[8])
            tid = a.get("transcript_id")
            if f[2] == "transcript" and "cds_mut" in a:
                tx[tid] = {"contig": f[0], "strand": f[6], "gene": a.get("gene_name"), "consensus": a.get("consensus"), "cds_mut": a["cds_mut"], "cds": []}
            elif f[2] == "CDS" and tid in tx:
                tx[tid]["cds"].append((int(f[3]), int(f[4])))
        for t in tx.values():
            parts = t["cds_mut"].split("|")
            closest, muts = parts[0], parts[1] if len(parts) > 1 else ""
            aas = parts[2].split(":") if len(parts) > 2 else []
            # muts like ":27*cg:339*at:16*ct" = relative offsets, base known->assembly
            offs, pos = [], 0
            for m in [x for x in muts.split(":") if x]:
                mm = re.match(r"(\d+)\*([acgt-]+)", m)
                if not mm:
                    continue
                # cs-style: N matching bases, then the substitution at the next CDS position
                pos += int(mm.group(1)) + 1; offs.append((pos, mm.group(2)))
            cds = sorted(t["cds"])
            if t["strand"] == "-":
                cds = cds[::-1]
            for i, (off, bases) in enumerate(offs):
                # walk the CDS blocks to find the genomic coordinate of CDS offset `off` (1-based)
                rem, gpos = off, None
                for s, e in cds:
                    ln = e - s + 1
                    if rem <= ln:
                        gpos = s + rem - 1 if t["strand"] == "+" else e - rem + 1
                        break
                    rem -= ln
                w.writerow([sample, hap, t["gene"], t["contig"], gpos, t["strand"], t["consensus"], closest, off, bases[0].upper(), bases[1:].upper(), aas[i] if i < len(aas) else ""])

if __name__ == "__main__":
    main(sys.argv[1:], sys.stdout)
