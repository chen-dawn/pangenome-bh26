# pangenome-bh26

Work from the BioHackathon 2026 (Japan) on human pangenomes and the HLA region,
run on the NIG supercomputer BioHackathon node.

- `hla/` CWL workflow: extract the MHC from haplotype assemblies (pgr-tk),
  call HLA/KIR/C4 alleles (Immuannot), fetch every HLA gene across haplotypes
  (pgr-tk) and visualise them (pgr-tk principal bundles, pggb + odgi).
  See `hla/README.md`.
- `data/` provenance notes for the pangenome datasets uploaded to
  `/home/asianhla/data/upload/` on the NIG node (Arab Pangenome Reference,
  Korean pangenome K-PanRef).
- `1000g_ground_truth/` experimental HLA gold-standard labels for validating
  `hla/` pipeline calls: Gourraud et al. 2014 Sanger typing of 1,267 1000
  Genomes samples, and Lai et al. 2024 4-field HPRC labels for 44 samples.
  See `1000g_ground_truth/README.md`.
