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
