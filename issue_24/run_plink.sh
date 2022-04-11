#!/bin/bash
~/.local/share/plinkr/plink_1_9_unix/plink \
  --bfile issue_24 \
  --pheno issue_24.phe \
  --snp snp_5 --window 2 \
  --make-bed --out issue_24_subset
