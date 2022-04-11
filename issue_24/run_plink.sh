#!/bin/bash
~/.local/share/plinkr/plink_1_9_unix/plink --bfile issue_24 --snp snp_5 --window 0.002 --make-bed --out issue_24_subset ; cat issue_24_subset.bim


# plink --file data --extract mysnps.txt 
