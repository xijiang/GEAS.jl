#!/usr/bin/env bash
echo Remove duplicates, including the first one
plink --bfile for_CMSEdit \
      --chr-set 30 \
      --list-duplicate-vars \
      --out t
plink --bfile for_CMSEdit \
      --chr-set 30 \
      --chr 1-29 \
      --exclude t.dupvar \
      --recode vcf-iid bgz \
      --out ../../run/a
