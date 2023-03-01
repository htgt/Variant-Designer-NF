#!/bin/bash

cd /project/intermediate_outputs
primer_file=$1
gene_name=$2

fwd=($(awk '{if (NR!=1) {print $1}}' ${primer_file}))
rev=($(awk '{if (NR!=1) {print $2}}' ${primer_file}))
chrom=($(awk '{if (NR!=1) {print $3}}' ${primer_file}))


rm -f "${gene_name}_internal.bed"
rm -f "${gene_name}_amplicon.bed"

for index in ${!fwd[*]}; do 
   is in 
  perl /home/dispr/dispr.pl \
    --pf "tag:F:${fwd[$index]}" \
    --pr  "tag:R:${rev[$index]}" \
    --ref "/resources/chroms/${chrom}.fa" \
    --bed 'output.bed' \
    --internal-bed 'int.bed' \
    --skip-count \
    --verbose

    cat 'int.bed' >> "${gene_name}_internal.bed"
    cat 'output.bed' >> "${gene_name}_amplicon.bed"
    rm 'int.bed'
    rm 'output.bed'
done
