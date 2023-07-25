#!/bin/bash
grna_ids=("1139540371" "1139540424" "1139540434" "1139540471" "1139540505" "1139541055" "1139541475" "1139541557" "1139542201" "1139542458" "1139542593" "1139543542" "1139543681")
for grna_id in ${!grna_ids[@]}; do
  echo "element $grna_id is ${grna_ids[$grna_id]}"
  
  url="https://wge.stemcell.sanger.ac.uk/api/crispr_by_id?id=${grna_ids[$grna_id]}&id=${grna_ids[$grna_id]}&species=Grch38"
  echo $url
  response=$(curl -s -o "grna_${grna_ids[$grna_id]}.txt" "$url")
  
   if [[ $? -eq 0 ]]; then
    # File download successful
    echo "File downloaded for gRNA ID: ${grna_id}"
    aws s3 cp "grna_${grna_ids[$grna_id]}.txt" s3://${1}/project/GeneEditingPipeline/resources/VariantDesigner/other/ 
    
    if [[ $? -eq 0 ]]; then
      echo "File uploaded to S3 bucket for gRNA ID: ${grna_ids[$grna_id]}"
      rm grna_${grna_ids[$grna_id]}.txt
      echo "local file deleted"
    else
      echo "File upload failed for gRNA ID: ${grna_id}"
    fi
  else
    echo "File download failed for gRNA ID: ${grna_id}"
  fi
  
done