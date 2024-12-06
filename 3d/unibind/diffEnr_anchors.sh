for t1 in 0
#for t1 in 0 12 30
do
  for t2 in 12
  #for t2 in 0 12 30 60
  do
    if [ "$t1" != "$t2" ]; then
      bed1="unibind_data/anchors_specific_anchors1.tsv"
      bed2="unibind_data/anchors_specific_anchors2.tsv"


      echo ########
      echo $bed1
      echo $bed2
      bedtools subtract -a $bed1 -b $bed2 | bedtools sort -i stdin | bedtools merge -i stdin > bed1.bed
      bedtools subtract -b $bed1 -a $bed2 | bedtools sort -i stdin | bedtools merge -i stdin > bed2.bed
      mkdir -p anchor_res/t${t1}_t${t2}_1_2
      Rscript lola.R bed1.bed bed2.bed anchor_res/t${t1}_t${t2}_1_2
      cat anchor_res/t${t1}_t${t2}_1_2/extracted_regions/* | sort -k1,1 -k2,2n -u > anchor_res/t${t1}_t${t2}_1_2/extracted_regions_merged.bed
      rm bed1.bed bed2.bed
      

      bed1="unibind_data/anchors_specific_anchors2.tsv"
      bed2="unibind_data/anchors_specific_anchors1.tsv"

      echo ########
      echo $bed1
      echo $bed2
      bedtools subtract -a $bed1 -b $bed2 | bedtools sort -i stdin | bedtools merge -i stdin > bed1.bed
      bedtools subtract -b $bed1 -a $bed2 | bedtools sort -i stdin | bedtools merge -i stdin > bed2.bed
      mkdir -p anchor_res/t${t1}_t${t2}_2_1
      Rscript lola.R bed1.bed bed2.bed anchor_res/t${t1}_t${t2}_2_1
      rm bed1.bed bed2.bed
      cat anchor_res/t${t1}_t${t2}_2_1/extracted_regions/* | sort -k1,1 -k2,2n -u > anchor_res/t${t1}_t${t2}_2_1/extracted_regions_merged.bed
    fi
  done
done