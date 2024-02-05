bed1="/home/carlos/Desktop/manuscripts/notebooks/gnn/t0_t12_results_0_1.tsv"
bed2="/home/carlos/Desktop/manuscripts/notebooks/gnn/t0_t12_results_1_0.tsv"
echo ########
echo $bed1
echo $bed2
bedtools subtract -a $bed1 -b $bed2 | bedtools sort -i stdin | bedtools merge -i stdin > bed1.bed
bedtools subtract -b $bed1 -a $bed2 | bedtools sort -i stdin | bedtools merge -i stdin > bed2.bed
mkdir -p gnn_res/gnn_0_1_vs_1_0
Rscript lola.R bed1.bed bed2.bed gnn_res/gnn_0_1_vs_1_0
cat gnn_res/gnn_0_1_vs_1_0/extracted_regions/* | sort -k1,1 -k2,2n -u > gnn_res/gnn_0_1_vs_1_0/extracted_regions_merged.bed
rm bed1.bed bed2.bed

bed1="/home/carlos/Desktop/manuscripts/notebooks/gnn/t0_t12_results_1_0.tsv"
bed2="/home/carlos/Desktop/manuscripts/notebooks/gnn/t0_t12_results_0_1.tsv"
echo ########
echo $bed1
echo $bed2
bedtools subtract -a $bed1 -b $bed2 | bedtools sort -i stdin | bedtools merge -i stdin > bed1.bed
bedtools subtract -b $bed1 -a $bed2 | bedtools sort -i stdin | bedtools merge -i stdin > bed2.bed
mkdir -p gnn_res/gnn_1_0_vs_0_1
Rscript lola.R bed1.bed bed2.bed gnn_res/gnn_1_0_vs_0_1
cat gnn_res/gnn_1_0_vs_0_1/extracted_regions/* | sort -k1,1 -k2,2n -u > gnn_res/gnn_1_0_vs_0_1/extracted_regions_merged.bed
rm bed1.bed bed2.bed

bed1="/home/carlos/Desktop/manuscripts/notebooks/gnn/t0_t12_results_0_0.tsv"
bed2="/home/carlos/Desktop/manuscripts/notebooks/gnn/t0_t12_results_1_1.tsv"
echo ########
echo $bed1
echo $bed2
bedtools subtract -a $bed1 -b $bed2 | bedtools sort -i stdin | bedtools merge -i stdin > bed1.bed
bedtools subtract -b $bed1 -a $bed2 | bedtools sort -i stdin | bedtools merge -i stdin > bed2.bed
mkdir -p gnn_res/gnn_0_0_vs_1_1
Rscript lola.R bed1.bed bed2.bed gnn_res/gnn_0_0_vs_1_1
cat gnn_res/gnn_0_0_vs_1_1/extracted_regions/* | sort -k1,1 -k2,2n -u > gnn_res/gnn_0_0_vs_1_1/extracted_regions_merged.bed
rm bed1.bed bed2.bed

bed1="/home/carlos/Desktop/manuscripts/notebooks/gnn/t0_t12_results_1_1.tsv"
bed2="/home/carlos/Desktop/manuscripts/notebooks/gnn/t0_t12_results_0_0.tsv"
echo ########
echo $bed1
echo $bed2
bedtools subtract -a $bed1 -b $bed2 | bedtools sort -i stdin | bedtools merge -i stdin > bed1.bed
bedtools subtract -b $bed1 -a $bed2 | bedtools sort -i stdin | bedtools merge -i stdin > bed2.bed
mkdir -p gnn_res/gnn_1_1_vs_0_0
Rscript lola.R bed1.bed bed2.bed gnn_res/gnn_1_1_vs_0_0
cat gnn_res/gnn_1_1_vs_0_0/extracted_regions/* | sort -k1,1 -k2,2n -u > gnn_res/gnn_1_1_vs_0_0/extracted_regions_merged.bed
rm bed1.bed bed2.bed



