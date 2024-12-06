motif_file="/home/carlos/oldies/projects/ctcf_orient/CTCF_hg38.reformat.bed"
for loop_file in /home/carlos/oldies/manuscripts/review/v3/unibind_data/loops_*.tsv
do
awk 'NR>1' $loop_file > temp.loop
out_name=$(basename $loop_file ".tsv")
python /home/carlos/oldies/projects/ctcf_orient/CTCF_orientation/ctcf_orientation.py -l temp.loop -m $motif_file -o orient_svgs/${out_name}.svg
rm temp.loop
done
