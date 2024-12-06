for loop_file in /home/carlos/oldies/manuscripts/review/v3/unibind_data/loops_*.tsv
do
awk 'NR>1' $loop_file > temp.loop
out_name=$(basename $loop_file ".tsv")
python /home/carlos/oldies/projects/ctcf_orient/LoopWidth/LoopWidth_piechart.py  --loop temp.loop --res 10000 --output width_svgs/${out_name}.svg
rm temp.loop
done
