project=MIR184

ls ~/MarcosCosta/projects/RNAseq/$project/ > samples.tmp

if [[ off = on ]]; then 

	for sample in $(cat samples.array.tmp); do 
		echo $sample next
		echo "Name,$sample" > $sample.csv
		awk 'NR > 1 { print $1","$5}' ~/MarcosCosta/projects/RNAseq/$project/$sample/star_salmon/$sample/quant.genes.sf >> $sample.csv
		head $sample.csv

	done

	else
	echo "prep is off"
fi 

first=$(ls *.csv | head -1)
second=$(ls *.csv | head -2 | tail -1)
rest=$(ls *.csv | tail -n+3)
echo $first
echo $second 
join -t',' $first $second > counts_data_$project.tmp
head -5 counts_data_$project.tmp

for csv in $rest; do
join -t',' counts_data_$project.tmp $csv > joined.tmp
mv joined.tmp counts_data_$project.tmp 
done
rm *.csv
rm *.tmp
mv counts_data_$project.tmp counts_data_$project.csv

awk -F ',' 'NR == 1 {print; next} {printf "%s", $1; for(i=2; i<=NF; i++) printf ",%d", int($i); printf "\n"}' counts_data_$project.csv > interger_counts_data_$project.csv
