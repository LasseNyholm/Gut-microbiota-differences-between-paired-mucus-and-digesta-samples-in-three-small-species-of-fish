#Modified from protocol developed by Ostaizka Aizpurua
#TERMINAL: prepare the sequences for DADA2
#Demultiplexing

module load AdapterRemoval/2.2.2

FASTQ='/path/to/rawdate'
INDEX='/path/to/indexfiles'
OUT='/output/path'

cd ${FASTQ}
find *_[1-2].fastq.gz > list
sed 's/_[1-2].fastq.gz//g' list > list2
uniq list2 > sample_list
rm -f list*
sample_list=$(cat sample_list)
echo ${sample_list[@]}

for a in $sample_list
  do
    AdapterRemoval --file1 ${FASTQ}/"$a"_1.fastq.gz --file2 ${FASTQ}/"$a"_2.fastq.gz --basename ${OUT}/ --minquality 30 --trimns --maxns 5 --trimqualities --threads 40 --gzip --barcode-list ${INDEX}/"$a".txt --barcode-mm-r1 2 --barcode-mm-r2 2
done

cd ${OUT}
rm .*.discarded.gz
rm .*.settings
rm .*.singleton.truncated.gz
rm .*unidentified*

# Renaming files
find .*pair1.truncated.gz > list
sed 's/.//' list > list2
sed 's/.pair1.truncated.gz//' list2 > list3
uniq list3 > sample_list
rm -f list*
sample_list=$(cat sample_list)
echo ${sample_list[@]}

for a in $sample_list
  do
    mv ."$a".pair1.truncated.gz  "$a"_1.fq.gz
    mv ."$a".pair2.truncated.gz  "$a"_2.fq.gz
done


# Merge replicates
find *1.fq.gz > temp
sed 's/_[1-3].1.fq.gz//g' temp > temp2
uniq temp2 > sample_list.txt
rm temp2
sample_list=$(cat sample_list.txt)
for sample in $sample_list
  do
    cat ${sample}_1.1.fq.gz ${sample}_2.1.fq.gz ${sample}_3.1.fq.gz > ${sample}.1.fq.gz
    cat ${sample}_1.2.fq.gz ${sample}_2.2.fq.gz ${sample}_3.2.fq.gz > ${sample}.2.fq.gz
done

mkdir Single_Replicates
for sample in $sample_list
  do
    mv ${sample}_1.1.fq.gz Single_Replicates/
    mv ${sample}_2.1.fq.gz Single_Replicates/
    mv ${sample}_3.1.fq.gz Single_Replicates/
    mv ${sample}_1.2.fq.gz Single_Replicates/
    mv ${sample}_2.2.fq.gz Single_Replicates/
    mv ${sample}_3.2.fq.gz Single_Replicates/
done


# Get samples name
ls *.1.fq.gz | sed "s/\.1\.fq\.gz//g" > samples

#REMOVE PRIMERS
grep -v "_S" samples > samples_tagsteady
wc -l samples
wc -l samples_tagsteady

#samples_tagsteady
for sample in $(cat samples_tagsteady)
do
echo "On sample: $sample"
cutadapt -e 0.15 -g ^CTANGGGNNGCANCAG -G ^GACTACNNGGGTATCTAAT \
--discard-untrimmed \
-o ${sample}_1a_trimmed.fq.gz -p ${sample}_2a_trimmed.fq.gz \
${sample}.1.fq.gz ${sample}.2.fq.gz \
>> cutadapt_primer_trimming_stats.txt 2>&1


echo "On sample: $sample"
cutadapt -e 0.15 -g ^GGACTACNNGGGTATCTAAT -G ^CCTANGGGNNGCANCAG \
--discard-untrimmed \
-o ${sample}_1b_trimmed.fq.gz -p ${sample}_2b_trimmed.fq.gz \
${sample}.1.fq.gz ${sample}.2.fq.gz \
>> cutadapt_primer_trimming_stats.txt 2>&1

#Merge both files
cat ${sample}_1a_trimmed.fq.gz ${sample}_2b_trimmed.fq.gz > ${sample}_1_trimmed.fq.gz
cat ${sample}_2a_trimmed.fq.gz ${sample}_1b_trimmed.fq.gz > ${sample}_2_trimmed.fq.gz

#Remove intermediate files
rm ${sample}_1a_trimmed.fq.gz ${sample}_2b_trimmed.fq.gz ${sample}_2a_trimmed.fq.gz ${sample}_1b_trimmed.fq.gz
done
