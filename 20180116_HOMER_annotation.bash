# This code uses HOMER to annotate the delta 20 methylation regions and then analyses the results in R

cd risk

PATH=$PATH:/data/seqtools/homer/bin/
PATH=$PATH:/data/seqtools/weblogo/
PATH=$PATH:/data/seqtools/samtools-1.1/


/data/seqtools/bowtie2-2.2.3/bowtie2 -x mm9 -U $filename > ./sam/"`basename "$filename" .trimmed.fastq`.sam" -p 8;
samtools view -bS ./sam/"`basename "$filename" .trimmed.fastq`.sam" > ./bam/"`basename "$filename" .trimmed.fastq`.bam";
done

find . -name "*.txt" | while read filename; 
do 
annotatePeaks.pl $filename hg19 -annStats "`basename "$filename" .txt`_homer_annStats.txt" > "`basename "$filename" .txt`_homer.txt"
done