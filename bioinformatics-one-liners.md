# Bioinformatics-one-liners
This scripts stores commonly used bash command for Bioinformatics.   
Orignated from Ming Tang and add more examples.

### 1. Handle fasta/fastq file
#### Get the sequences length distribution from a fastq file
```bash
zcat file.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}'
```
#### Reverse complement a sequence
```bash
echo 'ATTGCTATGCTNNNT' | rev | tr 'ACTG' 'TGAC'
```
#### Split a multifasta file into single ones
```bash
csplit -z -q -n 4 -f sequence_ sequences.fasta /\>/ {*}  
```
#### Split a multi-FASTA file into individual FASTA files
```bash
awk '/^>/{s=++d".fa"} {print > s}' multi.fasta
```
#### Number of reads in a fastq file
```bash
zcat file.fastq.gz | echo $((`wc -l`/4))
parallel "echo {} && gunzip -c {} | wc -l | awk '{d=\$1; print d/4;}'" ::: *.gz
```
#### Remove lines contain "N" (or anything) in fasta
```bash
awk -v RS=">" '!/N/{printf $0RT}' in.fasta
```
#### [seqtk](https://github.com/lh3/seqtk) collections
```bash
#Support for .gz file

#Filter fasta lines according to reads length
seqtk seq -L 51 in.fa
#drop sequences containing ambiguous bases
seqtk seq -N in.fa
#Convert FASTQ to FASTA
seqtk seq -A in.fq
#Convert multi-line FASTQ/FASTA to 4-line FASTQ/FASTA
seqtk seq -l0 in.fq/in.fa
seqtk seq -l 60 in.fq/in.fa
#Extract sequences in regions contained in file reg.bed
seqtk subseq in.fa reg.bed
#Subsample 10000 read pairs from two large paired FASTQ files (remember to use the same random seed to keep pairing)
seqtk sample -s100 read1.fq 10000
seqtk sample -s100 read2.fq 10000
#get chr, length, #A, #C, #G, #T, #2, #3, #4, #CpG, #tv, #ts, #CpG-ts
seqtk comp [-u] [-r in.bed] in.fa
#identify high- or low-GC regions
seqtk gc in.fa
```
#### Sequence length of every entry in a multifasta file
```bash
awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' in.fasta
```
#### Deinterleaving a FASTQ
```bash
cat file.fq | paste - - - - - - - - | tee >(cut -f1-4 | tr '\t'  
'\n' > out1.fq) | cut -f5-8 | tr '\t' '\n' 
```
#### grep fastq reads containing a pattern but maintain the fastq format
```bash
grep -A 2 -B 1 'AAGTTGATAACGGACTAGCCTTATTTT' in.fastq | sed '/^--$/d' 
# or
zcat reads.fq.gz | paste - - - - | awk -v FS="\t" -v OFS="\n" '$2 ~ "AAGTTGATAACGGACTAGCCTTATTTT" {print $1, $2, $3, $4}' 
```

### 2. Handle bed/bam/bw/vcf.gtf/bedGraph file
*Most manipulations of bed/fasta/fastq can be solved using bedtools/samtools/seqtk/bedtk/deeptools !!!*
#### [bedtk](https://github.com/lh3/bedtk) collections
```bash
#Support for .gz file

#intersection (no sorting needed; bedtools intersect)
bedtk isec [options] <loaded.bed> <streamed.bed>
#subtraction (bedtools subtract)
bedtk sub <loaded.bed> <streamed.bed>
#total region length
bedtk sum <in.bed>
# filter a BED or VCF file
bedtk flt test-anno.bed.gz test-iso.bed.gz
bedtk flt -v test/test-anno.bed.gz test-iso.bed.gz      # non-overlapping lines
bedtk flt -cw100 test-anno.bed.gz test-sub.vcf.gz  # with a window
# compute the breadth of coverage
bedtk cov test-anno.bed.gz test-iso.bed.gz
# sort a BED file
bedtk sort test-iso.bed.gz
bedtk sort -s chr_list.txt test-iso.bed.gz
# merge overlapping records (no sorting needed)
bedtk merge test-anno.bed.gz
```

#### Filter duplicated lines according to some lines (like dplyr::filter in r)
```bash
awk '!a[$1$2$3]++' in.file 
#or
sort -u -k1,2,3 in.file
```
#### bam2bed (use bedtools bamtobed/bedtobam for interchange)
```bash
samtools view file.bam | perl -F'\t' -ane '$strand=($F[1]&16)?"-":"+";$length=1;$tmp=$F[5];$tmp =~ s/(\d+)[MD]/$length+=$1/eg;print "$F[2]\t$F[3]\t".($F[3]+$length)."\t$F[0]\t0\t$strand\n";' > file.bed
```
#### bam2wig (deeptools bamCompare -> bedGraph/bigwig)
```bash
samtools mpileup -BQ0 file.sorted.bam | perl -pe '($c, $start, undef, $depth) = split;if ($c ne $lastC || $start != $lastStart+1) {print "fixedStep chrom=$c start=$start step=1 span=1\n";}$_ = $depth."\n";($lastC, $lastStart) = ($c, $start);' | gzip -c > file.wig.gz
```
#### Handle each chromosome separately and parallely running them on several cores
```bash
samtools view -H yourFile.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P 24 sh -c "samtools mpileup -BQ0 -d 100000 -uf yourGenome.fa -r {} yourFile.bam | bcftools view -vcg - > tmp.{}.vcf"
```
#### To merge the results afterwards
```bash
samtools view -H yourFile.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | perl -ane 'system("cat tmp.$F[0].bcf >> yourFile.vcf");'
```
#### split large file by id/label/column
```bash
awk '{print >> $1; close($1)}' input_file
```
#### split a bed file by chromosome
```bash
awk '{print $0 >> $1".bed"}' example.bed
```
#### sort vcf file with header
```bash
cat my.vcf | awk '$0~"^#" { print $0; next } { print $0 | "sort -k1,1V -k2,2n" }'
```
#### Merge all bed files and add a column for the filename
```bash
awk '{print $0 "\t" FILENAME}' *bed 
```
#### Add or remove chr from the start of each line
```bash
# add chr
sed 's/^/chr/' my.bed
# remove chr
sed 's/^chr//' my.bed
```
#### Cut out columns based on column names in another file
```bash
csvtk cut -t -f $(paste -s -d , list.txt) data.tsv
```
#### Parallelized samtools mpileup 
https://www.biostars.org/p/134331/
```bash
BAM="yourFile.bam"
REF="reference.fasta"
samtools view -H $BAM | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P 24 sh -c "samtools mpileup -BQ0 -d 100000 -uf $REF -r \"{}\" $BAM | bcftools call -cv > \"{}\".vcf"
```
#### Get the promoter regions from a gtf file
Create TSS bed from GTF in one line: 
```bash
zcat gencode.v29lift37.annotation.gtf.gz | awk '$3=="gene" {print $0}' | grep protein_coding | awk -v OFS="\t" '{if ($7=="+") {print $1, $4, $4+1} else {print $1, $5-1, $5}}' > tss.bed
```
### 3. Handle various
#### gnu sed print invisible characters
```bash
cat my_file | sed -n 'l'
cat -A
```
#### exit a dead ssh session
```bash
~.
```
#### copy large files, copy the from_dir directory inside the to_dir directory
```bash
rsync -av from_dir to_dir
## copy every file inside the frm_dir to to_dir
rsync -av from_dir/ to_dir
##re-copy the files avoiding completed ones:
rsync -avhP /from/dir /to/dir
```
#### mkdir using the current date
```bash
mkdir $(date +%F)
```
#### All the folders' size in the current folder (GNU du)
```bash
du -h --max-depth=1
```
#### Pretty output
```bash
fold -w 60
cat file.txt | column -t | less -S
```
#### Pass tab as delimiter
```bash
-t $'\t'
```
#### awk with the first line printed always
```bash
awk ' NR ==1 || ($10 > 1 && $11 > 0 && $18 > 0.001)'  input_file
```
#### Delete blank lines with sed
```bash
sed /^$/d
```
#### Delete the last line
```bash
sed $d
```
#### Select lines from a file based on columns in another file
```bash
awk -F"\t" 'NR==FNR{a[$1$2$3]++;next};a[$1$2$3] > 0' file2 file1 
```
#### Finally learned about the !$ in unix: take the last thing (word) from the previous command
```bash
echo hello, world; echo !$` gives 'world'
```
#### Create a script of the last executed command
```bash
echo "!!" > foo.sh
```
#### Reuse all parameter of the previous command line
```bash 
!*
```
#### find bam in current folder (search recursively) and copy it to a new directory using 5 CPUs  
```bash
find . -name "*bam" | xargs -P5 -I{} rsync -av {} dest_dir
ls -X will group files by extension
```
#### count how many columns in a tsv file
```bash
cat file.tsv | head -1 | tr "\t" "\n" | wc -l  
## from csvkit
csvcut -n -t file.
## emulate csvcut -n -t
less files.tsv | head -1| tr "\t" "\n" | nl
awk -F "\t" 'NR == 1 {print NF}' file.tsv
awk '{print NF; exit}'
```
#### mkdir and cd into that dir shortcut
```bash
mkdir blah && cd $_
```
#### Check if a tsv files have the same number of columns for all rows
```bash
awk '{print NF}' test.tsv | sort -nu | head -n 1
```
#### Convert multiple lines to a single line
This is better than `tr "\n" "\t"` because somtimes I do not want to convert the last newline to tab.
```bash
cat myfile.txt | paste -s 
```
#### Merge multiple files with same header by keeping the header of the first file
```bash
awk 'FNR==1 && NR!=1{next;}{print}' *.csv 
# or
awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' file*.txt
```
#### Insert a field into the first line
```bash
cut -f1-4 F5.hg38.enhancers.expression.usage.matrix | head
CNhs11844	CNhs11251	CNhs11282	CNhs10746
chr10:100006233-100006603	1	0	0
chr10:100008181-100008444	0	0	0
chr10:100014348-100014634	0	0	0
chr10:100020065-100020562	0	0	0
chr10:100043485-100043744	0	0	0
chr10:100114218-100114567	0	0	0
chr10:100148595-100148922	0	0	0
chr10:100182422-100182522	0	0	0
chr10:100184498-100184704	0	0	0

sed '1 s/^/enhancer\t/' F5.hg38.enhancers.expression.usage.matrix | cut -f1-4 | head
enhancer	CNhs11844	CNhs11251	CNhs11282
chr10:100006233-100006603	1	0	0
chr10:100008181-100008444	0	0	0
chr10:100014348-100014634	0	0	0
chr10:100020065-100020562	0	0	0
chr10:100043485-100043744	0	0	0
chr10:100114218-100114567	0	0	0
chr10:100148595-100148922	0	0	0
chr10:100182422-100182522	0	0	0
chr10:100184498-100184704	0	0	0

```
#### Replace a pattern in a specific column
```bash
## column5 
awk '{gsub(pattern,replace,$5)}1' in.file
## http://bioinf.shenwei.me/csvtk/usage/
csvtk replace -f 5 -p pattern -r replacement 
```
#### Move a process to a screen session
https://www.linkedin.com/pulse/move-running-process-screen-bruce-werdschinski/
```bash
1. Suspend: Ctrl+z
2. Resume: bg
3. Disown: disown %1
4. Launch screen
5. Find pid: prep BLAH
6. Reparent process: reptyr
```
#### Count uinque values in a column and put in a new 
```bash
# input
blabla_1 A,B,C,C
blabla_2 A,E,G
blabla_3 R,Q,A,B,C,R,Q

# output
blabla_1 3
blabla_2 3
blabla_3 5

awk '{split(x,C); n=split($2,F,/,/); for(i in F) if(C[F[i]]++) n--; print $1, n}' file
```
#### Reverse one column of a txt file
reverse column 3 and put it to column5
```bash
awk -v OFS="\t" '{"echo "$3 "| rev" | getline $5}{print $0}' 
```
#### Get the full path of a file
```bash
realpath file.txt
readlink -f file.txt 
```
#### Rename a file, bash string manipulation
```bash
rename 's/xx/yy/' *fa
```
#### Run simple command in one-line manner
```bash
ls * | while read i; do *** ;done
```
