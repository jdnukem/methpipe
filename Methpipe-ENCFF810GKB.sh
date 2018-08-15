# Methpipe

################################################################################
# building and installing
################################################################################
cd /home/segil_lab/tools/
git clone --recursive https://github.com/smithlabcode/methpipe.git
make all
    path: /home/segil_lab/tools/methpipe/
sudo make install

conda install -c bioconda pyfaidx
cd /media/segil_lab/SegilRaid/staging/Duc/mm10_bis/mm10
faidx -x mm10.fa
    # this splits mm10.fa into individual *.fa files

cd /home/segil_lab/tools/
git clone https://github.com/smithlabcode/walt.git
cd /home/segil_lab/tools/walt
make all
sudo make install
    #note: commented out checks for .fastq, .fa file extensions in src/walt.cpp


/home/segil_lab/tools/walt/bin/makedb \
-c /media/segil_lab/SegilRaid/staging/Duc/mm10_bis/mm10/ \
-o /media/segil_lab/SegilRaid/staging/Duc/mm10_bis/mm10/mm10.fa.dbindex
    outputs:
      mm10.fa.dbindex       # use this index file for walt
      mm10.fa.dbindex_CT00
      mm10.fa.dbindex_CT01
      mm10.fa.dbindex_GA10
      mm10.fa.dbindex_GA11

################################################################################
# ENCFF810GKB.fastq.gz
################################################################################
file="/media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep1/ENCFF810GKB.fastq.gz"
file_name="${file##*/}"
trimmed_file="${file%%.*}_trimmed.fq.gz"
split_file1="${trimmed_file%%.*}.split1.fq.gz"
split_file2="${trimmed_file%%.*}.split2.fq.gz"
mapped_file1="${split_file1%%.fq.gz}.mapped.sam"
mapped_file2="${split_file2%%.fq.gz}.mapped.sam"
merged_file="${mapped_file1%%.*}.splits.merged.mr"
sorted_file="${merged_file}.sorted_start"
dremove_stats="${sorted_bam}.dremove.stat.txt"
dremove_output="${sorted_bam}.dremove"

trim-galore \
--illumina \
ENCFF110EYU.fastq.gz
trimmed_file="${file%%.*}_trimmed.fq.gz"

# file: rep2/ENCFF110EYU_trimmed.fq.gz
#split files
split_file1="${trimmed_file%%.*}.split1.fq.gz"
split_file2="${trimmed_file%%.*}.split2.fq.gz"

lines=$(pigz -cd -p 4 "$file" | wc -l)
half_lines=$(echo "$lines / 2" | bc)
rm_lines=$(echo "$lines - $half_lines" | bc)

pigz -cd -p 4 "$file" | head -"${half_lines}" | \
pigz -p 4 > "$split_file1"

pigz -cd -p 4 "$file_name" | tail -"${rm_lines}" | \
pigz -p 4 > "$split_file2"

# map split files
/home/segil_lab/tools/walt/bin/walt \
-i /media/segil_lab/SegilRaid/staging/Duc/mm10_bis/mm10/mm10.fa.dbindex \
-r <(zcat "$split_file1") \
-N 5000000 \
-u -a \
-o "$mapped_file1"

walt -i <index> -r <reads> -o <output file> [options]
      # head ENCFF810GKB_trimmed.split1.mapped.sam.mapstats
      # [TOTAL NUMBER OF READS: 191369745]
      # [UNIQUELY MAPPED READS: 147336200 (76.99%)]
      # [AMBIGUOUS MAPPED READS: 18077699 (9.45%)]
      # [UNMAPPED READS: 25955846 (13.56%)]

/home/segil_lab/tools/walt/bin/walt \
-i /media/segil_lab/SegilRaid/staging/Duc/mm10_bis/mm10/mm10.fa.dbindex \
-r <(zcat "$split_file2") \
-N 5000000 \
-u -a \
-o "$mapped_file2"
      # head ENCFF810GKB_trimmed.split2.mapped.sam.mapstats
      # [TOTAL NUMBER OF READS: 191369745]
      # [UNIQUELY MAPPED READS: 148038177 (77.36%)]
      # [AMBIGUOUS MAPPED READS: 18346652 (9.59%)]
      # [UNMAPPED READS: 24984916 (13.06%)]

samtools sort -@ 8 ENCFF810GKB_trimmed.split1.mapped.sam -o ENCFF810GKB_trimmed.split1.mapped.sort.sam
samtools sort -@ 8 ENCFF810GKB_trimmed.split2.mapped.sam -o ENCFF810GKB_trimmed.split2.mapped.sort.sam

samtools merge -@ 8 ENCFF810GKB_trimmed.splits.merged.sam \
ENCFF810GKB_trimmed.split1.mapped.sort.sam \
ENCFF810GKB_trimmed.split2.mapped.sort.sam

java -Xmx4G -jar /home/segil_lab/miniconda3/envs/bds_atac/share/picard-1.126-4/picard.jar MarkDuplicates \
INPUT=ENCFF810GKB_trimmed.splits.merged.sam \
OUTPUT=ENCFF810GKB_trimmed.splits.merged.dupmark.sam \
METRICS_FILE=ENCFF810GKB_trimmed.splits.merged.dupmark.qc \
VALIDATION_STRINGENCY=LENIENT \
ASSUME_SORTED=true REMOVE_DUPLICATES=false
      # tail ENCFF810GKB_trimmed.splits.merged.dupmark.qc
      # ## htsjdk.samtools.metrics.StringHeader
      # # picard.sam.markduplicates.MarkDuplicates INPUT=[ENCFF810GKB_trimmed.splits.merged.sam] OUTPUT=ENCFF810GKB_trimmed.splits.merged.dupmark.sam METRICS_FILE=ENCFF810GKB_trimmed.splits.merged.dupmark.qc REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT    MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES READ_NAME_REGEX=[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).* OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 VERBOSITY=INFO QUIET=false COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false
      # ## htsjdk.samtools.metrics.StringHeader
      # # Started on: Tue Jul 31 10:51:50 PDT 2018

      # ## METRICS CLASS	picard.sam.DuplicationMetrics
      # LIBRARY	        UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE
      # Unknown Library	295374377	              0	                  50940762	      37168144	                0	                    0	                            0.125834

samtools view \
-@ 8 \
-F 1804 \
-f 2 \
-b ENCFF810GKB_trimmed.splits.merged.dupmark.sam \
> ENCFF810GKB_trimmed.splits.merged.nodup.sam

bsrate -c mm10 \
-o ENCFF810GKB_trimmed.splits.merged.sorted.nodup.sam.bsrate \
ENCFF810GKB_trimmed.splits.merged.sorted.nodup.sam

===========
samtools view -bS ENCFF810GKB_trimmed.splits.merged.sam \
> ENCFF810GKB_trimmed.splits.merged.bam

STAR --runThreadN 7 \
--limitBAMsortRAM 26843545600 \
--limitIObufferSize 536870912 \
--genomeLoad LoadAndKeep \
--runMode inputAlignmentsFromBAM \
--bamRemoveDuplicatesType UniqueIdentical \
--outBAMsortingThreadN 7 \
--inputBAMfile ENCFF810GKB_trimmed.splits.merged.chr1.bam \
--outFileNamePrefix "dupe_"

samtools rmdup ENCFF810GKB_trimmed.splits.merged.bam \
ENCFF810GKB_trimmed.splits.merged.samtools_dedup.bam



==================
samtools index ENCFF810GKB_trimmed.splits.merged.bam
samtools view -@ 8 -b ENCFF810GKB_trimmed.splits.merged.bam chr1 > \
ENCFF810GKB_trimmed.splits.merged.chr1.bam

STAR --runThreadN 7 \
--limitBAMsortRAM 26843545600 \
--limitIObufferSize 536870912 \
--genomeLoad LoadAndKeep \
--runMode inputAlignmentsFromBAM \
--bamRemoveDuplicatesType UniqueIdentical \
--outBAMsortingThreadN 7 \
--inputBAMfile ENCFF810GKB_trimmed.splits.merged.chr1.bam \
--outFileNamePrefix "dupe_"

STAR --runThreadN 7 \
--runMode inputAlignmentsFromBAM \
--bamRemoveDuplicatesType UniqueIdentical \
--inputBAMfile ENCFF810GKB_trimmed.splits.merged.chr1.bam \
--outFileNamePrefix "dupe_"






LC_ALL=C sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 -o ENCFF810GKB_trimmed.split1.mapped.mr.sorted_start ENCFF810GKB_trimmed.split1.mapped.mr
/home/segil_lab/tools/methpipe-3.4.3/bin/duplicate-remover -S ENCFF810GKB_trimmed.split1.mapped.mr.dremove_stat.txt -o ENCFF810GKB_trimmed.split1.mapped.mr.dremove /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep1/ENCFF810GKB_trimmed.split1.mapped.mr.sorted_start

LC_ALL=C sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 -o ENCFF810GKB_trimmed.split2.mapped.mr.sorted_start ENCFF810GKB_trimmed.split2.mapped.mr
/home/segil_lab/tools/methpipe-3.4.3/bin/duplicate-remover -S ENCFF810GKB_trimmed.split2.mapped.mr.dremove_stat.txt -o ENCFF810GKB_trimmed.split2.mapped.mr.dremove ENCFF810GKB_trimmed.split2.mapped.mr.sorted_start






def rm_dup_pe(dupmark_bam, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(dupmark_bam)))
    # strip extension appended in the previous step
    prefix = strip_ext(prefix,'dupmark')
    nodup_bam = '{}.nodup.bam'.format(prefix)

    cmd1 = 'samtools view -@ {} -F 1804 -f 2 -b {} > {}'
    cmd1 = cmd1.format(
        nth,
        dupmark_bam,
        nodup_bam)
    run_shell_cmd(cmd1)
    return nodup_bam

.
LC_ALL=C sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 \
-o $sorted_file \
$merged_file





#convert to bam
bam_file1="${mapped_file1%.*}.bam"
samtools view \
-S \
-b "$mapped_file1" \
> "$bam_file1"

bam_file2="${mapped_file2%.*}.bam"
samtools view \
-S \
-b "$mapped_file2" \
> "$bam_file2"

# merging split files
merged_bam="${bam_file1%_*}_trimmed.splits.mapped.merged.bam"
samtools merge \
$merged_bam \
$bam_file1 \
$bam_file2

# running in methpipe
merged_bam="${bam_file1%_*}_trimmed.splits.mapped.merged.bam"
sorted_bam="${merged_bam%.*}.sorted.bam"
LC_ALL=C sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 \
-o "$sorted_bam"
"$merged_bam" \

# remove dupilicates
dremove_stats="${sorted_bam}.dremove.stat.txt"
dremove_output="${sorted_bam}.dremove"
/home/segil_lab/tools/methpipe-3.4.3/bin/duplicate-remover \
-S $dremove_stats \
-o $dremove_output \
$sorted_bam



















trim-galore \
--illumina \
ENCFF810GKB.fastq.gz
    # output
    ENCFF810GKB_trimmed.fq.gz


# file: rep1/ENCFF810GKB_trimmed.fq.gz
#split files
file_name="/media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep1/ENCFF810GKB_trimmed.fq.gz"
lines=$(pigz -cd -p 4 "$file_name" | wc -l)
half_lines=$(echo "$lines / 2" | bc)
rm_lines=$(echo "$lines - $half_lines" | bc)

pigz -cd -p 4 "$file_name" | head -"${half_lines}" | \
pigz -p 4 > /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep1/ENCFF810GKB_trimmed_1.fq.gz

pigz -cd -p 4 "$file_name" | tail -"${rm_lines}" | \
pigz -p 4 > /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep1/ENCFF810GKB_trimmed_2.fq.gz

######### running in methpipe2











# map split files
/home/segil_lab/tools/walt/bin/walt \
-i /media/segil_lab/SegilRaid/staging/Duc/mm10_bis/mm10/mm10.fa.dbindex \
-r <(zcat /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep2/ENCFF110EYU_trimmed_1.fq.gz) \
-N 5000000 \
-u -a \
-o /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep2/ENCFF110EYU_trimmed_1_mapping.sam
        # [TOTAL NUMBER OF READS: 317339441]
        # [UNIQUELY MAPPED READS: 253781540 (79.97%)]
        # [AMBIGUOUS MAPPED READS: 32923495 (10.37%)]
        # [UNMAPPED READS: 30634406 (9.65%)]
        #
        #
        # [READS SHORTER THAN 38 ARE IGNORED]
        # [4594370 (1.45%) READS ARE SHORTER THAN 38]

/home/segil_lab/tools/walt/bin/walt \
-i /media/segil_lab/SegilRaid/staging/Duc/mm10_bis/mm10/mm10.fa.dbindex \
-r <(zcat /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep2/ENCFF110EYU_trimmed_2.fq.gz) \
-N 5000000 \
-u -a \
-o /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep2/ENCFF110EYU_trimmed_2_mapping.sam
          # [TOTAL NUMBER OF READS: 317339441]
          # [UNIQUELY MAPPED READS: 253097741 (79.76%)]
          # [AMBIGUOUS MAPPED READS: 33087329 (10.43%)]
          # [UNMAPPED READS: 31154371 (9.82%)]
          #
          #
          # [READS SHORTER THAN 38 ARE IGNORED]
          # [3194320 (1.01%) READS ARE SHORTER THAN 38]

          # .sam lines : 590266037

#convert to bam
samtools view \
-S \
-b /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep2/ENCFF110EYU_trimmed_1_mapping.sam \
> /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep2/ENCFF110EYU_trimmed_1_mapping.bam

samtools view \
-S \
-b /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep2/ENCFF110EYU_trimmed_2_mapping.sam \
> /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep2/ENCFF110EYU_trimmed_2_mapping.bam

samtools merge \
/media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep2/ENCFF110EYU_merged.bam \
/media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep2/ENCFF110EYU_trimmed_1_mapping.bam \
/media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep2/ENCFF110EYU_trimmed_2_mapping.bam




# did not run this yet
samtools sort /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep1/ENCFF627KYH_mapping.bam \
-o /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep1/ENCFF627KYH_mapping.sorted.bam

#running in tmux - methpipe_split
samtools view \
-S \
-b /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep2/ENCFF110EYU_trimmed_1_mapping.sam \
> /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep2/ENCFF110EYU_trimmed_1_mapping.bam

# running in tmux
samtools view \
-S \
-b /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep2/ENCFF110EYU_trimmed_2_mapping.sam \
> /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep2/ENCFF110EYU_trimmed_2_mapping.bam

samtools sort /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep1/ENCFF627KYH_mapping.bam \
-o /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep1/ENCFF627KYH_mapping.sorted.bam




bedtools bamtobed \
-i /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep1/ENCFF627KYH_mapping.sam \
| head -10








tail numberoflines

zcat /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10_5/rep2/ENCFF110EYU_trimmed.fq.gz | \
for i in rep2/ENCFF110EYU_trimmed.fq.gz; do \
    split -a 3 -d -l 12000000 ${i} reads_split/$(basename $i); done





scp /Users/johnnguyen/Downloads/Ecker_Protocol_Summary.png jd@128.125.247.142:/media/segil_lab/SegilRaid/staging/Duc/wgbs/
scp /Users/johnnguyen/Documents/1-Segil_Lab/Data/Sample_84_1/84_1_S8_R1_001.fastq.gz jd@128.125.247.142:/media/segil_lab/SegilRaid/staging/Duc/wgbs/ &&\
scp /Users/johnnguyen/Documents/1-Segil_Lab/Data/Sample_84_1/84_1_S8_R2_001.fastq.gz jd@128.125.247.142:/media/segil_lab/SegilRaid/staging/Duc/wgbs/



scp /dev/disk2s2/Volumes/SEQUENCING /WGBS-Seq/Variants_Remi_122017/Sample_84_1/84_1_S8_R1_001.fastq.gz jd@128.125.247.142:/media/segil_lab/SegilRaid/staging/Duc/wgbs/ &&\
scp /WGBS-Seq/Variants_Remi_122017/Sample_84_1/84_1_S8_R2_001.fastq.gz jd@128.125.247.142:/media/segil_lab/SegilRaid/staging/Duc/wgbs/

/dev/disk2s2/Volumes/SEQUENCING /WGBS-Seq/Variants_Remi_122017/Sample_84_1

/Users/johnnguyen/Downloads
