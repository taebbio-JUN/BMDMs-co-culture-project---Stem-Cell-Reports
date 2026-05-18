PATH=$PATH:$HOME/bin:/home/user/RNA-seq_data/samtools/program/bin
export BT2_HOME="/home/user/RNA-seq_data/bowtie2/bowtie2-2.5.4-linux-x86_64"
dir="/home/user/RNA-seq_data"
trimmomatic=${dir}/Trimmomatic-0.40
bowtie2=${dir}/bowtie2/bowtie2-2.5.4-linux-x86_64
samtools=${dir}/samtools/program/bin
stringtie=${dir}/stringtie_v3.0.3_git_202602062242/stringtie
gffcompare=${dir}/gffcompare-0.12.10.Linux_x86_64
source_dir=${dir}/20260428_Physiology
trimmed=${source_dir}/trimmed
trimmed_unpair=${source_dir}/trimmed/unpaired
trimmed_pair=${source_dir}/trimmed/paired

################For NCBI gene ID (Refseq) GTF################
#gtf_dir=${dir}/GRCm39/gtf/NCBI_RefSeq_annotation_GTF/GRCm39_NCBI_RefSeq_Jun_2025_GCF_00001635_27/ncbi_dataset/data/GCF_000001635.27
#annotation_gtf="GCF_000001635.27_GRCm39_genomic.gtf"
#bt2_name="GRCm39"
#gene_coord_sorted=${trimmed_pair}/${bt2_name}/gene_coord_sorted
#assembled=${gene_coord_sorted}/assembled_ncbi

#################For Gencode GTF ################
####- To use Bowtie 2, you must build the index files (.bt2) yourself.
gtf_dir=${dir}/GRCm39/gtf/Gencode_GTF
bed_dir=${dir}/GRCm39/bed
original_bt2_dir=${dir}/GRCm39/GRCm39_m38_Gencode_for_bowtie2
annotation_gtf="gencode.vM38.chr_patch_hapl_scaff.annotation.gtf"
#annotation_gtf="gencode.vM38.annotation.gtf"
#annotation_gtf="gencode.vM38.basic.annotation.gtf"
GRCMm39_bed="GRCm39_knowngene_from_genome_ucsc_edu.bed"
fasta_file_location=${dir}/GRCm39/gtf/Gencode_GTF
fasta_file="GRCm39.genome.fa"
bt2_name="GRCm39_gencode_m38_Sep_2025"
gene_coord_sorted=${trimmed_pair}/${bt2_name}/gene_coord_sorted
assembled=${gene_coord_sorted}/assembled_gencode
################################################################

#################For mm10 UCSC GTF ################
####- To use Bowtie 2, the index files must be downloaded
#gtf_dir=${dir}/mm10/gtf
#bed_dir=${dir}/mm10/bed
#annotation_gtf="mm10.ensGene.gtf"
#GRCMm39_bed="GRCm39_knowngene_from_genome_ucsc_edu.bed"
#fasta_file_location=${dir}/mm10/gtf
#fasta_file="GRCm39.genome.fa"
#bt2_name="mm10"
#gene_coord_sorted=${trimmed_pair}/${bt2_name}/gene_coord_sorted
#assembled=${gene_coord_sorted}/assembled_gencode
################################################################
mkdir ${trimmed}
mkdir ${trimmed_pair}
mkdir ${trimmed_unpair}
mkdir ${trimmed_pair}/${bt2_name}
cp ${original_bt2_dir}/*.bt2 ${trimmed_pair}/${bt2_name}

mkdir ${gene_coord_sorted}
mkdir ${assembled}

## Checksum 
echo "Checksum Check : md5sum_double_check.txt"
md5sum ${source_dir}/*.gz > md5sum_double_check.txt


## Trimmomatic trimming
#java -jar trimmomatic-0.40.jar PE hetero-1intact_1.fastq.gz hetero-1intact_2.fastq.gz #hetero-1intact_1_paired.fastq.gz hetero-1intact_1_unpaired.fastq.gz hetero-1intact_2_paired.fastq.gz #hetero-1intact_2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

#echo "==========Trimmomatic start=========="
for sample in `ls ${source_dir}/*_1.fastq.gz`
do
base=$(basename $sample "_1.fastq.gz")
echo "----------${base}_1.fastq.gz_START--------"
java -jar ${trimmomatic}/trimmomatic-0.40.jar PE ${source_dir}/${base}_1.fastq.gz ${source_dir}/${base}_2.fastq.gz ${source_dir}/trimmed/paired/${base}_1_paired.fastq.gz ${source_dir}/trimmed/unpaired/${base}_1_unpaired.fastq.gz ${source_dir}/trimmed/paired/${base}_2_paired.fastq.gz ${source_dir}/trimmed/unpaired/${base}_2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
echo "----------${base}_1.fastq.gz_END----------"
done
#echo "==========Trimmomatic end=========="
#./notify.sh


## Bowtie2 index building from Gencode FASTA file
#echo "==========Bowtie2-build bt2 start=========="
#echo "----------${fasta_file}_START--------"
#echo "----------$(date +"%Y-%m-%d_%H-%M-%S")--------"
#mkdir ${fasta_file_location}/bt2
#(cd ${fasta_file_location}; ${bowtie2}/bowtie2-build ${fasta_file_location}/${fasta_file} ${bt2_gencode_name})
#echo "----------${fasta_file}_END----------"
#echo "==========Bowtie2-build bt2 end=========="


## Bowtie2 mapping
echo "==========Bowtie2 mapping start=========="
mkdir ${trimmed_pair}/${bt2_name}
for sample in `ls ${trimmed_pair}/*_1_paired.fastq.gz`
do
base=$(basename $sample "_1_paired.fastq.gz")
echo "----------${base}_1_paired.fastq.gz_START--------"
echo "----------$(date +"%Y-%m-%d_%H-%M-%S")--------"
(cd ${trimmed_pair}; ${bowtie2}/bowtie2 -p 3 -x ${bt2_name} -1 ${base}_1_paired.fastq.gz -2 ${base}_2_paired.fastq.gz | ${samtools}/samtools view -@ 3 -bS -o ${trimmed_pair}/${bt2_name}/${base}.bam)
echo "----------${base}_1_paired.fastq.gz_END----------"
done
echo "==========Bowtie2 mapping end=========="
./notify.sh


## Samtools Sorting
#(cd ${trimmed_pair};  ${samtools}/samtools sort ${trimmed_pair}/GRCm39/sni_ngs_ryukyus_bam.bam -o ${trimmed_pair}/GRCm39/sni_ngs_ryukyus_bam.sorted.bam )
echo "==========Samtools sorting start=========="
for sample in `ls ${trimmed_pair}/${bt2_name}/*.bam`
do
base=$(basename $sample ".bam")
echo "----------${base}_sorted.bam_START--------"
echo "----------$(date +"%Y-%m-%d_%H-%M-%S")--------"
(cd ${trimmed_pair}; ${samtools}/samtools sort -@ 4 ${trimmed_pair}/${bt2_name}/${base}.bam -o ${trimmed_pair}/${bt2_name}/${base}_sorted.bam)
echo "----------${base}_sorted.bam_END----------"
done
echo "==========Samtools sorting end=========="
./notify.sh

(cd ${gene_coord_sorted}; mv ../*_sorted.bam .) 

## Samtools Indexing
#(cd ${gene_coord_sorted}; ${samtools}/samtools index ${gene_coord_sorted}/xxx_sorted.bam)
echo "==========Samtools indexing start=========="
for sample in `ls ${gene_coord_sorted}/*_sorted.bam`
do
base=$(basename $sample ".bam")
echo "----------${base}_sorted.bam_START--------"
echo "----------$(date +"%Y-%m-%d_%H-%M-%S")--------"
(cd ${gene_coord_sorted}; ${samtools}/samtools index ${gene_coord_sorted}/${base}.bam)
echo "----------${base}_sorted.bam_END----------"
done
echo "==========Samtools indexing end=========="
./notify.sh


## Stringtie assembly
cp ${gtf_dir}/${annotation_gtf} ${assembled}

echo "==========Stringtie assembly start=========="
echo "----------File and label list--------"
#declare -A file_label=(
["1K1R_sorted"]="1K1R"
["1K2R_sorted"]="1K2R"
["1K3R_sorted"]="1K3R"
["5K1R_sorted"]="5K1R"
["5K2R_sorted"]="5K2R"
["5K3R_sorted"]="5K3R"
["BMDM1_sorted"]="BMDM1"
["BMDM2_sorted"]="BMDM2"
["BMDM3_sorted"]="BMDM3"
["Con1R_sorted"]="Con1R"
["Con2R_sorted"]="Con2R"
["Con3R_sorted"]="Con3R"
["TNFa1K1_sorted"]="TNFa1K1"
["TNFa1K2_sorted"]="TNFa1K2"
["TNFa1K3_sorted"]="TNFa1K3"
["TNFa2_sorted"]="TNFa2"
["TNFa3_sorted"]="TNFa3"
["TNFa4_sorted"]="TNFa4"
["TNFa5K1_sorted"]="TNFa5K1"
["TNFa5K2_sorted"]="TNFa5K2"
["TNFa5K3_sorted"]="TNFa5K3"
["TNFaBMDM1_sorted"]="TNFaBMDM1"
["TNFaBMDM2_sorted"]="TNFaBMDM2"
["TNFaBMDM3_sorted"]="TNFaBMDM3"
)

for i in "${!file_label[@]}"
do
  echo "file  : $i.bam"
  echo "label: ${file_label[$i]}"
done
echo "-------------------------------------"

# create assembly per sample using 4 threads
(cd ${gene_coord_sorted};
for i in "${!file_label[@]}"
do
  echo "-----------${i}.bam - ${file_label[$i]}_START--------"
  echo "----------$(date +"%Y-%m-%d_%H-%M-%S")--------"
    ${stringtie}/stringtie -p 4 -N ${gene_coord_sorted}/${i}.bam -G ${assembled}/${annotation_gtf} -o ${assembled}/${i}_rm_nascent.gtf 
  echo "----------${i}.bam - ${file_label[$i]}_END----------"
done)
echo "==========Stringtie assembly end=========="

## Stringtie merge

### Create mergelist from gtf files
for i in "${!file_label[@]}"
do
  echo "${i}_rm_nascent.gtf" >> ${assembled}/mergelist.txt
done


echo "==========Stringtie merge start=========="
echo "-----------all_merged.gtf_START--------"
echo "----------$(date +"%Y-%m-%d_%H-%M-%S")--------"
(cd ${assembled}; ${stringtie}/stringtie --merge --rf -o ${assembled}/all_merged.gtf -G ${assembled}/${annotation_gtf} ${assembled}/mergelist.txt)
echo "-----------all_merged.gtf_END--------"
echo "==========Stringtie merge end=========="

#check out the transcripts
cat ${assembled}/all_merged.gtf | head

# how many transcripts?
cat ${assembled}/all_merged.gtf | grep -v "^#" | awk '$3=="transcript" {print}' | wc -l


## Gffcompare
echo "==========Gffcompare start=========="
echo "-----------all_merged.gtf_START--------"
echo "----------$(date +"%Y-%m-%d_%H-%M-%S")--------"
(cd ${assembled}; ${gffcompare}/gffcompare -r ${assembled}/${annotation_gtf} -o merged ${assembled}/all_merged.gtf)

(cd ${assembled}; cat merged.stats)
echo "-----------all_merged.gtf_END--------"
echo "==========Gffcompare end=========="

mkdir ${assembled}/ballgown4
echo "==========Stringtie output ballgown2 start=========="
for i in "${!file_label[@]}"
do
  echo "-----------${file_label[$i]}.gtf_START--------"
  echo "----------$(date +"%Y-%m-%d_%H-%M-%S")--------"
  ${stringtie}/stringtie -e -B -G ${assembled}/all_merged.gtf -o ${assembled}/ballgown4/${file_label[$i]}/${file_label[$i]}.gtf ${gene_coord_sorted}/${i}.bam
  echo "----------${file_label[$i]}.gtf_END----------"
done
echo "==========Stringtie output ballgown2 end=========="

## create CSV for DESeq2 or edgeR
mv ${assembled}/ballgown4/gtftocsv_list.txt ${assembled}/ballgown3/gtftocsv_list.old.txt
for i in "${!file_label[@]}"
do
 echo "${file_label[$i]} ./${file_label[$i]}/${file_label[$i]}.gtf" >> ${assembled}/ballgown4/gtftocsv_list.txt
done
echo "==========Stringtie prepDE.py start=========="
echo "----------$(date +"%Y-%m-%d_%H-%M-%S")--------"
# -l 202 : paired , readLength 101 bp
(cd ${assembled}/ballgown4; python3 ${stringtie}/prepDE.py3 -l 202 -i ${assembled}/ballgown4/gtftocsv_list.txt -v)
echo "==========Stringtie prepDE.py end=========="
