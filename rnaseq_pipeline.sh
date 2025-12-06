#!/bin/bash

# RNA-Seq 上游分析流程
# 包括：质控 -> 比对 -> 定量

set -e  # 遇到错误即退出

# ==================== 配置参数 ====================
# 工作目录
BASE_DIR="/data/mzwang/sheep"
RAW_DATA_DIR="${BASE_DIR}/PRJNA1293839"
REF_DIR="${BASE_DIR}/reference_genome"

# 输出目录
OUTPUT_DIR="${BASE_DIR}/analysis"
FASTP_DIR="${OUTPUT_DIR}/01_fastp"
QC_DIR="${OUTPUT_DIR}/02_qc_reports"
HISAT2_INDEX_DIR="${OUTPUT_DIR}/03_hisat2_index"
BAM_DIR="${OUTPUT_DIR}/04_bam"
COUNT_DIR="${OUTPUT_DIR}/05_counts"
LOG_DIR="${OUTPUT_DIR}/logs"

# 参考基因组文件
GENOME_FA="${REF_DIR}/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna.gz"
GTF_FILE="${REF_DIR}/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.gtf.gz"

# fastp 参数
FASTP_PARAMS="-f 3 -t 3 -l 36 -r -W 4 -M 15"

# 线程数
THREADS=8

# STAR 索引目录
STAR_INDEX_DIR="${OUTPUT_DIR}/03_star_index"

# ==================== 创建输出目录 ====================
echo "Creating output directories..."
mkdir -p ${FASTP_DIR} ${QC_DIR} ${STAR_INDEX_DIR} ${BAM_DIR} ${COUNT_DIR} ${LOG_DIR}

# ==================== 获取样本列表 ====================
cd ${RAW_DATA_DIR}
SAMPLES=($(ls *_1.fastq.gz | sed 's/_1.fastq.gz//'))
echo "Found ${#SAMPLES[@]} samples: ${SAMPLES[@]}"

# ==================== 1. fastp 质控 ====================
echo "================================"
echo "Step 1: Quality control with fastp"
echo "================================"

for sample in "${SAMPLES[@]}"; do
    echo "Processing sample: ${sample}"
    
    fastp \
        -i ${RAW_DATA_DIR}/${sample}_1.fastq.gz \
        -I ${RAW_DATA_DIR}/${sample}_2.fastq.gz \
        -o ${FASTP_DIR}/${sample}_clean_1.fastq.gz \
        -O ${FASTP_DIR}/${sample}_clean_2.fastq.gz \
        ${FASTP_PARAMS} \
        -h ${QC_DIR}/${sample}_fastp.html \
        -j ${QC_DIR}/${sample}_fastp.json \
        --thread ${THREADS} \
        2>&1 | tee ${LOG_DIR}/${sample}_fastp.log
    
    echo "Finished: ${sample}"
done

echo "fastp completed for all samples!"

# ==================== 2. 构建 STAR 索引 ====================
echo "================================"
echo "Step 2: Building STAR index"
echo "================================"

# 解压基因组和GTF文件（如果需要）
if [ ! -f "${STAR_INDEX_DIR}/genome.fa" ]; then
    echo "Decompressing genome file..."
    gunzip -c ${GENOME_FA} > ${STAR_INDEX_DIR}/genome.fa
fi

if [ ! -f "${STAR_INDEX_DIR}/annotation.gtf" ]; then
    echo "Decompressing GTF file..."
    gunzip -c ${GTF_FILE} > ${STAR_INDEX_DIR}/annotation.gtf
fi

# 构建STAR索引
if [ ! -f "${STAR_INDEX_DIR}/SA" ]; then
    echo "Building STAR index..."
    STAR \
        --runMode genomeGenerate \
        --runThreadN ${THREADS} \
        --genomeDir ${STAR_INDEX_DIR} \
        --genomeFastaFiles ${STAR_INDEX_DIR}/genome.fa \
        --sjdbGTFfile ${STAR_INDEX_DIR}/annotation.gtf \
        --sjdbOverhang 149 \
        2>&1 | tee ${LOG_DIR}/star_build.log
else
    echo "STAR index already exists, skipping..."
fi

# ==================== 3. STAR 双端比对（Two-pass mode）====================
echo "================================"
echo "Step 3: Alignment with STAR (Two-pass mode)"
echo "================================"

for sample in "${SAMPLES[@]}"; do
    echo "Aligning sample: ${sample}"
    
    # 创建样本专用输出目录
    SAMPLE_OUT_DIR="${BAM_DIR}/${sample}"
    mkdir -p ${SAMPLE_OUT_DIR}
    
    STAR \
        --genomeDir ${STAR_INDEX_DIR} \
        --sjdbGTFfile ${STAR_INDEX_DIR}/annotation.gtf \
        --twopassMode Basic \
        --readFilesIn ${FASTP_DIR}/${sample}_clean_1.fastq.gz ${FASTP_DIR}/${sample}_clean_2.fastq.gz \
        --outFileNamePrefix ${SAMPLE_OUT_DIR}/${sample}_ \
        --runThreadN ${THREADS} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterMismatchNmax 3 \
        --outSAMunmapped Within \
        --chimSegmentMin 10 \
        --chimOutType Junctions \
        --chimOutJunctionFormat 1 \
        --outFilterType BySJout \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outSAMattributes NH HI AS nM NM \
        --outSAMattrRGline ID:${sample} SM:${sample} \
        2>&1 | tee ${LOG_DIR}/${sample}_star.log
    
    # 重命名输出文件
    mv ${SAMPLE_OUT_DIR}/${sample}_Aligned.sortedByCoord.out.bam ${BAM_DIR}/${sample}.sorted.bam
    
    # 建立BAM索引
    samtools index ${BAM_DIR}/${sample}.sorted.bam
    
    # 移动其他输出文件到样本目录
    if [ -f "${SAMPLE_OUT_DIR}/${sample}_ReadsPerGene.out.tab" ]; then
        mv ${SAMPLE_OUT_DIR}/${sample}_ReadsPerGene.out.tab ${SAMPLE_OUT_DIR}/${sample}_star_counts.tab
    fi
    
    echo "Finished: ${sample}"
done

echo "STAR alignment completed for all samples!"

# ==================== 4. featureCounts 定量 ====================
echo "================================"
echo "Step 4: Gene-level quantification with featureCounts"
echo "================================"

# 准备所有 BAM 文件列表
BAM_FILES=$(ls ${BAM_DIR}/*.sorted.bam | tr '\n' ' ')

echo "Running featureCounts for gene-level quantification..."
featureCounts \
    -T ${THREADS} \
    -p \
    -B \
    -C \
    -t exon \
    -g gene_id \
    -a ${STAR_INDEX_DIR}/annotation.gtf \
    -o ${COUNT_DIR}/gene_counts.txt \
    ${BAM_FILES} \
    2>&1 | tee ${LOG_DIR}/featureCounts.log

echo "featureCounts completed!"

# 生成简化版计数矩阵（只保留基因ID和计数）
echo "Generating simplified count matrix..."
cut -f1,7- ${COUNT_DIR}/gene_counts.txt | \
    sed '1d' > ${COUNT_DIR}/gene_counts_matrix.txt

echo "Count matrix saved to: ${COUNT_DIR}/gene_counts_matrix.txt"

# ==================== 5. 生成比对统计 ====================
echo "================================"
echo "Step 5: Generate alignment statistics"
echo "================================"

for sample in "${SAMPLES[@]}"; do
    samtools flagstat ${BAM_DIR}/${sample}.sorted.bam > ${LOG_DIR}/${sample}_flagstat.txt
done

echo "================================"
echo "Pipeline completed successfully!"
echo "================================"
echo "Results location:"
echo "  - Clean reads: ${FASTP_DIR}"
echo "  - QC reports: ${QC_DIR}"
echo "  - STAR index: ${STAR_INDEX_DIR}"
echo "  - BAM files: ${BAM_DIR}"
echo "  - Gene counts: ${COUNT_DIR}/gene_counts.txt"
echo "  - Count matrix: ${COUNT_DIR}/gene_counts_matrix.txt"
echo "  - Logs: ${LOG_DIR}"
echo ""
echo "Next steps:"
echo "  - Check alignment statistics in ${LOG_DIR}/*_star.log"
echo "  - Review featureCounts summary in ${COUNT_DIR}/gene_counts.txt.summary"
echo "  - Use gene_counts_matrix.txt for downstream differential expression analysis"
