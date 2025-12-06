#!/bin/bash
# RNA-seq Analysis Pipeline for Pig Data
# 流程: Trimmomatic -> STAR -> featureCounts
# 使用方法: 将此脚本放在项目文件夹下（如PRJNA1139427），直接运行即可

set -e  # 遇到错误立即退出

# ===================== 配置部分 =====================

# 自动检测当前项目目录
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_NAME="$(basename ${PROJECT_DIR})"

# 原始数据目录（当前项目文件夹）
RAW_DATA_DIR="${PROJECT_DIR}"

# 参考基因组目录（上层目录的reference_genome）
REF_DIR="${PROJECT_DIR}/../reference_genome"

# 软件路径
TRIMMOMATIC_JAR="/data/mzwang/.conda/envs/pig/share/trimmomatic-0.39-2/trimmomatic.jar"
TRIMMOMATIC_ADAPTERS="/data/mzwang/.conda/envs/pig/share/trimmomatic-0.39-2/adapters"
STAR_BIN="STAR"
FEATURECOUNTS_BIN="featureCounts"
STRINGTIE_BIN="stringtie"

# 参考基因组文件
REF_GENOME_FA="${REF_DIR}/GCF_000003025.6_Sscrofa11.1_genomic.renamed.fna"
GTF_FILE="${REF_DIR}/Sus_scrofa.Sscrofa11.1.100.gtf"
STAR_INDEX_DIR="${REF_DIR}/star_index"

# 输出目录（统一放在results文件夹下）
RESULTS_DIR="${PROJECT_DIR}/results"
TRIMMED_DIR="${RESULTS_DIR}/trimmed_data"
MAPPED_DIR="${RESULTS_DIR}/mapped_data"
COUNTS_DIR="${RESULTS_DIR}/gene_counts"
STRINGTIE_DIR="${RESULTS_DIR}/stringtie_abundance"

# 线程数设置
THREADS_TRIMMOMATIC=8
THREADS_STAR=16
THREADS_FEATURECOUNTS=8

# 接头文件
ADAPTER_FILE="${TRIMMOMATIC_ADAPTERS}/TruSeq3-PE-2.fa"

# ===================== 验证文件存在性 =====================

echo "======================================"
echo "RNA-seq分析流程启动"
echo "项目名称: ${PROJECT_NAME}"
echo "项目目录: ${PROJECT_DIR}"
echo "开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "======================================"

# 检查参考基因组文件
if [ ! -f "${REF_GENOME_FA}" ]; then
    echo "错误: 参考基因组文件不存在: ${REF_GENOME_FA}"
    exit 1
fi

if [ ! -f "${GTF_FILE}" ]; then
    echo "错误: GTF注释文件不存在: ${GTF_FILE}"
    exit 1
fi

echo "参考基因组检查通过"

# ===================== 创建输出目录 =====================

mkdir -p ${TRIMMED_DIR}/logs
mkdir -p ${MAPPED_DIR}/logs
mkdir -p ${COUNTS_DIR}/logs
mkdir -p ${STRINGTIE_DIR}/logs

echo "结果目录创建完成: ${RESULTS_DIR}"

# ===================== 步骤1: 构建STAR索引 =====================

if [ ! -d "${STAR_INDEX_DIR}" ] || [ -z "$(ls -A ${STAR_INDEX_DIR})" ]; then
    echo ""
    echo "======================================"
    echo "步骤1: 构建STAR索引（预计30-60分钟）"
    echo "======================================"
    mkdir -p ${STAR_INDEX_DIR}
    
    ${STAR_BIN} --runThreadN 20 \
        --runMode genomeGenerate \
        --genomeDir ${STAR_INDEX_DIR} \
        --genomeFastaFiles ${REF_GENOME_FA} \
        --sjdbGTFfile ${GTF_FILE} \
        --sjdbOverhang 100 \
        2>&1 | tee ${REF_DIR}/star_index.log
    
    echo "STAR索引构建完成"
else
    echo ""
    echo "步骤1: STAR索引已存在，跳过构建"
fi

# 验证STAR索引完整性
if [ ! -f "${STAR_INDEX_DIR}/genomeParameters.txt" ]; then
    echo "错误: STAR索引不完整，缺少genomeParameters.txt文件"
    echo "请检查索引构建日志: ${REF_DIR}/star_index.log"
    exit 1
fi

# ===================== 步骤2: 获取样本列表 =====================

echo ""
echo "======================================"
echo "步骤2: 识别样本"
echo "======================================"

cd ${RAW_DATA_DIR}
SAMPLES=($(ls *_1.fastq.gz 2>/dev/null | sed 's/_1.fastq.gz//'))

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "错误: 在 ${RAW_DATA_DIR} 中未找到任何 *_1.fastq.gz 文件"
    exit 1
fi

echo "找到 ${#SAMPLES[@]} 个样本"
echo "样本列表（前10个）:"
printf '%s\n' "${SAMPLES[@]:0:10}"

# ===================== 步骤3-5: 处理每个样本 =====================

SAMPLE_COUNT=0
for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    echo ""
    echo "======================================"
    echo "处理样本 ${SAMPLE_COUNT}/${#SAMPLES[@]}: ${SAMPLE}"
    echo "时间: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "======================================"
    
    # 文件路径定义
    R1_RAW="${RAW_DATA_DIR}/${SAMPLE}_1.fastq.gz"
    R2_RAW="${RAW_DATA_DIR}/${SAMPLE}_2.fastq.gz"
    R1_CLEAN="${TRIMMED_DIR}/${SAMPLE}_1.clean.fq.gz"
    R2_CLEAN="${TRIMMED_DIR}/${SAMPLE}_2.clean.fq.gz"
    R1_UNPAIRED="${TRIMMED_DIR}/${SAMPLE}_1.unpaired.fq.gz"
    R2_UNPAIRED="${TRIMMED_DIR}/${SAMPLE}_2.unpaired.fq.gz"
    SAMPLE_DIR="${MAPPED_DIR}/${SAMPLE}"
    BAM_FILE="${SAMPLE_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    COUNT_FILE="${COUNTS_DIR}/${SAMPLE}.featureCounts.txt"
    STRINGTIE_GTF="${STRINGTIE_DIR}/${SAMPLE}.stringtie.gtf"
    STRINGTIE_ABUNDANCE="${STRINGTIE_DIR}/${SAMPLE}.gene_abundance.txt"
    
    # ---------- 步骤3: Trimmomatic质控 ----------
    if [ ! -f "${R1_CLEAN}" ] || [ ! -f "${R2_CLEAN}" ]; then
        echo "[${SAMPLE}] 步骤3: Trimmomatic质控..."
        
        java -jar ${TRIMMOMATIC_JAR} PE -phred33 \
            -threads ${THREADS_TRIMMOMATIC} \
            ${R1_RAW} ${R2_RAW} \
            ${R1_CLEAN} ${R1_UNPAIRED} \
            ${R2_CLEAN} ${R2_UNPAIRED} \
            ILLUMINACLIP:${ADAPTER_FILE}:2:30:10:2:TRUE \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36 \
            2>&1 | tee ${TRIMMED_DIR}/logs/${SAMPLE}.trimmomatic.log
        
        echo "[${SAMPLE}] 质控完成"
    else
        echo "[${SAMPLE}] 质控文件已存在，跳过"
    fi
    
    # ---------- 步骤4: STAR比对 ----------
    if [ ! -f "${BAM_FILE}" ]; then
        echo "[${SAMPLE}] 步骤4: STAR比对..."
        
        mkdir -p ${SAMPLE_DIR}
        
        ${STAR_BIN} --runThreadN ${THREADS_STAR} \
            --genomeDir ${STAR_INDEX_DIR} \
            --sjdbGTFfile ${GTF_FILE} \
            --readFilesIn ${R1_CLEAN} ${R2_CLEAN} \
            --readFilesCommand zcat \
            --outFileNamePrefix ${SAMPLE_DIR}/${SAMPLE}_ \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --outFilterMismatchNmax 3 \
            --quantMode TranscriptomeSAM \
            --outSAMattributes NH HI AS NM MD \
            2>&1 | tee ${MAPPED_DIR}/logs/${SAMPLE}.star.log
        
        # 验证BAM文件生成
        if [ ! -f "${BAM_FILE}" ]; then
            echo "错误: BAM文件未生成 - ${BAM_FILE}"
            echo "请检查STAR日志: ${MAPPED_DIR}/logs/${SAMPLE}.star.log"
            exit 1
        fi
        
        echo "[${SAMPLE}] 比对完成"
    else
        echo "[${SAMPLE}] BAM文件已存在，跳过"
    fi
    
    # ---------- 步骤5: featureCounts定量 ----------
    if [ ! -f "${COUNT_FILE}" ]; then
        echo "[${SAMPLE}] 步骤5: featureCounts定量..."
        
        ${FEATURECOUNTS_BIN} \
            -T ${THREADS_FEATURECOUNTS} \
            -p \
            -t exon \
            -g gene_id \
            -a ${GTF_FILE} \
            -o ${COUNT_FILE} \
            ${BAM_FILE} \
            2>&1 | tee ${COUNTS_DIR}/logs/${SAMPLE}.featureCounts.log
        
        # 验证计数文件生成
        if [ ! -f "${COUNT_FILE}" ]; then
            echo "错误: 计数文件未生成 - ${COUNT_FILE}"
            exit 1
        fi
        
        echo "[${SAMPLE}] 定量完成"
    else
        echo "[${SAMPLE}] 计数文件已存在，跳过"
    fi
    
    # ---------- 步骤6: StringTie TPM定量 ----------
    if [ ! -f "${STRINGTIE_ABUNDANCE}" ]; then
        echo "[${SAMPLE}] 步骤6: StringTie TPM定量..."
        
        ${STRINGTIE_BIN} \
            -e \
            -B \
            -p ${THREADS_FEATURECOUNTS} \
            -G ${GTF_FILE} \
            -o ${STRINGTIE_GTF} \
            -A ${STRINGTIE_ABUNDANCE} \
            ${BAM_FILE} \
            2>&1 | tee ${STRINGTIE_DIR}/logs/${SAMPLE}.stringtie.log
        
        # 验证输出文件生成
        if [ ! -f "${STRINGTIE_ABUNDANCE}" ]; then
            echo "错误: StringTie丰度文件未生成 - ${STRINGTIE_ABUNDANCE}"
            exit 1
        fi
        
        echo "[${SAMPLE}] StringTie定量完成"
    else
        echo "[${SAMPLE}] StringTie丰度文件已存在，跳过"
    fi
    
    echo "[${SAMPLE}] 处理完成！"
done

# ===================== 步骤6: 合并表达矩阵 =====================

echo ""
echo "======================================"
echo "步骤7: 合并所有样本的表达矩阵"
echo "======================================"

MERGED_FILE="${COUNTS_DIR}/all_samples_counts_matrix.txt"

# 提取第一个样本的基因信息列
FIRST_SAMPLE="${SAMPLES[0]}"
FIRST_COUNT_FILE="${COUNTS_DIR}/${FIRST_SAMPLE}.featureCounts.txt"

echo "提取基因信息..."
grep -v "^#" ${FIRST_COUNT_FILE} | cut -f1-6 > ${MERGED_FILE}.tmp

# 提取每个样本的计数列
echo "提取各样本计数数据..."
for SAMPLE in "${SAMPLES[@]}"; do
    COUNT_FILE="${COUNTS_DIR}/${SAMPLE}.featureCounts.txt"
    echo "  处理: ${SAMPLE}"
    grep -v "^#" ${COUNT_FILE} | tail -n +2 | awk '{print $NF}' > ${COUNTS_DIR}/${SAMPLE}.counts.tmp
done

# 合并所有列
echo "合并所有数据..."
paste ${MERGED_FILE}.tmp ${COUNTS_DIR}/*.counts.tmp > ${MERGED_FILE}

# 清理临时文件
rm -f ${MERGED_FILE}.tmp ${COUNTS_DIR}/*.counts.tmp

echo "表达矩阵生成完成: ${MERGED_FILE}"

# ===================== 步骤8: 合并StringTie TPM矩阵 =====================

echo ""
echo "======================================"
echo "步骤8: 合并StringTie TPM表达矩阵"
echo "======================================"

MERGED_TPM_FILE="${STRINGTIE_DIR}/all_samples_TPM_matrix.txt"
MERGED_FPKM_FILE="${STRINGTIE_DIR}/all_samples_FPKM_matrix.txt"

# 创建临时文件存储所有样本数据
TEMP_DIR="${STRINGTIE_DIR}/temp_merge"
mkdir -p ${TEMP_DIR}

echo "提取各样本TPM和FPKM数据..."

# 提取第一个样本的基因ID和名称作为基础
FIRST_SAMPLE="${SAMPLES[0]}"
FIRST_ABUNDANCE="${STRINGTIE_DIR}/${FIRST_SAMPLE}.gene_abundance.txt"

# 提取基因信息（Gene ID和Gene Name）
tail -n +2 ${FIRST_ABUNDANCE} | awk -F'\t' '{print $1"\t"$2}' > ${TEMP_DIR}/gene_info.txt

# 提取每个样本的TPM和FPKM值
for SAMPLE in "${SAMPLES[@]}"; do
    ABUNDANCE_FILE="${STRINGTIE_DIR}/${SAMPLE}.gene_abundance.txt"
    echo "  处理: ${SAMPLE}"
    
    # 提取TPM值（第9列）
    tail -n +2 ${ABUNDANCE_FILE} | awk -F'\t' '{print $9}' > ${TEMP_DIR}/${SAMPLE}.tpm.tmp
    
    # 提取FPKM值（第8列）
    tail -n +2 ${ABUNDANCE_FILE} | awk -F'\t' '{print $8}' > ${TEMP_DIR}/${SAMPLE}.fpkm.tmp
done

# 创建TPM矩阵表头
echo -ne "Gene_ID\tGene_Name" > ${MERGED_TPM_FILE}
for SAMPLE in "${SAMPLES[@]}"; do
    echo -ne "\t${SAMPLE}" >> ${MERGED_TPM_FILE}
done
echo "" >> ${MERGED_TPM_FILE}

# 创建FPKM矩阵表头
echo -ne "Gene_ID\tGene_Name" > ${MERGED_FPKM_FILE}
for SAMPLE in "${SAMPLES[@]}"; do
    echo -ne "\t${SAMPLE}" >> ${MERGED_FPKM_FILE}
done
echo "" >> ${MERGED_FPKM_FILE}

# 合并TPM数据
echo "合并TPM数据..."
paste ${TEMP_DIR}/gene_info.txt ${TEMP_DIR}/*.tpm.tmp >> ${MERGED_TPM_FILE}

# 合并FPKM数据
echo "合并FPKM数据..."
paste ${TEMP_DIR}/gene_info.txt ${TEMP_DIR}/*.fpkm.tmp >> ${MERGED_FPKM_FILE}

# 清理临时文件
rm -rf ${TEMP_DIR}

echo "TPM表达矩阵生成完成: ${MERGED_TPM_FILE}"
echo "FPKM表达矩阵生成完成: ${MERGED_FPKM_FILE}"

# ===================== 生成统计报告 =====================

echo ""
echo "======================================"
echo "分析流程完成！"
echo "完成时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "======================================"
echo ""
echo "输出文件位置："
echo "  质控数据: ${TRIMMED_DIR}"
echo "  比对数据: ${MAPPED_DIR}"
echo "  计数数据: ${COUNTS_DIR}"
echo "  StringTie定量: ${STRINGTIE_DIR}"
echo "  Read counts矩阵: ${MERGED_FILE}"
echo "  TPM表达矩阵: ${MERGED_TPM_FILE}"
echo "  FPKM表达矩阵: ${MERGED_FPKM_FILE}"
echo ""
echo "样本统计："
echo "  总样本数: ${#SAMPLES[@]}"
echo "  成功处理: $(ls ${COUNTS_DIR}/*.featureCounts.txt 2>/dev/null | wc -l)"
echo "  StringTie完成: $(ls ${STRINGTIE_DIR}/*.gene_abundance.txt 2>/dev/null | wc -l)"
echo ""

# 生成比对率统计
echo "比对率统计（前20个样本）："
echo "----------------------------------------"
printf "%-20s %s\n" "样本名" "唯一比对率"
echo "----------------------------------------"
for SAMPLE in "${SAMPLES[@]:0:20}"; do
    LOG_FILE="${MAPPED_DIR}/${SAMPLE}/${SAMPLE}_Log.final.out"
    if [ -f "${LOG_FILE}" ]; then
        UNIQUE_MAPPED=$(grep "Uniquely mapped reads %" ${LOG_FILE} | awk '{print $NF}')
        printf "%-20s %s\n" "${SAMPLE}" "${UNIQUE_MAPPED}"
    fi
done
echo "----------------------------------------"

echo ""
echo "查看完整结果："
echo "  Read counts矩阵: head ${MERGED_FILE}"
echo "  TPM表达矩阵: head ${MERGED_TPM_FILE}"
echo "  FPKM表达矩阵: head ${MERGED_FPKM_FILE}"
echo "  比对日志: cat ${MAPPED_DIR}/logs/<sample>.star.log"
echo "  定量日志: cat ${COUNTS_DIR}/logs/<sample>.featureCounts.log"
echo "  StringTie日志: cat ${STRINGTIE_DIR}/logs/<sample>.stringtie.log"
echo ""
echo "======================================"
echo "所有任务完成！"
echo "======================================"
