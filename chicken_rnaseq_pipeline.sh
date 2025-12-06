#!/bin/bash
# RNA-seq Analysis Pipeline for Chicken Data
# 流程: fastp -> HISAT2 -> StringTie -> featureCounts
# 使用方法: 将此脚本放在项目文件夹下，直接运行即可

set -e  # 遇到错误立即退出

# ===================== 配置部分 =====================

# 自动检测当前项目目录
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_NAME="$(basename ${PROJECT_DIR})"

# 原始数据目录（当前项目文件夹）
RAW_DATA_DIR="${PROJECT_DIR}"

# 参考基因组目录（上层目录的genome_reference）
REF_DIR="${PROJECT_DIR}/../genome_reference"

# 软件路径（假设都在PATH中可以直接调用）
FASTP_BIN="fastp"
HISAT2_BIN="hisat2"
SAMTOOLS_BIN="samtools"
STRINGTIE_BIN="stringtie"
FEATURECOUNTS_BIN="featureCounts"

# 参考基因组文件
GENOME_INDEX="${REF_DIR}/huxu_hisat2_index"
GTF_FILE="${REF_DIR}/20250528_chicken_final_noUTR.gtf"

# 输出目录（统一放在results文件夹下）
RESULTS_DIR="${PROJECT_DIR}/results"

# 线程数设置
THREADS=16

# ===================== 验证文件存在性 =====================

echo "======================================"
echo "RNA-seq分析流程启动"
echo "项目名称: ${PROJECT_NAME}"
echo "项目目录: ${PROJECT_DIR}"
echo "开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "======================================"

# 检查参考基因组索引和GTF文件
if [ ! -f "${GENOME_INDEX}.1.ht2" ] && [ ! -f "${GENOME_INDEX}.1.ht2l" ]; then
    echo "错误: HISAT2索引文件不存在: ${GENOME_INDEX}"
    exit 1
fi

if [ ! -f "${GTF_FILE}" ]; then
    echo "错误: GTF注释文件不存在: ${GTF_FILE}"
    exit 1
fi

echo "参考基因组检查通过"

# ===================== 创建输出目录 =====================

mkdir -p ${RESULTS_DIR}

echo "结果目录创建完成: ${RESULTS_DIR}"

# ===================== 步骤1: 获取样本列表 =====================

echo ""
echo "======================================"
echo "步骤1: 识别样本"
echo "======================================"

cd ${RAW_DATA_DIR}

# 检测配对端测序数据
SAMPLES=($(ls *_1.fastq.gz 2>/dev/null | sed 's/_1.fastq.gz//' || ls *_1.fq.gz 2>/dev/null | sed 's/_1.fq.gz//'))

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "错误: 在 ${RAW_DATA_DIR} 中未找到任何 *_1.fastq.gz 或 *_1.fq.gz 文件"
    exit 1
fi

echo "找到 ${#SAMPLES[@]} 个配对端样本"
echo "样本列表（前10个）:"
printf '%s\n' "${SAMPLES[@]:0:10}"

# ===================== 步骤2-5: 处理每个样本 =====================

SAMPLE_COUNT=0
COMPLETED_COUNT=0

for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    echo ""
    echo "======================================"
    echo "处理样本 ${SAMPLE_COUNT}/${#SAMPLES[@]}: ${SAMPLE}"
    echo "时间: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "======================================"
    
    # 样本输出目录
    SAMPLE_OUTDIR="${RESULTS_DIR}/${SAMPLE}"
    mkdir -p ${SAMPLE_OUTDIR}
    
    # 文件路径定义
    # 支持 .fastq.gz 和 .fq.gz 两种后缀
    if [ -f "${RAW_DATA_DIR}/${SAMPLE}_1.fastq.gz" ]; then
        R1_RAW="${RAW_DATA_DIR}/${SAMPLE}_1.fastq.gz"
        R2_RAW="${RAW_DATA_DIR}/${SAMPLE}_2.fastq.gz"
    else
        R1_RAW="${RAW_DATA_DIR}/${SAMPLE}_1.fq.gz"
        R2_RAW="${RAW_DATA_DIR}/${SAMPLE}_2.fq.gz"
    fi
    
    R1_TRIMMED="${SAMPLE_OUTDIR}/${SAMPLE}_trimmed_1.fq.gz"
    R2_TRIMMED="${SAMPLE_OUTDIR}/${SAMPLE}_trimmed_2.fq.gz"
    BAM_FILE="${SAMPLE_OUTDIR}/${SAMPLE}_aligned_sorted.bam"
    TRANSCRIPTS_GTF="${SAMPLE_OUTDIR}/${SAMPLE}_transcripts.gtf"
    GENE_ABUND="${SAMPLE_OUTDIR}/${SAMPLE}_gene_abund.tab"
    COUNTS_FILE="${SAMPLE_OUTDIR}/${SAMPLE}_counts.txt"
    
    # ---------- 步骤2: fastp质控 ----------
    if [ ! -f "${R1_TRIMMED}" ] || [ ! -f "${R2_TRIMMED}" ]; then
        echo "[${SAMPLE}] 步骤2: fastp质控和过滤..."
        
        ${FASTP_BIN} \
            -i ${R1_RAW} -I ${R2_RAW} \
            -q 20 -n 15 -u 50 -l 50 -e 20 \
            --thread ${THREADS} \
            -o ${R1_TRIMMED} \
            -O ${R2_TRIMMED} \
            -j ${SAMPLE_OUTDIR}/${SAMPLE}_fastp.json \
            -h ${SAMPLE_OUTDIR}/${SAMPLE}_fastp.html \
            2>&1 | tee ${SAMPLE_OUTDIR}/${SAMPLE}_fastp.log
        
        echo "[${SAMPLE}] 质控完成"
    else
        echo "[${SAMPLE}] 质控文件已存在，跳过"
    fi
    
    # ---------- 步骤3: HISAT2比对 ----------
    if [ ! -f "${BAM_FILE}" ]; then
        echo "[${SAMPLE}] 步骤3: HISAT2比对..."
        
        ${HISAT2_BIN} \
            -x ${GENOME_INDEX} \
            -p ${THREADS} \
            -1 ${R1_TRIMMED} \
            -2 ${R2_TRIMMED} \
            --dta \
            2> ${SAMPLE_OUTDIR}/${SAMPLE}_hisat2.log \
            | ${SAMTOOLS_BIN} view -bS - \
            | ${SAMTOOLS_BIN} sort -@ ${THREADS} \
            -o ${BAM_FILE}
        
        # 验证BAM文件生成
        if [ ! -f "${BAM_FILE}" ]; then
            echo "错误: BAM文件未生成 - ${BAM_FILE}"
            echo "请检查HISAT2日志: ${SAMPLE_OUTDIR}/${SAMPLE}_hisat2.log"
            exit 1
        fi
        
        # 索引BAM文件
        ${SAMTOOLS_BIN} index ${BAM_FILE}
        
        echo "[${SAMPLE}] 比对完成"
    else
        echo "[${SAMPLE}] BAM文件已存在，跳过"
        # 确保索引文件存在
        if [ ! -f "${BAM_FILE}.bai" ]; then
            ${SAMTOOLS_BIN} index ${BAM_FILE}
        fi
    fi
    
    # ---------- 步骤4: StringTie定量 ----------
    if [ ! -f "${GENE_ABUND}" ]; then
        echo "[${SAMPLE}] 步骤4: StringTie转录本组装和定量..."
        
        ${STRINGTIE_BIN} \
            ${BAM_FILE} \
            -G ${GTF_FILE} \
            -o ${TRANSCRIPTS_GTF} \
            -p ${THREADS} \
            -B \
            -e \
            -A ${GENE_ABUND} \
            2>&1 | tee ${SAMPLE_OUTDIR}/${SAMPLE}_stringtie.log
        
        # 验证输出文件生成
        if [ ! -f "${GENE_ABUND}" ]; then
            echo "错误: StringTie丰度文件未生成 - ${GENE_ABUND}"
            exit 1
        fi
        
        echo "[${SAMPLE}] StringTie定量完成"
    else
        echo "[${SAMPLE}] StringTie丰度文件已存在，跳过"
    fi
    
    # ---------- 步骤5: featureCounts定量 ----------
    if [ ! -f "${COUNTS_FILE}" ]; then
        echo "[${SAMPLE}] 步骤5: featureCounts基因计数..."
        
        ${FEATURECOUNTS_BIN} \
            -a ${GTF_FILE} \
            -o ${COUNTS_FILE} \
            -T ${THREADS} \
            -p \
            -t exon \
            -g gene_id \
            ${BAM_FILE} \
            2>&1 | tee ${SAMPLE_OUTDIR}/${SAMPLE}_featureCounts.log
        
        # 验证计数文件生成
        if [ ! -f "${COUNTS_FILE}" ]; then
            echo "错误: 计数文件未生成 - ${COUNTS_FILE}"
            exit 1
        fi
        
        echo "[${SAMPLE}] featureCounts定量完成"
    else
        echo "[${SAMPLE}] 计数文件已存在，跳过"
    fi
    
    # ---------- 清理中间文件 ----------
    echo "[${SAMPLE}] 清理中间文件..."
    rm -f ${R1_TRIMMED} ${R2_TRIMMED} ${BAM_FILE} ${BAM_FILE}.bai
    
    COMPLETED_COUNT=$((COMPLETED_COUNT + 1))
    echo "[${SAMPLE}] 处理完成！已完成样本: ${COMPLETED_COUNT}/${#SAMPLES[@]}"
done

# ===================== 步骤6: 合并featureCounts表达矩阵 =====================

echo ""
echo "======================================"
echo "步骤6: 合并featureCounts表达矩阵"
echo "======================================"

MERGED_COUNTS="${RESULTS_DIR}/all_samples_counts_matrix.txt"

# 提取第一个样本的基因信息列
FIRST_SAMPLE="${SAMPLES[0]}"
FIRST_COUNTS="${RESULTS_DIR}/${FIRST_SAMPLE}/${FIRST_SAMPLE}_counts.txt"

echo "提取基因信息..."
grep -v "^#" ${FIRST_COUNTS} | cut -f1-6 > ${MERGED_COUNTS}.tmp

# 提取每个样本的计数列
echo "提取各样本计数数据..."
for SAMPLE in "${SAMPLES[@]}"; do
    COUNTS_FILE="${RESULTS_DIR}/${SAMPLE}/${SAMPLE}_counts.txt"
    echo "  处理: ${SAMPLE}"
    grep -v "^#" ${COUNTS_FILE} | tail -n +2 | awk '{print $NF}' > ${RESULTS_DIR}/${SAMPLE}.counts.tmp
done

# 合并所有列
echo "合并所有数据..."
paste ${MERGED_COUNTS}.tmp ${RESULTS_DIR}/*.counts.tmp > ${MERGED_COUNTS}

# 清理临时文件
rm -f ${MERGED_COUNTS}.tmp ${RESULTS_DIR}/*.counts.tmp

echo "featureCounts表达矩阵生成完成: ${MERGED_COUNTS}"

# ===================== 步骤7: 合并StringTie TPM矩阵 =====================

echo ""
echo "======================================"
echo "步骤7: 合并StringTie TPM表达矩阵"
echo "======================================"

MERGED_TPM="${RESULTS_DIR}/all_samples_TPM_matrix.txt"
MERGED_FPKM="${RESULTS_DIR}/all_samples_FPKM_matrix.txt"
MERGED_COVERAGE="${RESULTS_DIR}/all_samples_Coverage_matrix.txt"

# 创建临时目录
TEMP_DIR="${RESULTS_DIR}/temp_merge"
mkdir -p ${TEMP_DIR}

echo "提取各样本TPM、FPKM和Coverage数据..."

# 提取第一个样本的基因ID和名称作为基础
FIRST_ABUND="${RESULTS_DIR}/${FIRST_SAMPLE}/${FIRST_SAMPLE}_gene_abund.tab"

# 提取基因信息（Gene ID和Gene Name）
tail -n +2 ${FIRST_ABUND} | awk -F'\t' '{print $1"\t"$2}' > ${TEMP_DIR}/gene_info.txt

# 提取每个样本的TPM、FPKM和Coverage值
for SAMPLE in "${SAMPLES[@]}"; do
    ABUND_FILE="${RESULTS_DIR}/${SAMPLE}/${SAMPLE}_gene_abund.tab"
    echo "  处理: ${SAMPLE}"
    
    # 提取TPM值（第9列）
    tail -n +2 ${ABUND_FILE} | awk -F'\t' '{print $9}' > ${TEMP_DIR}/${SAMPLE}.tpm.tmp
    
    # 提取FPKM值（第8列）
    tail -n +2 ${ABUND_FILE} | awk -F'\t' '{print $8}' > ${TEMP_DIR}/${SAMPLE}.fpkm.tmp
    
    # 提取Coverage值（第7列）
    tail -n +2 ${ABUND_FILE} | awk -F'\t' '{print $7}' > ${TEMP_DIR}/${SAMPLE}.cov.tmp
done

# 创建TPM矩阵表头
echo -ne "Gene_ID\tGene_Name" > ${MERGED_TPM}
for SAMPLE in "${SAMPLES[@]}"; do
    echo -ne "\t${SAMPLE}" >> ${MERGED_TPM}
done
echo "" >> ${MERGED_TPM}

# 创建FPKM矩阵表头
echo -ne "Gene_ID\tGene_Name" > ${MERGED_FPKM}
for SAMPLE in "${SAMPLES[@]}"; do
    echo -ne "\t${SAMPLE}" >> ${MERGED_FPKM}
done
echo "" >> ${MERGED_FPKM}

# 创建Coverage矩阵表头
echo -ne "Gene_ID\tGene_Name" > ${MERGED_COVERAGE}
for SAMPLE in "${SAMPLES[@]}"; do
    echo -ne "\t${SAMPLE}" >> ${MERGED_COVERAGE}
done
echo "" >> ${MERGED_COVERAGE}

# 合并TPM数据
echo "合并TPM数据..."
paste ${TEMP_DIR}/gene_info.txt ${TEMP_DIR}/*.tpm.tmp >> ${MERGED_TPM}

# 合并FPKM数据
echo "合并FPKM数据..."
paste ${TEMP_DIR}/gene_info.txt ${TEMP_DIR}/*.fpkm.tmp >> ${MERGED_FPKM}

# 合并Coverage数据
echo "合并Coverage数据..."
paste ${TEMP_DIR}/gene_info.txt ${TEMP_DIR}/*.cov.tmp >> ${MERGED_COVERAGE}

# 清理临时文件
rm -rf ${TEMP_DIR}

echo "TPM表达矩阵生成完成: ${MERGED_TPM}"
echo "FPKM表达矩阵生成完成: ${MERGED_FPKM}"
echo "Coverage矩阵生成完成: ${MERGED_COVERAGE}"

# ===================== 生成统计报告 =====================

echo ""
echo "======================================"
echo "分析流程完成！"
echo "完成时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "======================================"
echo ""
echo "输出文件位置："
echo "  各样本结果: ${RESULTS_DIR}/<sample_name>/"
echo "  Read counts矩阵: ${MERGED_COUNTS}"
echo "  TPM表达矩阵: ${MERGED_TPM}"
echo "  FPKM表达矩阵: ${MERGED_FPKM}"
echo "  Coverage矩阵: ${MERGED_COVERAGE}"
echo ""
echo "样本统计："
echo "  总样本数: ${#SAMPLES[@]}"
echo "  成功处理: $(ls ${RESULTS_DIR}/*/`basename ${SAMPLES[0]}`_counts.txt 2>/dev/null | wc -l) 个"
echo ""

# 生成比对率统计
echo "比对率统计（前20个样本）："
echo "----------------------------------------"
printf "%-25s %s\n" "样本名" "总体比对率"
echo "----------------------------------------"
for SAMPLE in "${SAMPLES[@]:0:20}"; do
    LOG_FILE="${RESULTS_DIR}/${SAMPLE}/${SAMPLE}_hisat2.log"
    if [ -f "${LOG_FILE}" ]; then
        ALIGN_RATE=$(grep "overall alignment rate" ${LOG_FILE} | awk '{print $1}')
        printf "%-25s %s\n" "${SAMPLE}" "${ALIGN_RATE}"
    fi
done
echo "----------------------------------------"

echo ""
echo "查看完整结果："
echo "  Read counts矩阵: head ${MERGED_COUNTS}"
echo "  TPM表达矩阵: head ${MERGED_TPM}"
echo "  FPKM表达矩阵: head ${MERGED_FPKM}"
echo "  样本质控报告: firefox ${RESULTS_DIR}/<sample>/<sample>_fastp.html"
echo "  比对日志: cat ${RESULTS_DIR}/<sample>/<sample>_hisat2.log"
echo "  StringTie日志: cat ${RESULTS_DIR}/<sample>/<sample>_stringtie.log"
echo "  featureCounts日志: cat ${RESULTS_DIR}/<sample>/<sample>_featureCounts.log"
echo ""
echo "======================================"
echo "所有任务完成！"
echo "======================================"
