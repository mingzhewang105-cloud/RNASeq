#!/bin/bash

# RNA-Seq 上游分析流程
# 包括：质控 -> 比对 -> 定量 (featureCounts + StringTie)
# 使用方法：直接放在项目文件夹（如 PRJDB29864）中运行

set -e  # 遇到错误即退出

# ==================== 配置参数 ====================
# 当前目录作为工作目录（项目文件夹）
WORK_DIR=$(pwd)
PROJECT_NAME=$(basename ${WORK_DIR})

# 参考基因组目录（上一级目录）
REF_DIR="../reference_genome"

# 输出目录（在当前项目文件夹下创建results目录）
RESULTS_DIR="${WORK_DIR}/results"
FASTP_DIR="${RESULTS_DIR}/fastp"
STAR_INDEX_DIR="${RESULTS_DIR}/star_index"
BAM_DIR="${RESULTS_DIR}/star"
COUNT_DIR="${RESULTS_DIR}/featureCounts"
STRINGTIE_DIR="${RESULTS_DIR}/stringtie"
LOG_DIR="${RESULTS_DIR}/logs"

# 参考基因组文件
GENOME_FA="${REF_DIR}/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna.gz"
GTF_FILE="${REF_DIR}/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.gtf.gz"

# fastp 参数（按照参考文档）
FASTP_PARAMS="-f 3 -t 3 -l 36 -r -W 4 -M 15"

# 线程数
THREADS=8

# ==================== 检查参考基因组 ====================
echo "================================"
echo "RNA-Seq Analysis Pipeline"
echo "Project: ${PROJECT_NAME}"
echo "================================"

if [ ! -f "${GENOME_FA}" ] || [ ! -f "${GTF_FILE}" ]; then
    echo "错误: 找不到参考基因组文件！"
    echo "请确认以下文件存在:"
    echo "  ${GENOME_FA}"
    echo "  ${GTF_FILE}"
    exit 1
fi

echo "参考基因组: ${GENOME_FA}"
echo "GTF注释文件: ${GTF_FILE}"

# ==================== 创建输出目录 ====================
echo ""
echo "创建输出目录..."
mkdir -p ${FASTP_DIR} ${STAR_INDEX_DIR} ${BAM_DIR} ${COUNT_DIR} ${STRINGTIE_DIR} ${LOG_DIR}

# ==================== 获取样本列表 ====================
SAMPLES=($(ls ${WORK_DIR}/*_1.fastq.gz 2>/dev/null | xargs -n1 basename | sed 's/_1.fastq.gz//'))

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "错误: 未找到 *_1.fastq.gz 文件！"
    echo "请确认当前目录中有配对的 fastq.gz 文件"
    exit 1
fi

echo "检测到 ${#SAMPLES[@]} 个样本"
echo "样本列表: ${SAMPLES[@]}"
echo ""

# ==================== 1. fastp 质控 ====================
echo "================================"
echo "步骤 1: 质量控制 (fastp)"
echo "================================"

for sample in "${SAMPLES[@]}"; do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] 处理样本: ${sample}"
    
    fastp \
        -i ${WORK_DIR}/${sample}_1.fastq.gz \
        -I ${WORK_DIR}/${sample}_2.fastq.gz \
        -o ${FASTP_DIR}/${sample}_1.clean.fastq.gz \
        -O ${FASTP_DIR}/${sample}_2.clean.fastq.gz \
        ${FASTP_PARAMS} \
        -h ${FASTP_DIR}/${sample}_fastp.html \
        -j ${FASTP_DIR}/${sample}_fastp.json \
        --thread ${THREADS} \
        2>&1 | tee ${LOG_DIR}/${sample}_fastp.log
    
    echo "完成: ${sample}"
done

echo "fastp 质控完成！"
echo ""

# ==================== 2. 构建 STAR 索引 ====================
echo "================================"
echo "步骤 2: 构建 STAR 索引"
echo "================================"

# 解压基因组和GTF文件（如果需要）
if [ ! -f "${STAR_INDEX_DIR}/genome.fa" ]; then
    echo "解压基因组文件..."
    gunzip -c ${GENOME_FA} > ${STAR_INDEX_DIR}/genome.fa
fi

if [ ! -f "${STAR_INDEX_DIR}/annotation.gtf" ]; then
    echo "解压GTF文件..."
    gunzip -c ${GTF_FILE} > ${STAR_INDEX_DIR}/annotation.gtf
fi

# 构建STAR索引（如果不存在）
if [ ! -f "${STAR_INDEX_DIR}/SA" ]; then
    echo "构建 STAR 索引（这可能需要30-60分钟）..."
    STAR \
        --runMode genomeGenerate \
        --runThreadN ${THREADS} \
        --genomeDir ${STAR_INDEX_DIR} \
        --genomeFastaFiles ${STAR_INDEX_DIR}/genome.fa \
        --sjdbGTFfile ${STAR_INDEX_DIR}/annotation.gtf \
        --sjdbOverhang 99 \
        2>&1 | tee ${LOG_DIR}/star_build.log
    echo "STAR 索引构建完成！"
else
    echo "STAR 索引已存在，跳过构建步骤"
fi
echo ""

# ==================== 3. STAR 比对（Two-pass mode）====================
echo "================================"
echo "步骤 3: STAR 比对 (Two-pass mode)"
echo "================================"

for sample in "${SAMPLES[@]}"; do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] 比对样本: ${sample}"
    
    # 创建样本专用输出目录
    SAMPLE_OUT_DIR="${BAM_DIR}/${sample}"
    mkdir -p ${SAMPLE_OUT_DIR}
    
    STAR \
        --genomeDir ${STAR_INDEX_DIR} \
        --readFilesIn ${FASTP_DIR}/${sample}_1.clean.fastq.gz ${FASTP_DIR}/${sample}_2.clean.fastq.gz \
        --outFileNamePrefix ${SAMPLE_OUT_DIR}/ \
        --runThreadN ${THREADS} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        2>&1 | tee ${LOG_DIR}/${sample}_star.log
    
    # 建立BAM索引
    echo "为 ${sample} 建立索引..."
    samtools index ${SAMPLE_OUT_DIR}/Aligned.sortedByCoord.out.bam
    
    echo "完成: ${sample}"
done

echo "STAR 比对完成！"
echo ""

# ==================== 4. featureCounts 定量 ====================
echo "================================"
echo "步骤 4: 基因计数 (featureCounts)"
echo "================================"

# 准备所有 BAM 文件列表
BAM_FILES=$(ls ${BAM_DIR}/*/Aligned.sortedByCoord.out.bam | tr '\n' ' ')

echo "运行 featureCounts..."
featureCounts \
    -T ${THREADS} \
    -p \
    -B \
    -C \
    -t exon \
    -g gene_id \
    -a ${STAR_INDEX_DIR}/annotation.gtf \
    -o ${COUNT_DIR}/all_samples_counts.txt \
    ${BAM_FILES} \
    2>&1 | tee ${LOG_DIR}/featureCounts.log

echo "featureCounts 完成！"

# 生成简化版计数矩阵（保留基因信息和计数）
echo "生成计数矩阵..."
cat ${COUNT_DIR}/all_samples_counts.txt | \
    grep -v '^#' > ${COUNT_DIR}/counts_matrix.txt

echo "计数矩阵已保存: ${COUNT_DIR}/counts_matrix.txt"
echo ""

# ==================== 5. StringTie TPM 定量 ====================
echo "================================"
echo "步骤 5: TPM 定量 (StringTie)"
echo "================================"

for sample in "${SAMPLES[@]}"; do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] StringTie 处理: ${sample}"
    
    mkdir -p ${STRINGTIE_DIR}/${sample}
    
    stringtie \
        ${BAM_DIR}/${sample}/Aligned.sortedByCoord.out.bam \
        -G ${STAR_INDEX_DIR}/annotation.gtf \
        -o ${STRINGTIE_DIR}/${sample}/${sample}.gtf \
        -A ${STRINGTIE_DIR}/${sample}/${sample}_abundance.txt \
        -p ${THREADS} \
        -e \
        2>&1 | tee ${LOG_DIR}/${sample}_stringtie.log
    
    echo "完成: ${sample}"
done

echo "StringTie 定量完成！"
echo ""

# ==================== 6. 合并 TPM 矩阵 ====================
echo "================================"
echo "步骤 6: 合并 TPM 矩阵"
echo "================================"

# 创建Python脚本来合并TPM值
cat > ${RESULTS_DIR}/merge_tpm.py << 'PYTHON_SCRIPT'
import sys
import pandas as pd
from pathlib import Path

stringtie_dir = sys.argv[1]
output_file = sys.argv[2]

# 获取所有样本目录
sample_dirs = sorted([d for d in Path(stringtie_dir).iterdir() if d.is_dir()])

if not sample_dirs:
    print("错误: 未找到样本目录")
    sys.exit(1)

# 读取第一个样本的abundance文件
first_file = list(sample_dirs[0].glob("*_abundance.txt"))[0]
merged = pd.read_csv(first_file, sep='\t')
sample_name = sample_dirs[0].name

# 保留基因信息列
gene_info = merged[['Gene ID', 'Gene Name', 'Reference', 'Strand', 'Start', 'End', 'Coverage', 'FPKM']]
tpm_data = pd.DataFrame({'Gene_ID': merged['Gene ID'], sample_name: merged['TPM']})

# 合并其他样本
for sample_dir in sample_dirs[1:]:
    abundance_file = list(sample_dir.glob("*_abundance.txt"))[0]
    df = pd.read_csv(abundance_file, sep='\t')
    sample_name = sample_dir.name
    tpm_data[sample_name] = df['TPM'].values

# 保存结果
tpm_data.to_csv(output_file, sep='\t', index=False)
print(f"TPM矩阵已保存: {output_file}")
print(f"维度: {tpm_data.shape[0]} 基因 × {tpm_data.shape[1]-1} 样本")
PYTHON_SCRIPT

# 运行Python脚本合并TPM
python3 ${RESULTS_DIR}/merge_tpm.py ${STRINGTIE_DIR} ${RESULTS_DIR}/merged_TPM_matrix.txt

echo ""

# ==================== 7. 生成比对统计 ====================
echo "================================"
echo "步骤 7: 生成比对统计"
echo "================================"

for sample in "${SAMPLES[@]}"; do
    samtools flagstat ${BAM_DIR}/${sample}/Aligned.sortedByCoord.out.bam > ${LOG_DIR}/${sample}_flagstat.txt
done

echo "统计完成！"
echo ""

# ==================== 8. MultiQC 汇总报告（如果安装了）====================
echo "================================"
echo "步骤 8: 生成 MultiQC 报告"
echo "================================"

if command -v multiqc &> /dev/null; then
    echo "运行 MultiQC..."
    multiqc ${RESULTS_DIR} -o ${RESULTS_DIR} -n multiqc_report.html --force
    echo "MultiQC 报告已生成: ${RESULTS_DIR}/multiqc_report.html"
else
    echo "未安装 MultiQC，跳过此步骤"
fi
echo ""

# ==================== 完成 ====================
echo "================================"
echo "流程完成！"
echo "================================"
echo "结果位置:"
echo "  📁 项目目录: ${WORK_DIR}"
echo "  📁 结果目录: ${RESULTS_DIR}"
echo ""
echo "重要输出文件:"
echo "  ⭐ 质控报告: ${FASTP_DIR}/*_fastp.html"
echo "  ⭐ 比对文件: ${BAM_DIR}/*/Aligned.sortedByCoord.out.bam"
echo "  ⭐ Counts矩阵: ${COUNT_DIR}/counts_matrix.txt"
echo "  ⭐ TPM矩阵: ${RESULTS_DIR}/merged_TPM_matrix.txt"
echo "  ⭐ 日志文件: ${LOG_DIR}/"
echo ""
echo "下一步分析:"
echo "  - 查看 MultiQC 报告检查整体质量"
echo "  - 使用 counts_matrix.txt 进行差异表达分析 (DESeq2/edgeR)"
echo "  - 使用 merged_TPM_matrix.txt 进行表达量比较"
echo ""
echo "运行时间: $(date)"
