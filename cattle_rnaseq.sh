#!/bin/bash

# Cattle RNA-Seq 上游分析流程 - 使用NCBI参考基因组
# 参考原有流程：Trimmomatic质控 -> STAR比对 -> StringTie + featureCounts定量
# 使用方法：直接放在项目文件夹（如 PRJNA1112374）中运行

set -e  # 遇到错误即退出

# ==================== 配置参数 ====================
# 当前目录作为工作目录（项目文件夹）
WORK_DIR=$(pwd)
PROJECT_NAME=$(basename ${WORK_DIR})

# 参考基因组目录（上一级目录）
REF_DIR="../genome_reference"

# 输出目录（在当前项目文件夹下创建results目录）
RESULTS_DIR="${WORK_DIR}/results"
TRIMMOMATIC_DIR="${RESULTS_DIR}/trimmomatic"
STAR_INDEX_DIR="${RESULTS_DIR}/star_index"
BAM_DIR="${RESULTS_DIR}/star"
COUNT_DIR="${RESULTS_DIR}/featureCounts"
STRINGTIE_DIR="${RESULTS_DIR}/stringtie"
LOG_DIR="${RESULTS_DIR}/logs"

# NCBI参考基因组文件
GENOME_FA="${REF_DIR}/GCF_002263795.1_ARS-UCD1.2_genomic.fna"
GTF_FILE="${REF_DIR}/GCF_002263795.1_ARS-UCD1.2_genomic.gtf"

# Trimmomatic 参数（按照参考代码）
# 设置 trimmomatic 环境变量
if [ -z "${trimmomatic}" ]; then
    # 尝试自动查找 Trimmomatic
    TRIMMOMATIC_SHARE=$(find $(conda info --base 2>/dev/null)/envs/*/share -maxdepth 1 -name "trimmomatic-*" -type d 2>/dev/null | head -1)
    if [ -n "${TRIMMOMATIC_SHARE}" ] && [ -f "${TRIMMOMATIC_SHARE}/trimmomatic.jar" ]; then
        export trimmomatic=${TRIMMOMATIC_SHARE}
        echo "自动检测到Trimmomatic: ${trimmomatic}"
    else
        # 尝试默认路径
        DEFAULT_PATH="/data/mzwang/.conda/envs/cattle/share/trimmomatic-0.39-2"
        if [ -f "${DEFAULT_PATH}/trimmomatic.jar" ]; then
            export trimmomatic=${DEFAULT_PATH}
            echo "使用默认Trimmomatic路径: ${trimmomatic}"
        fi
    fi
fi

# 线程数
THREADS=8

# ==================== 检查环境 ====================
echo "================================"
echo "Cattle RNA-Seq 分析流程 (NCBI参考基因组)"
echo "Project: ${PROJECT_NAME}"
echo "================================"

# 检查参考基因组
if [ ! -f "${GENOME_FA}" ]; then
    echo "错误: 找不到参考基因组文件！"
    echo "期望路径: ${GENOME_FA}"
    exit 1
fi

if [ ! -f "${GTF_FILE}" ]; then
    echo "错误: 找不到GTF注释文件！"
    echo "期望路径: ${GTF_FILE}"
    exit 1
fi

echo "参考基因组: ${GENOME_FA}"
echo "GTF注释: ${GTF_FILE}"

# 检查Trimmomatic
if [ -z "${trimmomatic}" ] || [ ! -f "${trimmomatic}/trimmomatic.jar" ]; then
    echo "错误: Trimmomatic未正确配置"
    echo "请设置环境变量: export trimmomatic=/path/to/trimmomatic"
    echo "或检查路径: ${trimmomatic}"
    exit 1
fi

if [ ! -f "${trimmomatic}/adapters/TruSeq3-PE.fa" ]; then
    echo "错误: 找不到adapter文件"
    echo "期望路径: ${trimmomatic}/adapters/TruSeq3-PE.fa"
    exit 1
fi

echo "Trimmomatic: ${trimmomatic}/trimmomatic.jar"
echo "Adapters: ${trimmomatic}/adapters/"

# 检查其他软件
for cmd in STAR stringtie featureCounts samtools java; do
    if ! command -v ${cmd} &> /dev/null; then
        echo "错误: ${cmd} 未安装"
        exit 1
    fi
done

echo "✓ 所有必需软件已就绪"
echo ""

# ==================== 创建输出目录 ====================
echo "创建输出目录..."
mkdir -p ${TRIMMOMATIC_DIR} ${STAR_INDEX_DIR} ${BAM_DIR} ${COUNT_DIR} ${STRINGTIE_DIR} ${LOG_DIR}

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

# ==================== 1. Trimmomatic 质控 ====================
echo "================================"
echo "步骤 1: 质量控制 (Trimmomatic)"
echo "================================"

for sample in "${SAMPLES[@]}"; do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] 处理样本: ${sample}"
    
    # 使用jar文件运行Trimmomatic
    java -jar ${trimmomatic}/trimmomatic.jar PE -phred33 \
        ${WORK_DIR}/${sample}_1.fastq.gz \
        ${WORK_DIR}/${sample}_2.fastq.gz \
        ${TRIMMOMATIC_DIR}/${sample}_1.clean.fq.gz \
        ${TRIMMOMATIC_DIR}/${sample}_1_unpaired.fastq.gz \
        ${TRIMMOMATIC_DIR}/${sample}_2.clean.fq.gz \
        ${TRIMMOMATIC_DIR}/${sample}_2_unpaired.fastq.gz \
        -threads ${THREADS} \
        ILLUMINACLIP:${trimmomatic}/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
        2>&1 | tee ${LOG_DIR}/${sample}_trimmomatic.log
    
    echo "完成: ${sample}"
done

echo "Trimmomatic 质控完成！"
echo ""

# ==================== 2. 构建 STAR 索引 ====================
echo "================================"
echo "步骤 2: 构建 STAR 索引"
echo "================================"

# 构建STAR索引（如果不存在）
if [ ! -f "${STAR_INDEX_DIR}/SA" ]; then
    echo "构建 STAR 索引（这可能需要30-60分钟）..."
    STAR \
        --runMode genomeGenerate \
        --runThreadN ${THREADS} \
        --genomeDir ${STAR_INDEX_DIR} \
        --genomeFastaFiles ${GENOME_FA} \
        --sjdbGTFfile ${GTF_FILE} \
        --sjdbOverhang 99 \
        2>&1 | tee ${LOG_DIR}/star_build.log
    echo "STAR 索引构建完成！"
else
    echo "STAR 索引已存在，跳过构建步骤"
fi
echo ""

# ==================== 3. STAR 比对 ====================
echo "================================"
echo "步骤 3: STAR 比对 (Two-pass mode)"
echo "================================"

for sample in "${SAMPLES[@]}"; do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] 比对样本: ${sample}"
    
    # 创建样本专用输出目录
    SAMPLE_OUT_DIR="${BAM_DIR}/${sample}"
    mkdir -p ${SAMPLE_OUT_DIR}
    
    # STAR比对（按照参考代码的参数）
    STAR \
        --runThreadN ${THREADS} \
        --genomeDir ${STAR_INDEX_DIR} \
        --sjdbGTFfile ${GTF_FILE} \
        --quantMode TranscriptomeSAM \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --readFilesCommand zcat \
        --outFilterMismatchNmax 3 \
        --readFilesIn ${TRIMMOMATIC_DIR}/${sample}_1.clean.fq.gz ${TRIMMOMATIC_DIR}/${sample}_2.clean.fq.gz \
        --outFileNamePrefix ${SAMPLE_OUT_DIR}/${sample}- \
        --twopassMode Basic \
        2>&1 | tee ${LOG_DIR}/${sample}_star.log
    
    # 建立BAM索引
    echo "为 ${sample} 建立索引..."
    samtools index ${SAMPLE_OUT_DIR}/${sample}-Aligned.sortedByCoord.out.bam
    
    echo "完成: ${sample}"
done

echo "STAR 比对完成！"
echo ""

# ==================== 4. StringTie TPM/FPKM 定量 ====================
echo "================================"
echo "步骤 4: TPM/FPKM 定量 (StringTie)"
echo "================================"

# 创建过滤后的GTF（移除gene行，避免StringTie错误）
echo "准备StringTie专用GTF（过滤gene行）..."
grep -v $'\tgene\t' ${GTF_FILE} > ${STAR_INDEX_DIR}/annotation_stringtie.gtf

for sample in "${SAMPLES[@]}"; do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] StringTie 处理: ${sample}"
    
    mkdir -p ${STRINGTIE_DIR}/${sample}
    
    # 按照参考代码的参数运行StringTie，使用过滤后的GTF
    stringtie \
        -p ${THREADS} \
        -e \
        -B \
        -G ${STAR_INDEX_DIR}/annotation_stringtie.gtf \
        -o ${STRINGTIE_DIR}/${sample}/${sample}.gtf \
        -A ${STRINGTIE_DIR}/${sample}/${sample}.tsv \
        ${BAM_DIR}/${sample}/${sample}-Aligned.sortedByCoord.out.bam \
        2>&1 | tee ${LOG_DIR}/${sample}_stringtie.log
    
    echo "完成: ${sample}"
done

echo "StringTie 定量完成！"
echo ""

# ==================== 5. featureCounts 定量 ====================
echo "================================"
echo "步骤 5: 基因计数 (featureCounts)"
echo "================================"

# 准备所有 BAM 文件列表
BAM_FILES=$(ls ${BAM_DIR}/*/${*}-Aligned.sortedByCoord.out.bam | tr '\n' ' ')

echo "运行 featureCounts..."
featureCounts \
    -T ${THREADS} \
    -p \
    -t exon \
    -g gene_id \
    -a ${GTF_FILE} \
    -o ${COUNT_DIR}/all_samples_counts.txt \
    ${BAM_FILES} \
    2>&1 | tee ${LOG_DIR}/featureCounts.log

echo "featureCounts 完成！"

# 生成简化版计数矩阵
echo "生成计数矩阵..."
cat ${COUNT_DIR}/all_samples_counts.txt | \
    grep -v '^#' > ${COUNT_DIR}/counts_matrix.txt

echo "计数矩阵已保存: ${COUNT_DIR}/counts_matrix.txt"
echo ""

# ==================== 6. 合并 StringTie 表达量 ====================
echo "================================"
echo "步骤 6: 合并 StringTie 表达量矩阵"
echo "================================"

# 创建Python脚本来合并表达量
cat > ${RESULTS_DIR}/merge_expression.py << 'PYTHON_SCRIPT'
import sys
import pandas as pd
from pathlib import Path

stringtie_dir = sys.argv[1]
output_tpm = sys.argv[2]
output_fpkm = sys.argv[3]

# 获取所有样本目录
sample_dirs = sorted([d for d in Path(stringtie_dir).iterdir() if d.is_dir()])

if not sample_dirs:
    print("错误: 未找到样本目录")
    sys.exit(1)

# 读取第一个样本
first_file = list(sample_dirs[0].glob("*.tsv"))[0]
merged = pd.read_csv(first_file, sep='\t')
sample_name = sample_dirs[0].name

# 提取基因信息
gene_info = merged[['Gene ID', 'Gene Name', 'Reference', 'Strand', 'Start', 'End', 'Coverage']].copy()

# 创建TPM和FPKM数据框
tpm_data = pd.DataFrame({'Gene_ID': merged['Gene ID']})
fpkm_data = pd.DataFrame({'Gene_ID': merged['Gene ID']})

# 添加第一个样本
tpm_data[sample_name] = merged['TPM'].values
fpkm_data[sample_name] = merged['FPKM'].values

# 合并其他样本
for sample_dir in sample_dirs[1:]:
    tsv_file = list(sample_dir.glob("*.tsv"))[0]
    df = pd.read_csv(tsv_file, sep='\t')
    sample_name = sample_dir.name
    tpm_data[sample_name] = df['TPM'].values
    fpkm_data[sample_name] = df['FPKM'].values

# 保存结果
tpm_data.to_csv(output_tpm, sep='\t', index=False)
fpkm_data.to_csv(output_fpkm, sep='\t', index=False)

print(f"TPM矩阵已保存: {output_tpm}")
print(f"FPKM矩阵已保存: {output_fpkm}")
print(f"维度: {tpm_data.shape[0]} 基因 × {tpm_data.shape[1]-1} 样本")
PYTHON_SCRIPT

# 运行Python脚本合并表达量
python3 ${RESULTS_DIR}/merge_expression.py \
    ${STRINGTIE_DIR} \
    ${RESULTS_DIR}/merged_TPM_matrix.txt \
    ${RESULTS_DIR}/merged_FPKM_matrix.txt

echo ""

# ==================== 7. 生成比对统计 ====================
echo "================================"
echo "步骤 7: 生成比对统计"
echo "================================"

for sample in "${SAMPLES[@]}"; do
    samtools flagstat ${BAM_DIR}/${sample}/${sample}-Aligned.sortedByCoord.out.bam \
        > ${LOG_DIR}/${sample}_flagstat.txt
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
echo "  ⭐ Clean reads: ${TRIMMOMATIC_DIR}/*_clean.fq.gz"
echo "  ⭐ 比对文件: ${BAM_DIR}/*/*-Aligned.sortedByCoord.out.bam"
echo "  ⭐ Counts矩阵: ${COUNT_DIR}/counts_matrix.txt"
echo "  ⭐ TPM矩阵: ${RESULTS_DIR}/merged_TPM_matrix.txt"
echo "  ⭐ FPKM矩阵: ${RESULTS_DIR}/merged_FPKM_matrix.txt"
echo "  ⭐ 日志文件: ${LOG_DIR}/"
echo ""
echo "参考基因组信息:"
echo "  - 来源: NCBI"
echo "  - 版本: ARS-UCD1.2"
echo "  - 染色体命名: NC_XXXXXX (NCBI格式)"
echo ""
echo "下一步分析:"
echo "  - 查看 MultiQC 报告检查整体质量"
echo "  - 使用 counts_matrix.txt 进行差异表达分析 (DESeq2/edgeR)"
echo "  - 使用 TPM/FPKM 矩阵进行表达量比较"
echo ""
echo "运行时间: $(date)"
