#!/bin/bash

# StringTie完美重新运行脚本 - 适用于所有sheep项目
# 功能：
# 1. 删除所有旧的StringTie结果
# 2. 过滤GTF文件（移除gene行）
# 3. 运行StringTie
# 4. 完美合并TPM矩阵（使用outer join，处理基因数不一致问题）

set -e  # 遇到错误即退出

# ==================== 配置 ====================
WORK_DIR=$(pwd)
PROJECT_NAME=$(basename ${WORK_DIR})
STRINGTIE_DIR="${WORK_DIR}/results/stringtie"
BAM_DIR="${WORK_DIR}/results/star"
GTF_ORIGINAL="${WORK_DIR}/results/star_index/annotation.gtf"
GTF_FILTERED="${WORK_DIR}/results/star_index/annotation_stringtie.gtf"
LOG_DIR="${WORK_DIR}/results/logs"
THREADS=8

echo "================================"
echo "StringTie 完美重新运行脚本"
echo "项目: ${PROJECT_NAME}"
echo "================================"
echo ""

# ==================== 第1步: 删除所有旧的StringTie结果 ====================
echo "[步骤 1/4] 删除旧的StringTie结果..."

if [ -d "${STRINGTIE_DIR}" ]; then
    echo "  删除目录: ${STRINGTIE_DIR}"
    rm -rf ${STRINGTIE_DIR}
fi

# 删除旧的合并结果
rm -f ${WORK_DIR}/results/merged_TPM_matrix.txt 2>/dev/null
rm -f ${WORK_DIR}/results/merge_tpm.py 2>/dev/null

# 删除StringTie相关日志
rm -f ${LOG_DIR}/*_stringtie.log 2>/dev/null

echo "  ✓ 旧文件已清理"
echo ""

# ==================== 第2步: 准备过滤后的GTF ====================
echo "[步骤 2/4] 准备StringTie专用GTF（过滤gene行）..."

if [ ! -f "${GTF_ORIGINAL}" ]; then
    echo "  错误: 找不到GTF文件: ${GTF_ORIGINAL}"
    exit 1
fi

# 过滤掉gene行，避免StringTie报错
echo "  过滤 gene 类型的行..."
grep -v $'\tgene\t' ${GTF_ORIGINAL} > ${GTF_FILTERED}

echo "  原始GTF行数: $(wc -l < ${GTF_ORIGINAL})"
echo "  过滤后行数: $(wc -l < ${GTF_FILTERED})"
echo "  ✓ GTF准备完成: ${GTF_FILTERED}"
echo ""

# ==================== 第3步: 运行StringTie ====================
echo "[步骤 3/4] 运行StringTie..."
echo ""

# 获取样本列表
SAMPLES=($(ls ${WORK_DIR}/*_1.fastq.gz 2>/dev/null | xargs -n1 basename | sed 's/_1.fastq.gz//'))

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "  错误: 未找到样本文件"
    exit 1
fi

echo "  检测到 ${#SAMPLES[@]} 个样本"
echo ""

# 创建StringTie输出目录
mkdir -p ${STRINGTIE_DIR}
mkdir -p ${LOG_DIR}

# 处理每个样本
for sample in "${SAMPLES[@]}"; do
    echo "  [$(date '+%Y-%m-%d %H:%M:%S')] 处理: ${sample}"
    
    # 检查BAM文件是否存在
    BAM_FILE="${BAM_DIR}/${sample}/Aligned.sortedByCoord.out.bam"
    if [ ! -f "${BAM_FILE}" ]; then
        echo "    警告: BAM文件不存在，跳过: ${BAM_FILE}"
        continue
    fi
    
    # 创建样本输出目录
    mkdir -p ${STRINGTIE_DIR}/${sample}
    
    # 运行StringTie
    stringtie \
        ${BAM_FILE} \
        -G ${GTF_FILTERED} \
        -o ${STRINGTIE_DIR}/${sample}/${sample}.gtf \
        -A ${STRINGTIE_DIR}/${sample}/${sample}_abundance.txt \
        -p ${THREADS} \
        -e \
        2>&1 | tee ${LOG_DIR}/${sample}_stringtie.log
    
    # 检查是否成功
    if [ $? -eq 0 ] && [ -f "${STRINGTIE_DIR}/${sample}/${sample}_abundance.txt" ]; then
        echo "    ✓ 完成"
    else
        echo "    ✗ 失败"
    fi
done

echo ""
echo "  ✓ StringTie运行完成"
echo ""
