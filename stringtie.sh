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

# ==================== 第4步: 完美合并TPM矩阵 ====================
echo "[步骤 4/4] 合并TPM和FPKM矩阵..."
echo ""

# 创建完美的合并脚本
cat > ${WORK_DIR}/results/merge_expression_perfect.py << 'PYTHON_SCRIPT'
import sys
import pandas as pd
from pathlib import Path

stringtie_dir = sys.argv[1]
output_tpm = sys.argv[2]
output_fpkm = sys.argv[3]

print("开始合并表达量矩阵...")

# 获取所有样本目录
sample_dirs = sorted([d for d in Path(stringtie_dir).iterdir() if d.is_dir()])

if not sample_dirs:
    print("错误: 未找到样本目录")
    sys.exit(1)

print(f"找到 {len(sample_dirs)} 个样本")
print("")

# 读取所有样本
tpm_list = []
fpkm_list = []

for i, sample_dir in enumerate(sample_dirs):
    abundance_file = list(sample_dir.glob("*_abundance.txt"))
    
    if not abundance_file:
        print(f"  警告: {sample_dir.name} 没有abundance文件，跳过")
        continue
    
    abundance_file = abundance_file[0]
    df = pd.read_csv(abundance_file, sep='\t')
    sample_name = sample_dir.name
    
    # 提取TPM和FPKM
    tpm_df = df[['Gene ID', 'TPM']].copy()
    tpm_df.columns = ['Gene_ID', sample_name]
    
    fpkm_df = df[['Gene ID', 'FPKM']].copy()
    fpkm_df.columns = ['Gene_ID', sample_name]
    
    tpm_list.append(tpm_df)
    fpkm_list.append(fpkm_df)
    
    print(f"  样本 {i+1}/{len(sample_dirs)}: {sample_name} - {len(df)} 基因")

print("")
print("合并所有样本...")

# 使用outer join逐步合并（保留所有基因）
tpm_merged = tpm_list[0]
for tpm_df in tpm_list[1:]:
    tpm_merged = pd.merge(tpm_merged, tpm_df, on='Gene_ID', how='outer')

fpkm_merged = fpkm_list[0]
for fpkm_df in fpkm_list[1:]:
    fpkm_merged = pd.merge(fpkm_merged, fpkm_df, on='Gene_ID', how='outer')

# 填充缺失值为0
tpm_merged = tpm_merged.fillna(0)
fpkm_merged = fpkm_merged.fillna(0)

# 保存结果
tpm_merged.to_csv(output_tpm, sep='\t', index=False)
fpkm_merged.to_csv(output_fpkm, sep='\t', index=False)

print("")
print("✓ 合并完成！")
print(f"  TPM矩阵: {output_tpm}")
print(f"    维度: {tpm_merged.shape[0]} 基因 × {tpm_merged.shape[1]-1} 样本")
print(f"  FPKM矩阵: {output_fpkm}")
print(f"    维度: {fpkm_merged.shape[0]} 基因 × {fpkm_merged.shape[1]-1} 样本")
print("")
print("说明: 使用outer join合并，不同样本的基因数差异已自动处理（缺失值填充为0）")
PYTHON_SCRIPT

# 运行合并脚本
python3 ${WORK_DIR}/results/merge_expression_perfect.py \
    ${STRINGTIE_DIR} \
    ${WORK_DIR}/results/merged_TPM_matrix.txt \
    ${WORK_DIR}/results/merged_FPKM_matrix.txt

echo ""

# ==================== 完成 ====================
echo "================================"
echo "✓ 全部完成！"
echo "================================"
echo ""
echo "输出文件："
echo "  📊 TPM矩阵:  results/merged_TPM_matrix.txt"
echo "  📊 FPKM矩阵: results/merged_FPKM_matrix.txt"
echo "  📁 单个样本: results/stringtie/{sample}/{sample}_abundance.txt"
echo ""
echo "验证结果："
head -5 ${WORK_DIR}/results/merged_TPM_matrix.txt
echo ""
echo "运行时间: $(date)"
