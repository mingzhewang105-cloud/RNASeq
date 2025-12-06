"""
简化版RNA-seq分析流程 - 直接在项目文件夹运行
放置位置: ~/sheep/PRJDB29864/Snakefile (或其他项目文件夹)
"""

import glob

# ============ 配置区域 - 根据需要修改 ============
REFERENCE_DIR = "../reference_genome"
GENOME_FA = f"{REFERENCE_DIR}/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna.gz"
GENOME_GTF = f"{REFERENCE_DIR}/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.gtf.gz"
THREADS = 8  # 根据服务器配置调整

# ============ 自动识别样本 ============
SAMPLES = []
for f in glob.glob("*_1.fastq.gz"):
    sample = f.replace("_1.fastq.gz", "")
    SAMPLES.append(sample)

print(f"检测到 {len(SAMPLES)} 个样本: {SAMPLES}")

# ============ 最终输出 ============
rule all:
    input:
        expand("results/featureCounts/{sample}_counts.txt", sample=SAMPLES),
        expand("results/stringtie/{sample}_abundance.txt", sample=SAMPLES),
        "results/multiqc_report.html"

# ============ 步骤1: 质量控制 (fastp) ============
rule fastp:
    input:
        r1 = "{sample}_1.fastq.gz",
        r2 = "{sample}_2.fastq.gz"
    output:
        r1 = "results/fastp/{sample}_1.clean.fastq.gz",
        r2 = "results/fastp/{sample}_2.clean.fastq.gz",
        html = "results/fastp/{sample}_fastp.html",
        json = "results/fastp/{sample}_fastp.json"
    log:
        "results/logs/fastp/{sample}.log"
    threads: 4
    shell:
        """
        mkdir -p results/fastp results/logs/fastp
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              -f 3 -t 3 -l 36 -r -W 4 -M 15 \
              -h {output.html} -j {output.json} \
              -w {threads} 2> {log}
        """

# ============ 步骤2: 建立STAR索引 (只运行一次) ============
rule star_index:
    input:
        fasta = GENOME_FA,
        gtf = GENOME_GTF
    output:
        dir = directory("results/star_index"),
        done = "results/star_index/genome_index.done"
    log:
        "results/logs/star_index.log"
    threads: THREADS
    shell:
        """
        mkdir -p {output.dir} results/logs
        
        echo "解压参考基因组文件..."
        gunzip -c {input.fasta} > {output.dir}/genome.fa
        gunzip -c {input.gtf} > {output.dir}/annotation.gtf
        
        echo "开始建立STAR索引..."
        STAR --runMode genomeGenerate \
             --genomeDir {output.dir} \
             --genomeFastaFiles {output.dir}/genome.fa \
             --sjdbGTFfile {output.dir}/annotation.gtf \
             --runThreadN {threads} \
             --sjdbOverhang 99 2> {log}
        
        touch {output.done}
        echo "STAR索引建立完成！"
        """

# ============ 步骤3: STAR比对 ============
rule star_align:
    input:
        r1 = "results/fastp/{sample}_1.clean.fastq.gz",
        r2 = "results/fastp/{sample}_2.clean.fastq.gz",
        index_done = "results/star_index/genome_index.done"
    output:
        bam = "results/star/{sample}/Aligned.sortedByCoord.out.bam",
        log_final = "results/star/{sample}/Log.final.out"
    params:
        index = "results/star_index",
        prefix = "results/star/{sample}/"
    log:
        "results/logs/star/{sample}.log"
    threads: THREADS
    shell:
        """
        mkdir -p results/star/{wildcards.sample} results/logs/star
        
        STAR --genomeDir {params.index} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outFileNamePrefix {params.prefix} \
             --outSAMtype BAM SortedByCoordinate \
             --twopassMode Basic \
             --runThreadN {threads} \
             --outSAMunmapped Within \
             --outSAMattributes Standard 2> {log}
        """

# ============ 步骤4: 索引BAM文件 ============
rule index_bam:
    input:
        "results/star/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "results/star/{sample}/Aligned.sortedByCoord.out.bam.bai"
    log:
        "results/logs/samtools/{sample}.log"
    shell:
        """
        mkdir -p results/logs/samtools
        samtools index {input} 2> {log}
        """

# ============ 步骤5: 基因计数 (featureCounts) ============
rule featurecounts:
    input:
        bam = "results/star/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "results/star/{sample}/Aligned.sortedByCoord.out.bam.bai"
    output:
        counts = "results/featureCounts/{sample}_counts.txt",
        summary = "results/featureCounts/{sample}_counts.txt.summary"
    params:
        gtf = "results/star_index/annotation.gtf"
    log:
        "results/logs/featureCounts/{sample}.log"
    threads: 4
    shell:
        """
        mkdir -p results/featureCounts results/logs/featureCounts
        
        featureCounts -T {threads} \
                      -p -B -C \
                      -a {params.gtf} \
                      -o {output.counts} \
                      {input.bam} 2> {log}
        """

# ============ 步骤6: TPM定量 (StringTie) ============
rule stringtie:
    input:
        bam = "results/star/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "results/star/{sample}/Aligned.sortedByCoord.out.bam.bai"
    output:
        gtf = "results/stringtie/{sample}.gtf",
        abundance = "results/stringtie/{sample}_abundance.txt"
    params:
        gtf = "results/star_index/annotation.gtf"
    log:
        "results/logs/stringtie/{sample}.log"
    threads: 4
    shell:
        """
        mkdir -p results/stringtie results/logs/stringtie
        
        stringtie {input.bam} \
                  -G {params.gtf} \
                  -o {output.gtf} \
                  -A {output.abundance} \
                  -p {threads} \
                  -e 2> {log}
        """

# ============ 步骤7: 合并所有样本的counts文件 (新增!) ============
rule merge_counts:
    input:
        expand("results/featureCounts/{sample}_counts.txt", sample=SAMPLES)
    output:
        "results/merged_counts_matrix.txt"
    log:
        "results/logs/merge_counts.log"
    run:
        import pandas as pd
        import sys
        
        try:
            all_counts = []
            
            # 读取每个样本的counts文件
            for i, count_file in enumerate(input):
                sample_name = count_file.split('/')[-1].replace('_counts.txt', '')
                
                # 读取文件，跳过注释行
                df = pd.read_csv(count_file, sep='\t', comment='#')
                
                # 第一个文件保留基因信息列
                if i == 0:
                    merged = df[['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']].copy()
                
                # 提取counts列并重命名
                counts_col = df.columns[-1]  # 最后一列是counts
                merged[sample_name] = df[counts_col]
            
            # 保存合并结果
            merged.to_csv(output[0], sep='\t', index=False)
            
            print(f"成功合并 {len(input)} 个样本的counts文件")
            print(f"合并矩阵大小: {merged.shape[0]} 基因 × {merged.shape[1]-6} 样本")
            print(f"输出文件: {output[0]}")
            
        except Exception as e:
            print(f"错误: {e}", file=sys.stderr)
            raise

# ============ 步骤7: 质量报告 (MultiQC) ============
rule multiqc:
    input:
        expand("results/fastp/{sample}_fastp.json", sample=SAMPLES),
        expand("results/star/{sample}/Log.final.out", sample=SAMPLES),
        expand("results/featureCounts/{sample}_counts.txt.summary", sample=SAMPLES)
    output:
        "results/multiqc_report.html"
    log:
        "results/logs/multiqc.log"
    shell:
        """
        multiqc results/ -o results/ -n multiqc_report.html --force 2> {log}
        """
