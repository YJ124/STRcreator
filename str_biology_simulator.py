import argparse
import numpy as np
from scipy.linalg import expm
import csv
import sys

def build_transition_rate_matrix(max_copies, alpha, beta):
    """
    构建逐步突变模型(SMM)的转移率矩阵 Q
    """
    Q = np.zeros((max_copies, max_copies))
    for i in range(max_copies):
        # 突变速率随长度非线性变化 (mu_i = alpha + beta * i)
        mu_i = alpha + beta * i
        
        # 假设扩张和收缩的基础概率对等（可根据实际测序数据调整权重）
        expansion_rate = mu_i 
        contraction_rate = mu_i 

        if i < max_copies - 1:
            Q[i, i + 1] = expansion_rate
        if i > 0:
            Q[i, i - 1] = contraction_rate
        
        # 对角线元素：保持不变的速率
        # q_{i,i} = -(mu_u + mu_d)
        if i == 0:
            Q[i, i] = -expansion_rate
        elif i == max_copies - 1:
            Q[i, i] = -contraction_rate
        else:
            Q[i, i] = -(expansion_rate + contraction_rate)
    return Q

def read_fasta(fasta_file):
    """读取参考基因组到内存字典中"""
    print(f"Loading reference genome from {fasta_file}...")
    genome = {}
    current_chrom = ""
    current_seq = []
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_chrom:
                    genome[current_chrom] = "".join(current_seq)
                current_chrom = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_chrom:
            genome[current_chrom] = "".join(current_seq)
    return genome

def simulate_str_mutations(genome, bed_file, time_t, Q_matrix, max_copies):
    """根据CTMC模型进行变异，并返回修改后的基因组和Ground Truth记录"""
    print("Simulating biological STR mutations...")
    P_matrix = expm(Q_matrix * time_t) # P(t) = exp(Qt)
    
    # 归一化处理，防止浮点数精度问题导致概率和不为1
    P_matrix = P_matrix / P_matrix.sum(axis=1, keepdims=True)
    
    ground_truth = []
    
    # 按染色体分组读取BED，并按位置逆序排序，防止插入/缺失导致的坐标偏移影响后续位点
    bed_regions = {}
    with open(bed_file, 'r') as f:
        for line in f:
            chrom, start, end, motif = line.strip().split('\t')
            if chrom not in bed_regions:
                bed_regions[chrom] = []
            bed_regions[chrom].append((int(start), int(end), motif))
            
    for chrom in bed_regions:
        bed_regions[chrom].sort(key=lambda x: x[0], reverse=True)
        
    modified_genome = {}
    
    for chrom, seq in genome.items():
        if chrom not in bed_regions:
            modified_genome[chrom] = seq
            continue
            
        seq_list = list(seq)
        for start, end, motif in bed_regions[chrom]:
            motif_len = len(motif)
            original_copies = (end - start) // motif_len
            
            # 防止拷贝数超出矩阵上限
            safe_orig_copies = min(original_copies, max_copies - 1)
            
            # 根据转移概率矩阵抽取新的拷贝数
            probabilities = P_matrix[safe_orig_copies, :]
            new_copies = np.random.choice(np.arange(max_copies), p=probabilities)
            
            # 计算差异
            delta_copies = new_copies - original_copies
            new_str_seq = motif * new_copies
            
            # 替换原始序列
            seq_list[start:end] = list(new_str_seq)
            
            # 记录 Ground Truth (存储标签)
            ground_truth.append({
                'Chromosome': chrom,
                'Start': start,
                'End': end,
                'Motif': motif,
                'Original_Copies': original_copies,
                'Mutated_Copies': new_copies,
                'Delta_Copies': delta_copies,
                'Biological_Variant_Type': 'Expansion' if delta_copies > 0 else ('Contraction' if delta_copies < 0 else 'None')
            })
            
        modified_genome[chrom] = "".join(seq_list)
        
    return modified_genome, ground_truth

def main():
    parser = argparse.ArgumentParser(description="Biological STR Mutator based on CTMC and SMM")
    parser.add_argument("-r", "--reference", required=True, help="Input reference FASTA file")
    parser.add_argument("-b", "--bed", required=True, help="Input target BED file")
    parser.add_argument("-o", "--out_fasta", required=True, help="Output mutated FASTA file")
    parser.add_argument("-g", "--ground_truth", required=True, help="Output Ground Truth TSV file")
    parser.add_argument("-t", "--time", type=float, default=1.0, help="Evolutionary time t (default: 1.0)")
    parser.add_argument("--alpha", type=float, default=0.001, help="Baseline mutation rate alpha (default: 0.001)")
    parser.add_argument("--beta", type=float, default=0.0005, help="Length effect coefficient beta (default: 0.0005)")
    parser.add_argument("--max_copies", type=int, default=150, help="Maximum STR copies limit for state matrix (default: 150)")
    
    args = parser.parse_args()
    
    # 1. 构建转移矩阵
    Q = build_transition_rate_matrix(args.max_copies, args.alpha, args.beta)
    
    # 2. 读取参考基因组
    genome = read_fasta(args.reference)
    
    # 3. 模拟变异过程
    mutated_genome, gt_records = simulate_str_mutations(genome, args.bed, args.time, Q, args.max_copies)
    
    # 4. 输出变异后的基因组序列
    print(f"Writing mutated genome to {args.out_fasta}...")
    with open(args.out_fasta, 'w') as f:
        for chrom, seq in mutated_genome.items():
            f.write(f">{chrom}\n")
            # 每80个字符换行，遵循标准Fasta格式
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")
                
    # 5. 输出 Ground Truth 标签文件
    print(f"Writing ground truth labels to {args.ground_truth}...")
    with open(args.ground_truth, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['Chromosome', 'Start', 'End', 'Motif', 'Original_Copies', 'Mutated_Copies', 'Delta_Copies', 'Biological_Variant_Type'], delimiter='\t')
        writer.writeheader()
        writer.writerows(gt_records)
        
    print("Simulation of biological scene completed successfully!")

if __name__ == "__main__":
    main()