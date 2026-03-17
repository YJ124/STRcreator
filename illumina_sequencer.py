import argparse
import random
import numpy as np
from collections import defaultdict

def reverse_complement(seq):
    mapping = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')
    return seq.translate(mapping)[::-1]

def introduce_sequencing_errors(seq, base_error_rate):
    bases = ['A', 'T', 'C', 'G']
    mutated_seq = []
    quals = []
    
    for i, base in enumerate(seq):
        current_error_rate = base_error_rate * (1 + (i / len(seq)) * 9) 
        if base.upper() in bases and random.random() < current_error_rate:
            mutated_seq.append(random.choice([b for b in bases if b != base.upper()]))
            quals.append(',')  
        else:
            mutated_seq.append(base.upper())
            quals.append('F')  
            
    return "".join(mutated_seq), "".join(quals)

def simulate_paired_end_sequencing(amplicon_fasta, out_prefix, coverage, read_len, error_rate):
    print(f"Loading amplicon pool from {amplicon_fasta}...")
    locus_amplicons = defaultdict(list)
    
    with open(amplicon_fasta, 'r') as f:
        current_header = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_header = line[1:]
            else:
                parts = current_header.split('_')
                locus_id = f"{parts[0]}_{parts[1]}_{parts[2]}" 
                locus_amplicons[locus_id].append((current_header, line))
                
    out_r1 = f"{out_prefix}_R1.fq"
    out_r2 = f"{out_prefix}_R2.fq"
    
    print(f"Simulating WGS-like fragmentation and Illumina PE{read_len} sequencing at ~{coverage}x coverage...")
    
    reads_generated = 0
    with open(out_r1, 'w') as f1, open(out_r2, 'w') as f2:
        for locus, amplicons in locus_amplicons.items():
            if not amplicons: continue
            
            # 动态计算需要多少对 Read 才能使整个扩增子达到目标覆盖度
            avg_amp_len = np.mean([len(seq) for _, seq in amplicons])
            num_pairs_needed = int(np.ceil((coverage * avg_amp_len) / (2 * read_len)))
            
            # 从扩增子池中有放回抽样，完美保留生物学突变和 Stutter 的比例分布
            sampled_amplicons = random.choices(amplicons, k=num_pairs_needed)
            
            for idx, (header, seq) in enumerate(sampled_amplicons):
                seq_len = len(seq)
                if seq_len < read_len: continue
                
                # 模拟真实建库：正态分布的插入片段大小 (均值 350, 标准差 50)
                insert_size = int(np.random.normal(350, 50))
                # 约束片段大小不超过序列总长，也不短于 read_len
                insert_size = max(read_len, min(insert_size, seq_len)) 
                
                # 在扩增子上随机选择打断的起始位点
                max_start = seq_len - insert_size
                start_pos = random.randint(0, max_start) if max_start > 0 else 0
                
                # 截取物理片段
                fragment = seq[start_pos : start_pos + insert_size]
                
                # 双端测序
                r1_seq, r1_qual = introduce_sequencing_errors(fragment[:read_len], error_rate)
                r2_seq, r2_qual = introduce_sequencing_errors(reverse_complement(fragment[-read_len:]), error_rate)
                
                read_name = f"@{header}_SimRead_{idx+1}"
                f1.write(f"{read_name}/1\n{r1_seq}\n+\n{r1_qual}\n")
                f2.write(f"{read_name}/2\n{r2_seq}\n+\n{r2_qual}\n")
                reads_generated += 1
                
    print(f"Success! Generated {reads_generated} paired-end read pairs.")

def main():
    parser = argparse.ArgumentParser(description="Simulate Illumina Sequencing")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--out_prefix", required=True)
    parser.add_argument("-c", "--coverage", type=int, default=50)
    parser.add_argument("-l", "--read_len", type=int, default=150)
    parser.add_argument("-e", "--error_rate", type=float, default=0.001)
    
    args = parser.parse_args()
    simulate_paired_end_sequencing(args.input, args.out_prefix, args.coverage, args.read_len, args.error_rate)

if __name__ == "__main__":
    main()