import argparse
import numpy as np
import csv
from collections import defaultdict

def read_fasta(fasta_file):
    print(f"Loading reference genome {fasta_file}...")
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

def simulate_pcr_cycle(molecule_pool, p_eff, p_down, p_up, p_geom):
    """
    终极修复版：加入极高衰减因子的几何分布，严格限制滑移步长
    """
    new_pool = molecule_pool.copy() 
    p_stay = 1.0 - p_down - p_up
    
    for copies, count in molecule_pool.items():
        if count <= 0: continue
        
        amplified_count = np.random.binomial(count, p_eff)
        
        if amplified_count > 0:
            stutter_results = np.random.multinomial(amplified_count, [p_down, p_up, p_stay])
            down_count, up_count, stay_count = stutter_results
            
            # 收缩 (Contraction) - 引入陡峭的几何衰减
            if down_count > 0:
                steps = np.arange(1, 6) 
                probs = p_geom * (1 - p_geom)**(steps - 1)
                probs = probs / probs.sum() 
                step_counts = np.random.multinomial(down_count, probs)
                
                for i, step in enumerate(steps):
                    if step_counts[i] > 0:
                        new_copies = max(1, copies - step) 
                        new_pool[new_copies] = new_pool.get(new_copies, 0) + step_counts[i]
                        
            # 扩张 (Expansion)
            if up_count > 0:
                steps = np.arange(1, 6)
                probs = p_geom * (1 - p_geom)**(steps - 1)
                probs = probs / probs.sum()
                step_counts = np.random.multinomial(up_count, probs)
                
                for i, step in enumerate(steps):
                    if step_counts[i] > 0:
                        new_pool[copies + step] = new_pool.get(copies + step, 0) + step_counts[i]
            
            # 正常复制 (No stutter)
            new_pool[copies] = new_pool.get(copies, 0) + stay_count
            
    return new_pool

def run_pcr_stutter_simulation(genome, gt_file, out_fasta, cycles, p_eff, p_down, p_up, p_geom, flank_len, total_reads):
    print(f"Simulating PCR with extremely steep geometric decay (p_geom={p_geom})...")
    
    with open(gt_file, 'r') as f, open(out_fasta, 'w') as out:
        reader = csv.DictReader(f, delimiter='\t')
        for row_idx, row in enumerate(reader):
            chrom = row['Chromosome']
            start = int(row['Start']) - 1  
            end = int(row['End'])          
            motif = row['Motif']
            bio_copies = int(row['Mutated_Copies']) 
            
            if 'Ref_Copies' in row:
                ref_copies = int(row['Ref_Copies'])
            else:
                ref_copies = round((end - start) / len(motif))
            
            if chrom not in genome: continue
                
            left_flank = genome[chrom][max(0, start - flank_len):start]
            right_flank = genome[chrom][end:min(len(genome[chrom]), end + flank_len)]
            ref_middle = genome[chrom][start:end]
            
            initial_molecules = 100 
            pool = {bio_copies: initial_molecules}
            
            for _ in range(cycles):
                pool = simulate_pcr_cycle(pool, p_eff, p_down, p_up, p_geom)
            
            total_amplicons = sum(pool.values())
            reads_written = 0
            for copies, count in pool.items():
                proportion = count / total_amplicons
                num_reads_for_this_variant = int(np.round(proportion * total_reads))
                
                if num_reads_for_this_variant > 0:
                    delta_copies = copies - ref_copies
                    if delta_copies > 0:
                        mutated_str_seq = ref_middle + (motif * delta_copies)
                    elif delta_copies < 0:
                        trim_length = abs(delta_copies) * len(motif)
                        mutated_str_seq = ref_middle[:-trim_length] if trim_length < len(ref_middle) else ""
                    else:
                        mutated_str_seq = ref_middle
                    
                    amplicon_seq = left_flank + mutated_str_seq + right_flank
                    for _ in range(num_reads_for_this_variant):
                        reads_written += 1
                        header = f">STR_Locus_{row_idx+1}_{chrom}_{start}_BioCopies_{bio_copies}_StutterCopies_{copies}_Read_{reads_written}"
                        out.write(f"{header}\n{amplicon_seq}\n")

    print(f"PCR Stutter simulation complete! Amplicon pool saved to {out_fasta}")

def main():
    parser = argparse.ArgumentParser(description="Simulate PCR Stutter with Strict Geometric Penalty")
    parser.add_argument("-r", "--reference", required=True)
    parser.add_argument("-g", "--ground_truth", required=True)
    parser.add_argument("-o", "--out_fasta", required=True)
    parser.add_argument("-c", "--cycles", type=int, default=28)
    parser.add_argument("--peff", type=float, default=0.85)
    
    # 🌟 核心参数下调：基础打滑率砍掉一个数量级
    parser.add_argument("--pdown", type=float, default=0.0001, help="Contraction rate (default: 0.0001)")
    parser.add_argument("--pup", type=float, default=0.00002, help="Expansion rate (default: 0.00002)")
    
    # 🌟 新增参数：几何衰减因子，越大越陡峭 (限制在 -1)
    parser.add_argument("--pgeom", type=float, default=0.95, help="Geometric decay penalty (default: 0.95)")
    
    parser.add_argument("-f", "--flank_len", type=int, default=200)
    parser.add_argument("-d", "--depth", type=int, default=500)
    
    args = parser.parse_args()
    
    genome_dict = read_fasta(args.reference)
    run_pcr_stutter_simulation(genome_dict, args.ground_truth, args.out_fasta, args.cycles, args.peff, args.pdown, args.pup, args.pgeom, args.flank_len, args.depth)

if __name__ == "__main__":
    main()