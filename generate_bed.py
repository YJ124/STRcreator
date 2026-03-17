import re
import argparse
import itertools

def generate_canonical_motifs(min_len=2, max_len=5):
    """
    生成所有可能的碱基组合 (去除单碱基重复如 AAAA)
    """
    bases = ['A', 'C', 'G', 'T']
    all_motifs = []
    
    for length in range(min_len, max_len + 1):
        # itertools.product 生成所有排列组合
        for p in itertools.product(bases, repeat=length):
            motif = "".join(p)
            # 过滤掉同聚物 (例如 'AA', 'CCCC')，因为我们主要关注微卫星STR
            if len(set(motif)) > 1:
                all_motifs.append(motif)
                
    return all_motifs

def get_min_repeats(motif_len):
    """根据 Motif 长度动态返回最低需要的连续重复次数"""
    if motif_len == 2:
        return 10  # 2bp 至少重复 10 次 (20bp)
    elif motif_len == 3:
        return 6   # 3bp 至少重复 6 次 (18bp)
    elif motif_len == 4:
        return 5   # 4bp 至少重复 5 次 (20bp)
    else:
        return 4   # 5bp及以上 至少重复 4 次

def find_all_strs(fasta_file, bed_file, min_motif_len=2, max_motif_len=5):
    print(f"Loading {fasta_file} to scan for all possible STRs...")
    chromosomes = {}
    
    with open(fasta_file, 'r') as f:
        curr_chr = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                curr_chr = line[1:].split()[0]
                chromosomes[curr_chr] = []
            else:
                chromosomes[curr_chr].append(line.upper())
                
    motifs_to_search = generate_canonical_motifs(min_len=min_motif_len, max_len=max_motif_len)
    print(f"Generated {len(motifs_to_search)} different motifs to search. Scanning...")
    
    bed_records = []
    
    for chrom, seq_list in chromosomes.items():
        full_seq = "".join(seq_list)
        
        for motif in motifs_to_search:
            min_rep = get_min_repeats(len(motif))
            pattern = f"({motif}){{{min_rep},}}"
            
            for match in re.finditer(pattern, full_seq):
                bed_records.append({
                    'chrom': chrom,
                    'start': match.start(),
                    'end': match.end(),
                    'motif': motif
                })
                
    # 按物理坐标排序 (非常重要！防止不同 Motif 的结果杂乱无章)
    print("Sorting bed records by genomic coordinates...")
    bed_records.sort(key=lambda x: (x['chrom'], x['start']))
    
    print(f"Writing to {bed_file}...")
    with open(bed_file, 'w') as out:
        for record in bed_records:
            out.write(f"{record['chrom']}\t{record['start']}\t{record['end']}\t{record['motif']}\n")
            
    print(f"Done! Successfully found {len(bed_records)} STR loci.")

def main():
    parser = argparse.ArgumentParser(description="Auto-generate comprehensive STR BED file")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output BED file")
    args = parser.parse_args()
    
    find_all_strs(args.input, args.output)

if __name__ == "__main__":
    main()