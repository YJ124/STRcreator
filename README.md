# STRcreator
STRcreator is an ultrafast, high-fidelity sequencing simulator for Short Tandem Repeats (STRs). Using a decoupled architecture and geometric decay model, it accurately simulates PCR stutter and optical noise. It achieves near-zero JS divergence. ideal for benchmarking somatic mosaicism and MSI detection.
## 1.1 Introduction
Accurate simulation of Short Tandem Repeats (STRs) is critical for evaluating downstream variant calling algorithms. A key application of STRcreator is generating ground-truth datasets for studying **Somatic Mosaicism**. Somatic mosaicism, the presence of genetically distinct cell lineages within a single organism, often manifests as low-frequency mutations in highly mutable regions like STRs. 

Traditional simulators (e.g., ART, DWGSIM) strictly model optical errors, completely omitting the biological PCR polymerase slippage (stutter) inherent to STRs. Advanced profile-based simulators (e.g., simTR) model these artifacts via Hidden Markov Models but suffer from massive computational bottlenecks. STRcreator bridges this gap by decoupling biological sequence evolution from optical noise, delivering state-of-the-art biological fidelity with unmatched computational efficiency, making it highly suitable for benchmarking Microsatellite Instability (MSI) detection and somatic variant callers.

## 🏗️ Architecture
STRcreator utilizes a decoupled three-tier architecture:
1. **Biological Evolution**: Simulates fundamental STR motif expansions/contractions based on a Continuous-Time Markov Chain (CTMC).
2. **PCR Biochemical Kinetics**: Accurately mimics polymerase slippage using a strict **geometric exponential decay penalty**, ensuring realistic `-1` contraction-biased stutter profiles.
3. **Optical Sequencing Errors**: Fragments the amplicon pool and introduces empirical, platform-specific optical noise (e.g., Illumina PE150).

## ⚙️ Installation

Clone the repository and install the minimal required dependencies:

```bash
git clone [https://github.com/YourUsername/STRcreator.git](https://github.com/YourUsername/STRcreator.git)
cd STRcreator
pip install numpy scipy matplotlib
```
STRcreator simulates data in sequential, highly-optimized steps. Below is an example of generating a 30x coverage Illumina dataset for STR loci on chromosome 19.

Step 1: Simulate PCR Stutter (Biochemical Stage)
Generate the amplicon pool with biologically realistic stutter noise. The parameters below enforce a strict geometric penalty to match empirical baseline mutation rates.

```Bash
python pcr_stutter_simulator.py \
    -r hg38_chr19.fa \
    -g chr19_ground_truth.tsv \
    -o STRcreator_amplicons.fa \
    -c 28 \
    --peff 0.85 \
    --pdown 0.0001 \
    --pup 0.00002 \
    --pgeom 0.95 \
    -f 200 \
    -d 500
```
Step 2: Simulate Illumina Sequencing (Optical Stage)
Fragment the generated amplicons and introduce optical sequencing errors to create paired-end FASTQ files.

```Bash
python illumina_sequencer.py \
    -i STRcreator_amplicons.fa \
    -o STRcreator_chr19_30x \
    -c 30 \
    -l 150 \
    -e 0.001
```

Step 3: Downstream Alignment
The resulting synthetic reads (STRcreator_chr19_30x_R1.fq and R2.fq) are ready for standard alignment pipelines:

```Bash
bwa mem -t 8 hg38_chr19.fa STRcreator_chr19_30x_R1.fq STRcreator_chr19_30x_R2.fq | \
samtools view -Sb - | samtools sort -@ 8 -o STRcreator_chr19_30x.bam
samtools index STRcreator_chr19_30x.bam
```
