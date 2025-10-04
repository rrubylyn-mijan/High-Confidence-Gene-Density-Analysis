# High-Confidence-Gene-Density-Analysis
This workflow counts and calculates the genome-wide density of high-confidence (HC) protein-coding genes in wheat genome assemblies. Outputs are formatted for downstream visualization in tools such as Circos.

## 1. Count High-Confidence Genes
```bash
grep -c -w "gene" wheat.high.gff3
# Total = 104620
```

## 2. Extract Gene Coordinates
```bash
# Extract HC genes (chromosome, start, end)
awk '$3 == "gene"' wheat.high.gff3 \
    | awk '{print $1, $4, $5}' OFS="\t" > high-confidence-genes-wheat.bed
```

## 3. Prepare Genome File
```bash
# Extract chromosome sizes
grep "^##sequence-region" wheat.high.gff3 \
    | awk '{print $2, $4}' OFS="\t" > genome-wheat-file.txt

# Keep only chr* contigs present in the GFF
grep "^chr" high-confidence-genes-wheat.bed > filtered-high-confidence-genes-wheat.bed
awk '{print $1}' filtered-high-confidence-genes-wheat.bed | sort | uniq > bed-contigs-wheat.txt
grep -Ff bed-contigs-wheat.txt genome-wheat-file.txt > filtered-genome-wheat-file.txt
```

## 4. Calculate Gene Density in 1Mb Windows
```bash
bedtools makewindows -g filtered-genome-wheat-file.txt -w 1000000 \
  | bedtools intersect -a - -b filtered-high-confidence-genes-wheat.bed -c \
  > gene-density-per-bin-wheat.bed
```

## 5. Standardize Chromosome Names (Circos compatibility)
```bash
# Change chr* → ta* format
sed 's/^chr/ta/' gene-density-per-bin-wheat.bed > ta-gene-density-per-bin-wheat.bed

# Apply same renaming to ncRNA density if needed
sed 's/^chr/ta/' ncrna-density-wheat.bed > x-ncrna-density-wheat.bed
```

## 6. Identify Minimum and Maximum Density
```bash
awk 'NR==1 {min=$4; max=$4} 
     $4 < min {min=$4} 
     $4 > max {max=$4} 
     END {print "Minimum Density:", min; print "Maximum Density:", max}' \
     x-ncrna-density-wheat.bed
```

## 7. Frequency Distribution of Gene Density
```bash
awk '{counts[$NF]++} END {for (val in counts) print val, counts[val]}' \
  x-ncrna-density-wheat.bed | sort -n
```
