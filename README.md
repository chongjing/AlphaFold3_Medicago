# AlphaFold3_Medicago
This project is to predict interactions between 222 nodulation regulators with 653 rhizobial secreted proteins of *Sinorhizobium meliloti* 2011.

## Table of Contents

- [00. AlphaFold3 set-up](https://github.com/chongjing/AlphaFold3_Medicago#00-alphafold3-set-up)
- [01. Secreted Protein Characterization of *Sinorhizobium meliloti* 2011](https://github.com/chongjing/AlphaFold3_Medicago#01-screted-protein-of-sinorhizobium-meliloti-2011)
  - [01.1 Putative effectors](https://github.com/chongjing/AlphaFold3_Medicago?tab=readme-ov-file#011-putative-effectors)
  - [01.2 Function Annotation](https://github.com/chongjing/AlphaFold3_Medicago?tab=readme-ov-file#012-function-annotation)
  - [01.3 Subcellular Localization](https://github.com/chongjing/AlphaFold3_Medicago?tab=readme-ov-file#013-subcellular-location)
  - [01.4 KEGG annotation](https://github.com/chongjing/AlphaFold3_Medicago?tab=readme-ov-file#014-kegg-annotation)
  - [01.5 Enrichment visualization](https://github.com/chongjing/AlphaFold3_Medicago?tab=readme-ov-file#015-enrichment-visualization)
- [02. Sinorhizobium RNAseq analysis](https://github.com/chongjing/AlphaFold3_Medicago?tab=readme-ov-file#02-sinorhizobiumrnaseq)
  - [02.1 Data Download and Trimming](https://github.com/chongjing/AlphaFold3_Medicago?tab=readme-ov-file#021-data-download-and-trimming)
  - [02.2 Mapping and TPM counts](https://github.com/chongjing/AlphaFold3_Medicago?tab=readme-ov-file#022-mapping-and-tpm-counts)
  - [02.3 Expression Cluster](https://github.com/chongjing/AlphaFold3_Medicago?tab=readme-ov-file#023-expression-cluster)
- [03. AlphaFold3 Prediction](https://github.com/chongjing/AlphaFold3_Medicago?tab=readme-ov-file#03-alphafold3-prediction)
  - [03.1 MSA](https://github.com/chongjing/AlphaFold3_Medicago?tab=readme-ov-file#031-msa)
  - [03.2 Inference](https://github.com/chongjing/AlphaFold3_Medicago?tab=readme-ov-file#032-inference)
  - [03.3 AF3 summary & visualization](https://github.com/chongjing/AlphaFold3_Medicago?tab=readme-ov-file#033-af3-summary--visualization)
  - [03.4 Distribution of number of targets](https://github.com/chongjing/AlphaFold3_Medicago?tab=readme-ov-file#034-distribution-of-number-of-targets)


## 00. AlphaFold3 set-up
```bash
# build singularity image from a docker 
singularity build AlphaFold3.0.1.sif docker://alanfwilliams/alphafold3:latest

## download genetic databases
## after searching, i find can download from nju mirror
## script provided by alphafold3 not works
wget https://mirror.nju.edu.cn/alphafold/data-3/bfd-first_non_consensus_sequences.fasta.zst --no-check-certificate
wget https://mirror.nju.edu.cn/alphafold/data-3/mgy_clusters_2022_05.fa.zst --no-check-certificate
wget https://mirror.nju.edu.cn/alphafold/data-3/nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta.zst --no-check-certificate
wget https://mirror.nju.edu.cn/alphafold/data-3/pdb_2022_09_28_mmcif_files.tar.zst --no-check-certificate
wget https://mirror.nju.edu.cn/alphafold/data-3/pdb_seqres_2022_09_28.fasta.zst --no-check-certificate
wget https://mirror.nju.edu.cn/alphafold/data-3/rnacentral_active_seq_id_90_cov_80_linclust.fasta.zst --no-check-certificate
wget https://mirror.nju.edu.cn/alphafold/data-3/rfam_14_9_clust_seq_id_90_cov_80_rep_seq.fasta.zst --no-check-certificate
wget https://mirror.nju.edu.cn/alphafold/data-3/uniprot_all_2021_04.fa.zst --no-check-certificate
wget https://mirror.nju.edu.cn/alphafold/data-3/uniref90_2022_05.fa.zst --no-check-certificate
## unzip databases
for NAME in mgy_clusters_2022_05.fa \
            bfd-first_non_consensus_sequences.fasta \
            uniref90_2022_05.fa uniprot_all_2021_04.fa \
            pdb_seqres_2022_09_28.fasta \
            rnacentral_active_seq_id_90_cov_80_linclust.fasta \
            nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta \
            rfam_14_9_clust_seq_id_90_cov_80_rep_seq.fasta ; do
  echo "Start Fetching '${NAME}'"
  zstd --decompress ${NAME}.zst > "./public_databases/${NAME}" 
  zstd -d pdb_2022_09_28_mmcif_files.tar.zst &
done
```

## 01. Screted Protein of *Sinorhizobium meliloti* 2011
### 01.1 Putative effectors
```bash
mamba activate python3.7.12
LD_LIBRARY_PATH=/data/pathology/program/Miniforge3/envs/python3.7.12/lib:$LD_LIBRARY_PATH

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/346/065/GCF_000346065.1_ASM34606v1/GCF_000346065.1_ASM34606v1_protein.faa.gz
unzip GCF_000346065.1_ASM34606v1_protein.faa.gz

/data/pathology/program/Miniforge3/envs/python3.7.12/bin/signalp6 --fastafile GCF_000346065.1_ASM34606v1_protein.faa --output_dir ./signalp  --format all --organism other --mode slow --model_dir /data/pathology/program/SignalP/signalp6_slow_sequential/signalp-6-package/models/
awk -F"\t" '$6 > "0.8" {print $1}' signalp/output.gff3 | awk -F" " '{print $1}' > 03.SP.list
seqkit grep --by-name --pattern-file 03.SP.list --use-regexp --out-file 03.SP.protein.faa GCF_000346065.1_ASM34606v1_protein.faa

/data/pathology/program/Miniforge3/envs/python3.7.12/bin/biolib run --local 'DTU/DeepTMHMM:1.0.24' --fasta 03.SP.protein.faa > 04.TMHMM
 ```

### 01.2 Function Annotation
```bash
# interproscan 6.0.0 was used to annotate putative function
#export JAVA_HOME=/rds/user/cx264/rds-csc_programmes-FTKWLWDeHys/programs/java/jdk-21.0.2/
#export PATH=$JAVA_HOME/bin:/rds/project/rds-FTKWLWDeHys/programs/mambaforge/bin:/rds/project/rds-FTKWLWDeHys/programs/mambaforge/bin:/rds/user/cx264/rds-csc_programmes-FTKWLWDeHys/programs/java/jdk-21.0.2/bin:$PATH
export JAVA_HOME=/home/cx264/program/jdk-25.0.2/
export PATH=/home/cx264/program/ruby-3.1.2/install/bin:$HOME/local/bin:/home/cx264/program/jdk-25.0.2/bin:$PATH
curl -s https://get.nextflow.io | bash
cd /home/cx264/rds/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/04.AlphaFold3/06.interproscan
nextflow run ebi-pf-team/interproscan6 \
      -profile singularity \
      --input GCF_000346065.1_ASM34606v1_protein.faa  \
      --datadir interproscan --goterms --pathways
```
### 01.3 Subcellular Localization
```bash
# prokaryotic protein using deeplocpro 1.0.0
git clone https://github.com/Jaimomar99/deeplocpro.git
cd deeplocpro && pip install .
/home/cx264/program/anaconda3/bin/deeplocpro --fasta 01.Rhizobium.faa --output 01.Rhizobium

# eukaryotic protein using deeploc 2.1
wget https://services.healthtech.dtu.dk/download/91b90ec3-2fea-41ae-8de2-9973684473a8/deeploc-2.1.All.tar.gz
tar -xvzf deeploc-2.1.All.tar.gz && cd deeploc2_package/ && pip install .
deeploc2 --fasta 02.Medicago.faa --output 02.Medicago
```
### 01.4 KEGG annotation
```bash
/home/cx264/rds/rds-csc_programmes-FTKWLWDeHys/programs/kofamScan/kofam_scan-1.3.0/exec_annotation -o GCF_000346065.1_ASM34606v1_protein.KoFamScan.txt ../GCF_000346065.1_ASM34606v1_protein.faa --profile=/home/cx264/rds/rds-csc_programmes-FTKWLWDeHys/programs/kofamScan/profiles/prokaryote.hal --cpu=2
```
### 01.5 Enrichment visualization
```R
# GO and KEGG enrichment were performed using TBtools, outputs were used for visualization.
# GO enrichment
library(ggplot2)
library(dplyr)
library(patchwork)

# Read tab-delimited GO enrichment data
df <- read.delim("01.GO.enrich.tsv", header = TRUE, stringsAsFactors = FALSE)
head(df)

# Process p-values: handle zeros for log transformation
df$p_adj <- as.numeric(df$p_adj)
df$p_adj[df$p_adj == 0] <- 1e-20  # Avoid log(0)
df$neg_log_pval <- -log10(df$p_adj)

df <- df %>%
  group_by(Class) %>%
  arrange(desc(EnrichmentScore)) %>%
  mutate(GO_Name = factor(GO_Name, levels = rev(GO_Name))) %>%  # Reverse for top-to-bottom ordering
  ungroup()

# Create main combined plot (all classes)
main_plot <- ggplot(df, aes(
  x = EnrichmentScore,
  y = GO_Name,
  size = HitsGenesCountsInSelectedSet,
  color = neg_log_pval
)) +
  geom_point(alpha = 0.85, stroke = 0.3) +
  scale_color_gradient(
    low = "#228B22",   # Green for low significance
    high = "#D73027",  # Red for high significance
    name = "-log10(p_adj)",
    limits = c(0, max(df$neg_log_pval, na.rm = TRUE))
  ) +
  scale_size(
    name = "Gene count",
    range = c(3, 10),
    breaks = pretty(range(df$HitsGenesCountsInSelectedSet), n = 4)
  ) +
  labs(x = "Enrichment score", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),                     # Remove ALL internal grid lines
    panel.border = element_rect(fill = NA, colour = "black", size = 0.5),  # Keep outer frame
    axis.text.y = element_text(size = 9, face = "italic"),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(5, 5, 5, 5)
  )

# Create class-specific subplots (BP, CC, MF
create_class_plot <- function(class_name, title) {
  subset_df <- df %>% filter(Class == class_name) %>%
    arrange(desc(EnrichmentScore)) %>%
    mutate(GO_Name = factor(GO_Name, levels = rev(GO_Name)))

  if (nrow(subset_df) == 0) return(NULL)

  ggplot(subset_df, aes(
    x = EnrichmentScore,
    y = GO_Name,
    size = HitsGenesCountsInSelectedSet,
    color = neg_log_pval
  )) +
    geom_point(alpha = 0.85, stroke = 0.3) +
    scale_color_gradient(
      low = "#228B22",
      high = "#D73027",
      limits = c(0, max(df$neg_log_pval, na.rm = TRUE))
    ) +
    scale_size(range = c(3, 10)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", size = 0.4),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      legend.position = "none"
    )
}

bp_plot <- create_class_plot("BP", "BP")
cc_plot <- create_class_plot("CC", "CC")
mf_plot <- create_class_plot("MF", "MF")

# Combine plots: main plot (left) + 3 class panels stacked (right)
combined_plot <- main_plot +
  plot_layout(widths = c(4, 1)) +
  inset_element(bp_plot, left = 1, bottom = 0.66, right = 1, top = 1) +
  inset_element(cc_plot, left = 1, bottom = 0.33, right = 1, top = 0.66) +
  inset_element(mf_plot, left = 1, bottom = 0, right = 1, top = 0.33)
combined_plot <- main_plot | (bp_plot / cc_plot / mf_plot) +
  plot_layout(widths = c(4, 1), guides = "collect")
ggsave("go_enrichment_bubble_mainplot.pdf",
       main_plot,
       width = 6,
       height = 9,
       dpi = 1200,
       bg = "white")
ggsave("go_enrichment_bubble_mainplot.svg",
       main_plot,
       width = 6,
       height = 9,
       dpi = 1200,
       bg = "white")
#install.packages("svglite", repos = "https://cloud.r-project.org/")
ggsave("go_enrichment_bubble_mainplot.jpeg",
       main_plot,
       width = 6,
       height = 9,
       dpi = 1200,
       bg = "white")
```
```R
# KEGG enrichment
library(ggplot2)
library(dplyr)
df <- read.delim("02.KEGG.enrich.tsv", header = TRUE, stringsAsFactors = FALSE)
df$p_adj <- as.numeric(df$p_adj)
df$p_adj[df$p_adj == 0] <- 1e-16
df$neg_log_pval <- -log10(df$p_adj)
df <- df %>%
  arrange(desc(enrichFactor)) %>%
  mutate(`TermName` = factor(`TermName`, levels = rev(`TermName`)))
kegg_plot <- ggplot(df, aes(
  x = enrichFactor,
  y = `TermName`,
  size = GeneHitsInSelectedSet,
  color = neg_log_pval
)) +
  geom_point(alpha = 0.85, stroke = 0.3) +
  scale_color_gradient(
    low = "#228B22",   # Green for low significance
    high = "#D73027",  # Red for high significance
    name = "-log10(p-adj)",
    limits = c(0, max(df$neg_log_pval, na.rm = TRUE))
  ) +
  scale_size(
    name = "Gene count",
    range = c(3, 12),
    breaks = pretty(range(df$GeneHitsInSelectedSet), n = 4)
  ) +
  scale_y_discrete(position = "right") +  # CRITICAL: Move y-axis labels to RIGHT
  labs(
    x = "Enrichment factor",
    y = NULL,
    title = "KEGG Pathway Enrichment"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),                     # Remove ALL internal grids
    panel.border = element_rect(fill = NA, colour = "black", size = 0.6),  # Outer frame only
    axis.text.y = element_text(size = 9, face = "italic", hjust = 0),      # Align right-side labels
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.margin = margin(10, 15, 10, 5)  # Extra right margin for labels
  )

ggsave("02.KEGG_enrichment_bubble_mainplot.pdf",
       kegg_plot,
       width = 6,
       height = 4,
       dpi = 1200,
       bg = "white")
ggsave("02.KEGG_enrichment_bubble_mainplot.jpeg",
       kegg_plot,
       width = 6,
       height = 4,
       dpi = 1200,
       bg = "white")
ggsave("02.KEGG_enrichment_bubble_mainplot.svg",
       kegg_plot,
       width = 6,
       height = 4,
       dpi = 1200,
       bg = "white")
```

GO and KEGG enrichment:
<table>
  <tr>
    <td><img src="https://github.com/chongjing/AlphaFold3_Medicago/blob/main/enrichment/go_enrichment_bubble_mainplot.jpeg" alt="Image 1" width="400"/></td>
    <td><img src="https://github.com/chongjing/AlphaFold3_Medicago/blob/main/enrichment/02.KEGG_enrichment_bubble_mainplot.jpeg" alt="Image 2" width="400"/></td>
  </tr>
</table>

## 02. Sinorhizobium.RNAseq
### 02.1 Data Download and Trimming

```bash
cd /data/pathology/cxia/projects/Giles/Jinpeng/05.Sinorhizobium.RNAseq/01.Data/

#SRR17176337 SRR17176355 SRR17176359: 4DPI
for i in {18299090..18299092} {18299142..18299149} 17176337 17176355 17176359; do
        echo "Processing SRR${i} "
        /data/pathology/cxia/program/SRAToolkit/sratoolkit.3.0.10-ubuntu64/bin/prefetch SRR${i} && /data/pathology/cxia/program/SRAToolkit/sratoolkit.3.0.10-ubuntu64/bin/fastq-dump --split-files --gzip SRR${i} && rm -rf ./SRR${i}
        output_fwd_paired="SRR${i}_1P.fq.gz"
        output_fwd_unpaired="SRR${i}_1U.fq.gz"
        output_rev_paired="SRR${i}_2P.fq.gz"
        output_rev_unpaired="SRR${i}_2U.fq.gz"
        java -jar /data/pathology/program/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 32 -summary "SRR${i}.summary" "SRR${i}_1.fastq.gz" "SRR${i}_2.fastq.gz" \
        "$output_fwd_paired" "$output_fwd_unpaired" "$output_rev_paired" "$output_rev_unpaired" \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:60
        rm -rf SRR${i}_1.fastq.gz SRR${i}_2.fastq.gz
        echo "Done SRR${i}"
done
```
### 02.2 Mapping and TPM counts 

```bash
output_dir="/data/pathology/cxia/projects/Giles/Jinpeng/05.Sinorhizobium.RNAseq/02.mapping/"

# Run STAR iterate through each sample
for i in {18299090..18299092} {18299142..18299147}; do
    cd "${output_dir}"

    # Extract sample name
    sample_name="SRR${i}"

    # Set directory
    mkdir -p "${sample_name}"
    echo "Directory $sample_name created"
    # Change to the newly created directory
    cd "${sample_name}"
    echo "Changed directory to $sample_name"

    # Run STAR
    echo "Processing ${sample_name}"
    ulimit -n 20480
    /data/pathology/program/STAR/bin/Linux_x86_64/STAR \
      --runThreadN 32 \
      --genomeDir /data/pathology/cxia/projects/0.ref/10.Sinorhizobium_meliloti/index \
      --readFilesIn /data/pathology/cxia/projects/Giles/Jinpeng/05.Sinorhizobium.RNAseq/01.Data/SRR${i}_1P.fq.gz /data/pathology/cxia/projects/Giles/Jinpeng/05.Sinorhizobium.RNAseq/01.Data/SRR${i}_2P.fq.gz\
      --readFilesCommand zcat \
      --outFilterMultimapNmax 1 \
      --outSAMmultNmax 1 \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix ${sample_name}

    samtools index ${sample_name}Aligned.sortedByCoord.out.bam

    # Get expression FPKM TPM
    /data/pathology/program/stringtie-3.0.0.Linux_x86_64/stringtie -p 32 -G /data/pathology/cxia/projects/0.ref/10.Sinorhizobium_meliloti/GCF_000346065.1_ASM34606v1_genomic.gtf -e -B -A 4.sorted.FPKM.tsv ${sample_name}Aligned.sortedByCoord.out.bam

   cd "${output_dir}"
   echo "${sample_name} finished"
done
```
### 02.3 Expression Cluster 

```bash
# to get TPM 
python3 00.get_protein_TPM.py 01.list /data/pathology/cxia/projects/0.ref/10.Sinorhizobium_meliloti/GCF_000346065.1_ASM34606v1_genomic.tsv ../02.mapping/SRR18299092/4.sorted.FPKM.tsv 02.28d-3.protein_tpm.tsv
paste -d"\t" 02.14d-1.protein_tpm.tsv 02.14d-2.protein_tpm.tsv 02.14d-3.protein_tpm.tsv 02.21d-2.protein_tpm.tsv 02.21d-3.protein_tpm.tsv 02.21d-4.protein_tpm.tsv 02.28d-1.protein_tpm.tsv 02.28d-2.protein_tpm.tsv 02.28d-3.protein_tpm.tsv | awk -F"\t" '{print $1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12"\t"$14"\t"$16"\t"$18}' > 03.TPM4cluster.tsv
```

```r
library(TOmicsVis)
library(readr)

group <- read_tsv("groups.tsv")
group <- as.data.frame(group)
counts <- read_tsv("03.TPM4cluster.tsv")
pdf("042.pca.pdf")
pca_plot(
 sample_gene = counts,
 group_sample = group,
 multi_shape = FALSE,
 xPC = 1,
 yPC = 2,
 point_size = 4,
 text_size = 2,
 fill_alpha = 0.10,
 border_alpha = 0.00,
 legend_pos = "right",
 legend_dir = "vertical",
 ggTheme = "theme_bw"
)
dev.off()

library(ClusterGVis)
library(Biobase)
library(Mfuzz)

counts <- read.csv("03.TPM4cluster.tsv", sep="\t", row.names =1)
pdf("04.getCluster.pdf")
getClusters(counts)
cm <- clusterData(counts, cluster.method = "mfuzz", cluster.num = 4)
ct <- clusterData(counts, cluster.method = "TCseq", cluster.num = 4)
ck <- clusterData(counts, cluster.method = "kmeans", cluster.num = 4)
visCluster(cm, plot.type = "line", ms.col = c("green", "orange", "red"))
visCluster(ct, plot.type = "line", ms.col = c("green", "orange", "red"))
visCluster(ck, plot.type = "line")
dev.off()

#line plot
pdf("04.line_heatmap.pdf")
visCluster(object = cm,
           plot.type = "both",
           ms.col = c("green","orange","red"),
           column_names_rot = 45)
dev.off()

# save cluster information and membership (if exit)
write.csv(cm$wide.res, "04.cm.4clusters.csv",row.names = TRUE, quote = F)
write.csv(ct$wide.res, "04.ct.4clusters.csv",row.names = TRUE, quote = F)
write.csv(ck$wide.res, "04.ck.4clusters.csv",row.names = TRUE, quote = F)
```
![Manhattan plot](https://github.com/chongjing/AlphaFold3_Medicago/blob/main/plot/03.TPM4cluster.Sinorhizobium.jpg)

## 03. AlphaFold3 prediction
### 03.1 MSA
```bash
for i in {001..973}; do

        sed -n "${i}p" 02.list > "${i}.list"; seqkit grep --by-name --pattern-file "${i}.list" --use-regexp --out-file "${i}.fasta" 02.faa && python3 /data/lab/chen/software/alphafold3/fasta2json.py "${i}.fasta"; rm -rf ${i}.list ${i}.fasta;

        echo "Processing ${i}";
        singularity exec --nv --bind /scratch/user/chongjing.xia/20250318_095558/05.Medicago/tmp:/mnt/af_input \
                --bind /scratch/user/chongjing.xia/20250318_095558/05.Medicago/tmp/01.Medicago:/mnt/af_pipeline \
                --bind /data/lab/chen/software/alphafold3/public_databases:/mnt/public_databases \
                --bind /data/lab/chen/software/alphafold3/models:/mnt/models \
                /data/lab/chen/software/alphafold3/AlphaFold3.0.1.sif \
                python /app/alphafold/run_alphafold.py \
                --run_data_pipeline \
                --norun_inference \
                --json_path=/mnt/af_input/${i}.json \
                --model_dir=/mnt/models \
                --db_dir=/mnt/public_databases \
                --output_dir=/mnt/af_pipeline;
        echo "Done ${i}";
        rm -rf ./${i}.json
done
```

### 03.2 Inference
```bash
#!/bin/bash
#SBATCH --partition=cahnrs,cahnrs_bigmem,camas,kamiak
#SBATCH --job-name=301
#SBATCH --output=301.o
#SBATCH --error=301.e
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:1
#SBATCH --mem=64000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chongjing.xia@wsu.edu

cd /scratch/user/chongjing.xia/20250510_234922/05.Medicago/af_inference/input

module load singularity/3.8.0

# Inference
for i in {001..653}; do
        echo "Processing 301_${i}"
        cd /scratch/user/chongjing.xia/20250510_234922/05.Medicago/af_inference/input
        ## prepare JSON input from AlphaFold3 processed pipeline for inference
        sed -n "301p;$((10#${i} + 333))p" /scratch/user/chongjing.xia/20250510_234922/05.Medicago/01.input/02.list > 301_${i}.list
        seqkit grep --by-name --pattern-file 301_${i}.list --use-regexp --out-file 301_${i}.fasta /scratch/user/chongjing.xia/20250510_234922/05.Medicago/01.input/02.faa
        python3 /data/lab/chen/software/alphafold3/fasta2json.py 301_${i}.fasta
        python3 /data/lab/chen/software/alphafold3/AF3_MSA_formating.py 301_${i}.json /scratch/user/chongjing.xia/20250510_234922/05.Medicago/tmp/01.Medicago/301/301_data.json /scratch/user/chongjing.xia/20250510_234922/05.Medicago/tmp/01.Medicago/$((10#${i} + 333))/$((10#${i} + 333))_data.json
        rm -rf 301_${i}.list 301_${i}.fasta

        singularity exec --nv --bind /scratch/user/chongjing.xia/20250510_234922/05.Medicago/af_inference/input:/mnt/input \
                --bind /scratch/user/chongjing.xia/20250510_234922/05.Medicago/02.output:/mnt/output \
                --bind /data/lab/chen/software/alphafold3/public_databases:/mnt/public_databases \
                --bind /data/lab/chen/software/alphafold3/models:/mnt/models \
                --env 'XLA_PYTHON_CLIENT_PREALLOCATE=false' \
                --env 'TF_FORCE_UNIFIED_MEMORY=true' \
                --env 'XLA_CLIENT_MEM_FRACTION=3.2' \
                --env 'XLA_FLAGS=--xla_gpu_enable_triton_gemm=false' \
                /data/lab/chen/software/alphafold3/AlphaFold3.0.1.sif \
                python /app/alphafold/run_alphafold.py \
                --norun_data_pipeline \
                --run_inference \
                --json_path=/mnt/input/301_${i}.json \
                --model_dir=/mnt/models \
                --db_dir=/mnt/public_databases \
                --output_dir=/mnt/output
        echo "Done 301_${i}";
done
```

### 03.3 AF3 summary & visualization
```bash
cd /scratch/user/chongjing.xia/20250427_013151/05.Medicago
#extract iPTM, PTM, for all interactions
python3 03.extract_iptm_ptm_from_summary.json.py 001.20250618.csv

#filter using iPTM >= 0.6 && PTM >= 0.6
awk -F"," '$2 >= 0.6 && $3 >= 0.6 {print $1}' 001.20251006.csv | sort -k1,1n > 002.interaction.threshold_0.6.list
# make MasterTable, with MSA info (directory structure matters)
python3 03.extract_iptm_ptm_from_summary.json.py 002.interaction.threshold_0.6.list
# add protein names and sequences
python3 03.4.add_ProteinSeq.py

/data/pathology/program/Miniforge3/envs/R4.2.3/bin/R
```

```r
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
data <- read.delim("002.20250508.tsv", header = TRUE)
str(data)

# Adapted from a GWAS manhatten plot
# Calculate chromosome sizes (max position for each chromosome)
chromosome_sizes <- data %>%
  group_by(ID) %>%
  summarize(Size = max(Position, na.rm = TRUE))
chromosome_sizes


# Calculate cumulative chromosome sizes for proper X-axis scaling
chromosome_sizes <- chromosome_sizes %>%
  mutate(Cumulative_Size = cumsum(Size) - Size)

# Merge cumulative sizes back into the main data
data <- data %>%
  left_join(chromosome_sizes, by = "ID") %>%
  mutate(Cumulative_Position = Position + Cumulative_Size)

# Create a color palette for the chromosomes
chromosome_colors <- rep(rainbow(7), length.out = length(unique(data$ID)))

# Calculate midpoints for chromosome labels
chromosome_midpoints <- data %>%
  group_by(ID) %>%
  summarize(Midpoint = mean(Cumulative_Position))

pdf("002.iPTM.pdf", 20,9)
ggplot(data, aes(x = Cumulative_Position, y = 0.5*iptm+0.5*ptm, color = factor(ID))) +
  geom_point(size = 1.0, alpha = 1.0) +  # Adjust point size and transparency
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", linewidth = 0.5) +  # Threshold line at 0.8
  scale_color_manual(values = chromosome_colors) +
  labs(y = "AlphaFold3 0.5*iPTM+0.5*PTM", x = "gene ID") +
  theme_minimal() +
  theme(legend.position = "none",  # Remove legend
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.ticks = element_line(color = "black"),  # Add axis ticks 
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate chromosome names
    ) +
  scale_y_continuous(limits = c(0, 1)) +  # Set Y-axis limits
  scale_x_continuous(breaks = chromosome_midpoints$Midpoint,  # Add chromosome labels at midpoints
    labels = chromosome_midpoints$ID,  # Use chromosome names as labels
    expand = c(0.02, 0.02)  # Reduce padding around X-axis
    )
dev.off()
```
A manhattan plot of iPTM and PTM.
![Manhattan plot](https://github.com/chongjing/AlphaFold3_Medicago/blob/main/plot/002_0.5iPTM_0.5PTM.2.jpeg)

### 03.4 Distribution of number of targets

```R
library(ggplot2)
library(dplyr)

setwd("/home/cx264/rds/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/04.AlphaFold3/03.analysis/03.Frequency")

data <- read.table("002.20251006.iPTM_PTM.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Filter data
score_threshold <- 0.8
data_filtered <- data %>%
  filter(Score >= score_threshold) %>%
  # Ensure IDs are treated correctly
  mutate(Medicago = as.character(Medicago),
         Rhizobium = as.character(Rhizobium))

# Count how many UNIQUE Rhizobium proteins each Medicago protein targets
medicago_degree <- data_filtered %>%
  group_by(Medicago) %>%
  summarise(
    Num_Rhizobium_Targets = n_distinct(Rhizobium),
    .groups = 'drop'
  )

#Distribution of Medicago Protein Targets: How many Medicago proteins target X number of Rhizobium proteins?
p5 <- ggplot(medicago_degree, aes(x = Num_Rhizobium_Targets)) +
  geom_histogram(binwidth = 1, fill = "purple", color = "white", alpha = 0.8) + geom_text(stat = "bin", binwidth = 1,
            aes(label = after_stat(count)),
            vjust = -0.5, size = 10, color = "black") +
#  geom_vline(xintercept = mean(medicago_degree$Num_Rhizobium_Targets),
#             color = "red", linetype = "dashed", size = 1) +
  scale_x_continuous(breaks = 1:max(medicago_degree$Num_Rhizobium_Targets)) +
  labs(
#    title = "Distribution of Medicago Protein Targets",
#    subtitle = paste0("Number of Rhizobium proteins targeted per Medicago protein\n(Score cutoff >= ", score_threshold, ")"),
    x = "Number of Rhizobium Targets (per Medicago protein)",
    y = "Number of Medicago Proteins"
  ) +
  theme_classic(base_size = 25) +   # base_size controls overall font size
  theme(
#    plot.title = element_text(face = "bold", size = 16),
#    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 30),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black", size = 3),  # ensure axis lines are visible
    panel.grid = element_blank(),                           # no grid (already absent in classic)
    plot.margin = margin(5, 5, 5, 5)                    # add some margin
  )

ggsave("005.Medicago_Distribution.0.8.pdf", plot = p5, width = 15, height = 10)
svg("005.Medicago_Distribution.0.8.2.svg", 15, 10)
p5
dev.off()

#Top 10 Medicago Hub Proteins that interact with the most Rhizobium proteins
# Let's look at the top 10 Medicago proteins that interact with the most Rhizobium proteins
top_hubs <- medicago_degree %>%
  arrange(desc(Num_Rhizobium_Targets)) %>%
  head(10)
print(top_hubs)

p6 <- ggplot(top_hubs, aes(x = reorder(Medicago, -Num_Rhizobium_Targets), y = Num_Rhizobium_Targets)) +
  geom_bar(stat = "identity", fill = "coral", alpha = 0.8) +
  geom_text(aes(label = Num_Rhizobium_Targets), vjust = -0.5, size = 5) +
  labs(
#    title = "Top 10 Medicago Hub Proteins",
    x = "Medicago Protein",
    y = "Count of Rhizobium Targets"
  ) +
    theme_classic(base_size = 25) +   # base_size controls overall font size
  theme(
#    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 25),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black", size = 2),  # ensure axis lines are visible
    panel.grid = element_blank(),                           # no grid (already absent in classic)
    plot.margin = margin(5, 5, 5, 5)                    # add some margin
  )

ggsave("005.Top_Hubs.0.8.pdf", plot = p6, width = 8, height = 6)
svg("005.Top_Hubs.0.8.2.svg")
p6
dev.off()
```

<table>
  <tr>
    <td><img src="https://github.com/chongjing/AlphaFold3_Medicago/blob/main/plot/005.Medicago_Distribution.0.8.jpeg" alt="Image 1" width="400"/></td>
    <td><img src="https://github.com/chongjing/AlphaFold3_Medicago/blob/main/plot/005.Top_Hubs.0.8.jpeg" alt="Image 2" width="400"/></td>
  </tr>
</table>



