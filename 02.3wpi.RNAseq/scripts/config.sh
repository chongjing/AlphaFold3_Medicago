#!/bin/bash
# ============================================================================
# config.sh — Environment configuration for the Medicago RNA-seq pipeline
# ============================================================================
# Source this file at the top of each .srun script, or set these variables
# in your shell before running:
#   source config.sh
#
# Adjust paths below to match your system. Tool names (STAR, htseq-count,
# stringtie, Rscript, python3) must be on PATH.
# ============================================================================

# Conda environment with R 4.3.2 (used for matrix build)
export CONDA_R432_BIN="${CONDA_R432_BIN:-/home/cx264/program/anaconda3/envs/R4.3.2/bin}"

# Conda environment with R 4.4.3 (used for DE + enrichment)
export CONDA_R443_BIN="${CONDA_R443_BIN:-/home/cx264/program/anaconda3/envs/R4.4.3/bin}"

# Conda base (for python3)
export CONDA_BIN="${CONDA_BIN:-/home/cx264/program/anaconda3/bin}"

# Java (for qualimap)
export JAVA_HOME="${JAVA_HOME:-/home/cx264/program/jdk1.8.0_491}"

# Reference genome directory (must contain fasta/ and genes/ subdirs)
export REF_GENOME_DIR="${REF_GENOME_DIR:-/home/cx264/project/00.ref/M.truncatula_genome_v5}"

# Functional annotation directory (for GO/KEGG enrichment)
export ANNOTATION_DIR="${ANNOTATION_DIR:-/home/cx264/project/00.ref/M.truncatula_genome_v5/20220708_MtrunA17r5.0-ANR-EGN-r1.9_FunctionalAnnotation}"

# Add tool directories to PATH
export PATH="${CONDA_R432_BIN}:${CONDA_BIN}:${JAVA_HOME}/bin:$PATH"
