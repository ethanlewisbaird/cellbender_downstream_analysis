#!/bin/bash
#SBATCH --job-name=CB_dwnstrm_analysis_EB
#SBATCH --output=dwnstrm_CB_EB.out
#SBATCH --error=dwnstrm_CB_EB.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=72:00:00
#SBATCH --mem=240G


INPUT_DIR="/data/ebaird/scRNAseq/CB"
OUTPUT_DIR="/data/ebaird/scRNAseq/dwnstrm_CB"


INPUT_FILE="${INPUT_DIR}/2196_CB_output_filtered.h5"
SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/2196_CB_output_filtered.h5"
Rscript --vanilla /data/ebaird/scRNAseq/dwnstrm_CB/dwnstrm_CB.r $INPUT_FILE $SAMPLE_OUTPUT_DIR $2196

INPUT_FILE="${INPUT_DIR}/2197_CB_output_filtered.h5"
SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/2197_CB_output_filtered.h5"
Rscript --vanilla /data/ebaird/scRNAseq/dwnstrm_CB/dwnstrm_CB.r $INPUT_FILE $SAMPLE_OUTPUT_DIR $2197

INPUT_FILE="${INPUT_DIR}/2198_CB_output_filtered.h5"
SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/2198_CB_output_filtered.h5"
Rscript --vanilla /data/ebaird/scRNAseq/dwnstrm_CB/dwnstrm_CB.r $INPUT_FILE $SAMPLE_OUTPUT_DIR $2198

INPUT_FILE="${INPUT_DIR}/2199_CB_output_filtered.h5"
SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/2199_CB_output_filtered.h5"
Rscript --vanilla /data/ebaird/scRNAseq/dwnstrm_CB/dwnstrm_CB.r $INPUT_FILE $SAMPLE_OUTPUT_DIR $2199