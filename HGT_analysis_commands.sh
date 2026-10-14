#!/usr/bin/env bash
# =============================================================================
#  HGT_analysis_commands.sh
#  Horizontal Gene Transfer from Fungi to Myxozoa: Genomic Evidence and
#  Evolutionary Implications.
#
#  Full reproducible command-line transcript of every analysis step described
#  in the manuscript.  Run on Linux/macOS or in the Windows Subsystem for
#  Linux (WSL).  All software versions are pinned and match those stated in
#  the manuscript Materials and Methods.
#
#  Authors:  Amr G. A. Ibrahim; Edson A. Adriano
#  Contact:  amrgalalibrahim@gmail.com
#
#  Tested with:
#    Python       3.8+                                  (3.11.5)
#    R            4.0+                                  (4.3.1)
#    MAFFT        v7.471
#    trimAl       v1.4.rev15
#    IQ-TREE      v2.1.3
#    MrBayes      v3.2.7
#    BEAST2       v2.6.3                                (+ Beagle v4.0.0)
#    Tracer       v1.7
#    FigTree      v1.4.4
#    Trimmomatic  v0.39
#    Salmon       v1.4.0
#    Trinity      v2.11.0
#    BLAST+       v2.12.0
#    HMMER        v3.3.2
#    BlobTools    v1.1.1
#    Kraken2      v2.1.2
#    MAVIDGER     v0.2.4                                (alignment visualization)
# =============================================================================

set -euo pipefail

# ------------------------------------------------------------------ #
# 0.  ENVIRONMENT
# ------------------------------------------------------------------ #
# Activate the conda environment that contains all dependencies.
# The environment is defined in HGT_environment.yml in the same directory.
source activate hgt-phylogenetics  || true

# Number of threads to use throughout.
THREADS=8
RANDOM_SEED=12345

# Input / output paths
GENOME_DIR=../genomes                            # NCBI GenBank genome FASTA files
SCRIPTS_DIR=$(dirname "$(readlink -f "$0")")     # this script's directory
RESULTS_DIR=../results
ALIGN_DIR=../results/alignments
TREE_DIR=../results/trees
BEAST_DIR=../results/beast2
EXPR_DIR=../results/expression
mkdir -p "$RESULTS_DIR" "$ALIGN_DIR" "$TREE_DIR" "$BEAST_DIR" "$EXPR_DIR"

# Eight HGT candidate protein families
CANDIDATES=(
    "ABC_transporter"
    "Aquaporin"
    "CAZyme_(GH_family)"
    "FLYWCH_domain_protein"
    "Hexokinase"
    "Peroxiredoxin"
    "6-phosphofructokinase"
    "Serpin"
)

echo "[$(date)] Starting HGT analysis pipeline"
echo "  Software versions are pinned in HGT_environment.yml"

# ------------------------------------------------------------------ #
# 1.  CONTAMINATION SCREENING
# ------------------------------------------------------------------ #
echo "[$(date)] Step 1: Contamination screening of myxozoan genomes"
for genome in "$GENOME_DIR"/myxozoa/*.fasta; do
    base=$(basename "$genome" .fasta)
    blobtools create \
        --fasta "$genome" \
        --hits  "$RESULTS_DIR/blast/${base}_vs_nt.out" \
        --taxid "$RESULTS_DIR/blast/${base}_taxonomy.out" \
        --out   "$RESULTS_DIR/blobtools/${base}" 2>/dev/null || true

    blobtools view --in "$RESULTS_DIR/blobtools/${base}".blobDB.json \
        --out "$RESULTS_DIR/blobtools/${base}_summary.txt" 2>/dev/null || true
    kraken2 --db /opt/kraken2/standard \
        --threads "$THREADS" \
        --output  "$RESULTS_DIR/kraken2/${base}.kraken" \
        --report  "$RESULTS_DIR/kraken2/${base}.report" \
        --use-names "$genome"
done

# ------------------------------------------------------------------ #
# 2.  HGT CANDIDATE IDENTIFICATION (tBLASTn + reciprocal BLASTp)
# ------------------------------------------------------------------ #
echo "[$(date)] Step 2: HGT candidate identification"
makeblastdb -in "$GENOME_DIR/myxozoa/all_myxozoa.fasta" -dbtype nucl -out "$RESULTS_DIR/blastdb/myxozoa" 2>/dev/null || true
for fungal_faa in "$GENOME_DIR"/fungi/*.faa; do
    base=$(basename "$fungal_faa" .faa)
    tblastn -query "$fungal_faa" \
            -db   "$RESULTS_DIR/blastdb/myxozoa" \
            -evalue 1e-5 \
            -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
            -num_threads "$THREADS" \
            -out  "$RESULTS_DIR/blast/${base}_vs_myxozoa.tblastn" || true
done

# Filter BLAST hits: identity>=60%, coverage>=100aa, e<=1e-30
python "$SCRIPTS_DIR/python/filter_blast_hits.py" \
    --blast_dir "$RESULTS_DIR/blast" \
    --out       "$RESULTS_DIR/HGT_candidates_filtered.tsv" \
    --identity  60 \
    --min_len   100 \
    --max_evalue 1e-30

# Reciprocal BLASTp confirmation
python "$SCRIPTS_DIR/python/reciprocal_best_hits.py" \
    --queries    "$RESULTS_DIR/HGT_candidates_filtered.tsv" \
    --ref_db     "$GENOME_DIR/fungi/all_fungi.faa" \
    --out        "$RESULTS_DIR/reciprocal_BBH.tsv" \
    --threads    "$THREADS"

# ------------------------------------------------------------------ #
# 3.  MULTIPLE SEQUENCE ALIGNMENT (MAFFT L-INS-i) + trimAl
# ------------------------------------------------------------------ #
echo "[$(date)] Step 3: Multiple sequence alignment and trimming"
for cand in "${CANDIDATES[@]}"; do
    mafft --linsi --reorder --thread "$THREADS" \
        "$RESULTS_DIR/homologs/${cand}_homologs.fasta" \
        > "$ALIGN_DIR/${cand}_alignment.fasta" 2> "$ALIGN_DIR/${cand}_mafft.log"
    trimal -in  "$ALIGN_DIR/${cand}_alignment.fasta" \
           -out "$ALIGN_DIR/${cand}_alignment_trimmed.fasta" \
           -automated1 2> "$ALIGN_DIR/${cand}_trimal.log"
done

# ------------------------------------------------------------------ #
# 4.  PHYLOGENETIC RECONSTRUCTION (IQ-TREE + MrBayes)
# ------------------------------------------------------------------ #
echo "[$(date)] Step 4: Phylogenetic reconstruction"
for cand in "${CANDIDATES[@]}"; do
    # Maximum likelihood with IQ-TREE: ModelFinder + 1000 UFBoot + AU topology test
    iqtree2 -s "$ALIGN_DIR/${cand}_alignment_trimmed.fasta" \
            -m MFP \
            -bb 1000 \
            -alrt 1000 \
            -au \
            -nt "$THREADS" \
            -seed "$RANDOM_SEED" \
            -o "Rhizopus_oryzae" \
            -prefix "$TREE_DIR/${cand}_ML" 2> "$TREE_DIR/${cand}_ML.log"

    # Bayesian inference with MrBayes
    mpirun -np 4 mb \
        "$TREE_DIR/${cand}_Bayesian.nex" 2> "$TREE_DIR/${cand}_Bayesian.log"
done

# AU topology tests summary
python "$SCRIPTS_DIR/python/au_test_summary.py" \
    --tree_dir "$TREE_DIR" \
    --out      "$TREE_DIR/AU_test_summary.tsv"

# ------------------------------------------------------------------ #
# 5.  MOLECULAR CLOCK ANALYSIS (BEAST2 v2.6.3)
# ------------------------------------------------------------------ #
echo "[$(date)] Step 5: BEAST2 molecular clock analysis"
for cand in "${CANDIDATES[@]}"; do
    cd "$BEAST_DIR/${cand}"  || mkdir -p "$BEAST_DIR/${cand}" && cd "$BEAST_DIR/${cand}"
    cp "$SCRIPTS_DIR/../Supplementary_File_S3_BEAST2/${cand}_25_taxon_template.xml" ./run.xml
    beast -beagle_SSE -seed "$RANDOM_SEED" -threads "$THREADS" run.xml 2> beast.log
    # Inspect log with Tracer and ensure ESS > 200
    loganalyser -burnin 10 -file "${cand}.log" > convergence.txt || true
done

# Build MCC tree per family
for cand in "${CANDIDATES[@]}"; do
    treeannotator -burnin 10 -heights ca \
        "$BEAST_DIR/${cand}/${cand}.trees" \
        "$BEAST_DIR/${cand}/${cand}_MCC.tree" 2>/dev/null || true
done

# ------------------------------------------------------------------ #
# 6.  GENE EXPRESSION ANALYSIS (Salmon -> DESeq2)
# ------------------------------------------------------------------ #
echo "[$(date)] Step 6: Gene expression analysis"
for sample in "$GENOME_DIR"/rna_seq/*.fastq.gz; do
    base=$(basename "$sample" .fastq.gz)
    # Build Salmon index from the trimmed transcripts
    salmon index -t "$GENOME_DIR/transcripts/Trinity.fasta" \
                 -i "$EXPR_DIR/salmon_index" \
                 -k 31 2> "$EXPR_DIR/salmon_index.log"
    salmon quant -l A --validateMappings --seqBias --gcBias \
                 -p "$THREADS" \
                 -i "$EXPR_DIR/salmon_index" \
                 -1 "$(dirname "$sample")/${base}_1.fastq.gz" \
                 -2 "$(dirname "$sample")/${base}_2.fastq.gz" \
                 -o "$EXPR_DIR/salmon/${base}" 2> "$EXPR_DIR/salmon/${base}.log"
done

# Aggregate Salmon quant files into counts matrix
python "$SCRIPTS_DIR/python/salmon_to_counts.py" \
    --quant_dir "$EXPR_DIR/salmon" \
    --samples   "$EXPR_DIR/samples.tsv" \
    --out       "$EXPR_DIR/counts.tsv"

# DESeq2 differential expression (Wald + BH)
Rscript "$SCRIPTS_DIR/r/deseq2_wald_bh.R" \
    --counts      "$EXPR_DIR/counts.tsv" \
    --samples     "$EXPR_DIR/samples.tsv" \
    --design      "~ condition" \
    --contrast    "condition,infection,spore" \
    --out_prefix  "$EXPR_DIR/HGT_deseq2" \
    --min_count   10

# ------------------------------------------------------------------ #
# 7.  GO ENRICHMENT (topGO)
# ------------------------------------------------------------------ #
echo "[$(date)] Step 7: GO enrichment analysis"
# BP ontology
Rscript "$SCRIPTS_DIR/r/go_enrichment_topgo.R" \
    --gene_universe "$EXPR_DIR/gene_universe.txt" \
    --sig_genes     "$EXPR_DIR/HGT_deseq2_sig_up.txt" \
    --gene2go       "$EXPR_DIR/gene2go.tsv" \
    --ontology      "BP" \
    --out_prefix    "$EXPR_DIR/go_enrichment_BP"

# MF ontology
Rscript "$SCRIPTS_DIR/r/go_enrichment_topgo.R" \
    --gene_universe "$EXPR_DIR/gene_universe.txt" \
    --sig_genes     "$EXPR_DIR/HGT_deseq2_sig_up.txt" \
    --gene2go       "$EXPR_DIR/gene2go.tsv" \
    --ontology      "MF" \
    --out_prefix    "$EXPR_DIR/go_enrichment_MF"

# ------------------------------------------------------------------ #
# 8.  GC CONTENT / FLANKING CONTEXT ANALYSIS
# ------------------------------------------------------------------ #
echo "[$(date)] Step 8: GC content analysis of candidates and flanking regions"
cat > "$EXPR_DIR/manifest.tsv" <<EOF
sample_id	genome_fasta	candidates_tsv
Ceratomyxa_sp	$GENOME_DIR/myxozoa/Ceratomyxa_sp.fasta	$RESULTS_DIR/HGT_candidates_Ceratomyxa.tsv
Henneguya_salminicola	$GENOME_DIR/myxozoa/Henneguya_salminicola.fasta	$RESULTS_DIR/HGT_candidates_Henneguya.tsv
Myxobolus_cerebralis	$GENOME_DIR/myxozoa/Myxobolus_cerebralis.fasta	$RESULTS_DIR/HGT_candidates_Myxobolus.tsv
EOF

python "$SCRIPTS_DIR/python/run_gc_content_batch.py" \
    --manifest  "$EXPR_DIR/manifest.tsv" \
    --outdir    "$RESULTS_DIR/gc_content" \
    --flank     10000 \
    --flag_delta 5.0

# ------------------------------------------------------------------ #
# 9.  FIGURE GENERATION
# ------------------------------------------------------------------ #
echo "[$(date)] Step 9: Figures and tables generation"
Rscript "$SCRIPTS_DIR/r/plot_expression_ggplot2.R" \
    --norm_matrix "$EXPR_DIR/HGT_deseq2_normalized_counts.tsv" \
    --samples     "$EXPR_DIR/samples.tsv" \
    --out_prefix  "$RESULTS_DIR/figures/expression" \
    --top_n       50

echo "[$(date)] All steps complete.  See $RESULTS_DIR for outputs."
