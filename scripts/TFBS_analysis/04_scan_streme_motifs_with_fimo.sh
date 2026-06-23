
# ============================================================
# 04) Scan selected STREME motifs with FIMO
#
# Goal:
#   Extract STREME motifs 1-4 and scan them in positive and
#   background promoter-window FASTA files using FIMO.
#
# Run after script 04. Script 05 depends on
# the outputs generated here.
#
# Run:
#   bash 04_scan_streme_motifs_with_fimo.sh
# ============================================================



# ============================================================
# 1) Paths
# ============================================================

PROJECT_DIR="${PROJECT_DIR:-/path/to/project}"
RESULTS_DIR="${RESULTS_DIR:-$PROJECT_DIR/results}"

MEME_SIF="${MEME_SIF:-/path/to/meme.sif}"

STREME_DISCOVERY_DIR="${STREME_DISCOVERY_OUTPUT_DIR:-$RESULTS_DIR/streme_motif_discovery}"
INPUT_DIR="${INPUT_DIR:-$STREME_DISCOVERY_DIR/fasta}"
POS_FASTA="${POSITIVE_FASTA:-$INPUT_DIR/positive_shift_promoter_windows_50_150bp_up.fa}"
BG_FASTA="${BACKGROUND_FASTA:-$INPUT_DIR/background_all_cage_promoter_windows_50_150bp_up.fa}"

STREME_OUT_FINAL="${STREME_OUT_DIR:-$STREME_DISCOVERY_DIR/streme}"
STREME_FINAL="$STREME_OUT_FINAL/streme_FINAL.txt"

if [[ ! -f "$STREME_FINAL" && -f "$STREME_OUT_FINAL/streme.txt" ]]; then
  echo "WARNING: $STREME_FINAL not found."
  echo "Using fallback: $STREME_OUT_FINAL/streme.txt"
  STREME_FINAL="$STREME_OUT_FINAL/streme.txt"
fi

OUT_ROOT="${STREME_MOTIF_OUTPUT_DIR:-$RESULTS_DIR/streme_motif_analysis}"
MOTIF_DIR="$OUT_ROOT/motifs"
FIMO_DIR="$OUT_ROOT/fimo"
LOG_DIR="$OUT_ROOT/logs"
TABLE_DIR="$OUT_ROOT/tables"

mkdir -p "$MOTIF_DIR" "$FIMO_DIR" "$LOG_DIR" "$TABLE_DIR"


# ============================================================
# 2) FIMO parameters
# ============================================================

# Same conservative threshold as before. If you want a sensitivity run, try 1e-3.
FIMO_THRESH="1e-4"
FIMO_MAX_STORED_SCORES=10000000


# ============================================================
# 3) Checks
# ============================================================

for f in "$MEME_SIF" "$POS_FASTA" "$BG_FASTA" "$STREME_FINAL"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing file: $f"
    exit 1
  fi
done

echo "MEME_SIF: $MEME_SIF"
echo "STREME_OUT_FINAL: $STREME_OUT_FINAL/"
echo "STREME_FINAL: $STREME_FINAL"
echo "POS_FASTA: $POS_FASTA"
echo "BG_FASTA: $BG_FASTA"


# ============================================================
# 4) Extract only STREME motifs 1, 2, 3
# ============================================================

STREME_MOTIFS_1_2_3_4="$MOTIF_DIR/STREME_motifs_1_2_3_4.meme"

awk '
BEGIN {in_motif=0; keep=1}
(/^MOTIF /) {
  in_motif=1
  if ($2 ~ /^1-/ || $2 ~ /^2-/ || $2 ~ /^3-/ || $2 ~ /^4-/) {
    keep=1
    print
  } else {
    keep=0
  }
  next
}
{
  if (!in_motif || keep) print
}
' "$STREME_FINAL" > "$STREME_MOTIFS_1_2_3_4"

echo
echo "Extracted motif file:"
echo "$STREME_MOTIFS_1_2_3_4"
echo "Motifs:"
grep '^MOTIF ' "$STREME_MOTIFS_1_2_3_4" || true

N_MOTIFS=$(grep -c '^MOTIF ' "$STREME_MOTIFS_1_2_3_4" || true)
if [[ "$N_MOTIFS" -ne 4 ]]; then
  echo "ERROR: expected exactly 4 motifs, found $N_MOTIFS"
  exit 1
fi


# ============================================================
# 5) Run FIMO on positive and background FASTA
# ============================================================

echo
echo "Running FIMO on positive FASTA"
singularity exec "$MEME_SIF" fimo \
  --oc "$FIMO_DIR/fimo_streme_motifs_positive" \
  --thresh "$FIMO_THRESH" \
  --max-stored-scores "$FIMO_MAX_STORED_SCORES" \
  "$STREME_MOTIFS_1_2_3_4" \
  "$POS_FASTA" \
  > "$LOG_DIR/fimo_streme_motifs_positive.log" 2>&1

echo
echo "Running FIMO on background FASTA"
singularity exec "$MEME_SIF" fimo \
  --oc "$FIMO_DIR/fimo_streme_motifs_background" \
  --thresh "$FIMO_THRESH" \
  --max-stored-scores "$FIMO_MAX_STORED_SCORES" \
  "$STREME_MOTIFS_1_2_3_4" \
  "$BG_FASTA" \
  > "$LOG_DIR/fimo_streme_motifs_background.log" 2>&1


# ============================================================
# 6) Manifest
# ============================================================

cat > "$TABLE_DIR/00_streme_motif_fimo_manifest.tsv" <<EOF
parameter	value
MEME_SIF	$MEME_SIF
STREME_OUT_FINAL	$STREME_OUT_FINAL/
STREME_FINAL	$STREME_FINAL
STREME_MOTIFS_1_2_3_4	$STREME_MOTIFS_1_2_3_4
POSITIVE_FASTA	$POS_FASTA
BACKGROUND_FASTA	$BG_FASTA
FIMO_THRESH	$FIMO_THRESH
FIMO_MAX_STORED_SCORES	$FIMO_MAX_STORED_SCORES
EOF

echo
echo "DONE 04_scan_streme_motifs_with_fimo"
echo "Positive FIMO:   $FIMO_DIR/fimo_streme_motifs_positive/fimo.tsv"
echo "Background FIMO: $FIMO_DIR/fimo_streme_motifs_background/fimo.tsv"
echo
echo "Next:"
echo "Rscript 05_streme_motif_enrichment_and_variant_specific.R"
