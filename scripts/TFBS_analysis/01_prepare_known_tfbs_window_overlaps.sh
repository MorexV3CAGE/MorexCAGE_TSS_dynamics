
# ============================================================
# 01) Prepare known TFBS overlaps with promoter windows
#
# Goal:
#   Convert positive/background FASTA headers to BED windows and intersect
#   them with a genome-wide known TFBS BED file using a Python coordinate
#   overlap sweep.
#
# Run:
#   bash 01_prepare_known_tfbs_window_overlaps.sh
# ============================================================



# ============================================================
# 1) Paths
# ============================================================

PROJECT_DIR="${PROJECT_DIR:-/path/to/project}"
RESULTS_DIR="${RESULTS_DIR:-$PROJECT_DIR/results}"

INPUT_DIR="${INPUT_DIR:-$PROJECT_DIR/input/promoter_windows}"
POS_FASTA="${POSITIVE_FASTA:-$INPUT_DIR/positive_shift_promoter_windows_50_150bp_up.fa}"
BG_FASTA="${BACKGROUND_FASTA:-$INPUT_DIR/background_all_cage_promoter_windows_50_150bp_up.fa}"

TFBS_BED="${KNOWN_TFBS_BED:-$PROJECT_DIR/input/known_tfbs/known_tfbs_genomewide.bed.gz}"
FAMILY_MAP="${TF_FAMILY_MAP:-$PROJECT_DIR/input/known_tfbs/tf_family_map.tsv}"

OUT_ROOT="${KNOWN_TFBS_OUTPUT_DIR:-$RESULTS_DIR/known_tfbs_module_analysis}"
BED_DIR="$OUT_ROOT/windows"
INTERSECT_DIR="$OUT_ROOT/overlaps"
LOG_DIR="$OUT_ROOT/logs"
TABLE_DIR="$OUT_ROOT/tables"

mkdir -p "$BED_DIR" "$INTERSECT_DIR" "$LOG_DIR" "$TABLE_DIR"

for f in "$POS_FASTA" "$BG_FASTA" "$TFBS_BED"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing required file: $f"
    exit 1
  fi
done

echo "PROJECT_DIR: $PROJECT_DIR"
echo "RESULTS_DIR: $RESULTS_DIR"
echo "TFBS_BED:    $TFBS_BED"
echo "FAMILY_MAP:  $FAMILY_MAP"
echo "POS_FASTA:   $POS_FASTA"
echo "BG_FASTA:    $BG_FASTA"
echo "OUT_ROOT:    $OUT_ROOT"


# ============================================================
# 2) Detect TFBS BED column count
# ============================================================

if [[ "$TFBS_BED" == *.gz ]]; then
  TFBS_NCOL=$(zcat "$TFBS_BED" | awk 'BEGIN{FS=OFS="\t"} !/^#/ && NF>0 {print NF; exit}')
else
  TFBS_NCOL=$(awk 'BEGIN{FS=OFS="\t"} !/^#/ && NF>0 {print NF; exit}' "$TFBS_BED")
fi

if [[ -z "$TFBS_NCOL" ]]; then
  echo "ERROR: could not detect TFBS BED columns."
  exit 1
fi

echo "Detected TFBS BED columns: $TFBS_NCOL"
echo "First 3 TFBS BED rows:"
if [[ "$TFBS_BED" == *.gz ]]; then
  zcat "$TFBS_BED" | awk 'BEGIN{FS=OFS="\t"} !/^#/ && NF>0 {print; n++; if(n==3) exit}'
else
  awk 'BEGIN{FS=OFS="\t"} !/^#/ && NF>0 {print; n++; if(n==3) exit}' "$TFBS_BED"
fi


# ============================================================
# 3) Convert FASTA headers to BED windows
# ============================================================

POS_WINDOWS_BED="$BED_DIR/positive_promoter_windows.bed"
BG_WINDOWS_BED="$BED_DIR/background_promoter_windows.bed"

python3 - "$POS_FASTA" "positive" "$POS_WINDOWS_BED" > "$LOG_DIR/positive_windows_bed.log" 2>&1 <<'PY'
import sys, re
fasta, set_name, out_bed = sys.argv[1], sys.argv[2], sys.argv[3]
coord_re = re.compile(r"::([^:]+):(\d+)-(\d+)(?:\(([+-])\))?")
fallback_re = re.compile(r"([^:]+):(\d+)-(\d+)(?:\(([+-])\))?")
tss_re = re.compile(r"TSS(\d+)")
n = bad = 0
with open(out_bed, "w") as out, open(fasta) as f:
    for line in f:
        if not line.startswith(">"):
            continue
        h = line[1:].strip()
        fasta_name = h.split("::")[0]
        m = coord_re.search(h) or fallback_re.search(h)
        if not m:
            bad += 1
            continue
        chrom, start0, end0 = m.group(1), int(m.group(2)), int(m.group(3))
        strand = m.group(4) if m.group(4) else "."
        parts = fasta_name.split("|")
        gene_id = parts[1] if len(parts) > 1 else "NA"
        side = parts[2] if len(parts) > 2 else "NA"
        mtss = tss_re.search(fasta_name)
        tss = mtss.group(1) if mtss else "NA"
        out.write("\t".join(map(str, [chrom, start0, end0, fasta_name, set_name, side, gene_id, tss, strand])) + "\n")
        n += 1
print(f"Wrote {n} BED windows to {out_bed}", file=sys.stderr)
if bad:
    print(f"WARNING skipped {bad} headers without coordinates", file=sys.stderr)
PY

python3 - "$BG_FASTA" "background" "$BG_WINDOWS_BED" > "$LOG_DIR/background_windows_bed.log" 2>&1 <<'PY'
import sys, re
fasta, set_name, out_bed = sys.argv[1], sys.argv[2], sys.argv[3]
coord_re = re.compile(r"::([^:]+):(\d+)-(\d+)(?:\(([+-])\))?")
fallback_re = re.compile(r"([^:]+):(\d+)-(\d+)(?:\(([+-])\))?")
tss_re = re.compile(r"TSS(\d+)")
n = bad = 0
with open(out_bed, "w") as out, open(fasta) as f:
    for line in f:
        if not line.startswith(">"):
            continue
        h = line[1:].strip()
        fasta_name = h.split("::")[0]
        m = coord_re.search(h) or fallback_re.search(h)
        if not m:
            bad += 1
            continue
        chrom, start0, end0 = m.group(1), int(m.group(2)), int(m.group(3))
        strand = m.group(4) if m.group(4) else "."
        parts = fasta_name.split("|")
        gene_id = parts[1] if len(parts) > 1 else "NA"
        side = parts[2] if len(parts) > 2 else "NA"
        mtss = tss_re.search(fasta_name)
        tss = mtss.group(1) if mtss else "NA"
        out.write("\t".join(map(str, [chrom, start0, end0, fasta_name, set_name, side, gene_id, tss, strand])) + "\n")
        n += 1
print(f"Wrote {n} BED windows to {out_bed}", file=sys.stderr)
if bad:
    print(f"WARNING skipped {bad} headers without coordinates", file=sys.stderr)
PY

cat "$LOG_DIR/positive_windows_bed.log"
cat "$LOG_DIR/background_windows_bed.log"


# ============================================================
# 4) Python overlap/intersect
# ============================================================

POS_INTERSECT="$INTERSECT_DIR/known_tfbs_positive_windows.tsv"
BG_INTERSECT="$INTERSECT_DIR/known_tfbs_background_windows.tsv"

python3 - "$TFBS_BED" "$POS_WINDOWS_BED" "$BG_WINDOWS_BED" "$POS_INTERSECT" "$BG_INTERSECT" > "$LOG_DIR/python_intersect.log" 2>&1 <<'PY'
import sys, gzip, bisect
from collections import defaultdict

tfbs_path, pos_windows, bg_windows, pos_out, bg_out = sys.argv[1:]

def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def load_windows(path):
    by_chr = defaultdict(list)
    n = 0
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            by_chr[chrom].append((start, end, fields))
            n += 1
    for chrom in by_chr:
        by_chr[chrom].sort(key=lambda x: x[0])
    starts = {chrom: [x[0] for x in arr] for chrom, arr in by_chr.items()}
    return by_chr, starts, n

def overlapping_windows(chrom, hit_start, hit_end, by_chr, starts_by_chr):
    arr = by_chr.get(chrom)
    if not arr:
        return []
    starts = starts_by_chr[chrom]
    i = bisect.bisect_left(starts, hit_end)
    out = []
    j = i - 1
    # Windows are 100 bp in this analysis. Walk backwards until window end <= hit_start.
    while j >= 0 and arr[j][1] > hit_start:
        wstart, wend, wfields = arr[j]
        if wstart < hit_end and wend > hit_start:
            out.append(wfields)
        j -= 1
    return out

pos_by_chr, pos_starts, n_pos_win = load_windows(pos_windows)
bg_by_chr, bg_starts, n_bg_win = load_windows(bg_windows)

n_tfbs = 0
n_pos = 0
n_bg = 0

with open(pos_out, "w") as op, open(bg_out, "w") as ob, open_maybe_gzip(tfbs_path) as tf:
    for line in tf:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            continue

        chrom = fields[0]
        try:
            start = int(fields[1])
            end = int(fields[2])
        except ValueError:
            continue

        n_tfbs += 1

        for w in overlapping_windows(chrom, start, end, pos_by_chr, pos_starts):
            op.write("\t".join(fields + w) + "\n")
            n_pos += 1

        for w in overlapping_windows(chrom, start, end, bg_by_chr, bg_starts):
            ob.write("\t".join(fields + w) + "\n")
            n_bg += 1

print(f"positive windows loaded: {n_pos_win}")
print(f"background windows loaded: {n_bg_win}")
print(f"TFBS rows scanned: {n_tfbs}")
print(f"positive overlaps written: {n_pos}")
print(f"background overlaps written: {n_bg}")

if n_pos == 0:
    print("WARNING: zero positive overlaps written.")
if n_bg == 0:
    print("WARNING: zero background overlaps written.")
PY

cat "$LOG_DIR/python_intersect.log"

echo
echo "Positive intersect rows:"
wc -l "$POS_INTERSECT"

echo "Background intersect rows:"
wc -l "$BG_INTERSECT"

echo
echo "First positive intersect rows:"
head -n 5 "$POS_INTERSECT" || true


# ============================================================
# 5) Manifest
# ============================================================

cat > "$TABLE_DIR/00_known_tfbs_window_overlap_manifest.tsv" <<EOF
parameter	value
PROJECT_DIR	$PROJECT_DIR
RESULTS_DIR	$RESULTS_DIR
INPUT_DIR	$INPUT_DIR
POSITIVE_FASTA	$POS_FASTA
BACKGROUND_FASTA	$BG_FASTA
KNOWN_TFBS_BED	$TFBS_BED
TFBS_NCOL	$TFBS_NCOL
TF_FAMILY_MAP	$FAMILY_MAP
POSITIVE_WINDOWS_BED	$POS_WINDOWS_BED
BACKGROUND_WINDOWS_BED	$BG_WINDOWS_BED
POSITIVE_TFBS_OVERLAPS	$POS_INTERSECT
BACKGROUND_TFBS_OVERLAPS	$BG_INTERSECT
MIN_GAP_BP	10
MAX_GAP_BP	50
MOTIF_SOURCE	PlantTFDB
INTERSECT_METHOD	python_coordinate_overlap
EOF

echo
echo "DONE 01_prepare_known_tfbs_window_overlaps"
echo "Manifest:"
echo "$TABLE_DIR/00_known_tfbs_window_overlap_manifest.tsv"
echo
echo "Next in R/Jupyter:"
echo "Rscript 02_known_tfbs_module_enrichment_and_variant_specific.R"
