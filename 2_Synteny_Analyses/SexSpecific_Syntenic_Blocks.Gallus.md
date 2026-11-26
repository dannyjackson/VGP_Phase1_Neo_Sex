# VGP Sex Chromosome Synteny Analysis

## Run after SexShared_Syntenic_Blocks.md

# Run genespace with the whole genomes 
## Remove X and Z chromosomes from bed files
```
source myconda
mamba activate genespace

BASE="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/YW_Synteny_Gallus_gallus/genespace_input"

mkdir -p ${BASE}
cd ${BASE}
cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input/usable_bed .
cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input/usable_peptide .

INBED="${BASE}/usable_bed"
INPEP="${BASE}/usable_peptide"

OUTBED="${BASE}/usable_beds_noXZ"
OUTPEP="${BASE}/usable_peptides_noXZ"

mkdir -p "$OUTBED" "$OUTPEP"

shopt -s nullglob

# 1) Filter BEDs: drop X/Z (by name AND by accession from the CSV), and
#    rename W/Y accessions to their chromosome name
SEX_CHR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.with_gff.csv"

for bed in "$INBED"/*.bed; do
  sp="$(basename "$bed" .bed)"
  outbed="${OUTBED}/${sp}.bed"

  # X chromosomes named by number (see SexShared md)
  extra=""
  case "$sp" in
    Narcine_bancroftii)    extra=12 ;;
    Scyliorhinus_canicula) extra=28 ;;
  esac

  awk -v sp="$sp" -v extra="$extra" '
    BEGIN{FS=OFS="\t"; if (extra!="") drop[extra]=1}
    FNR==NR {
      split($0, a, ","); gsub(/\r/, "", a[3])
      if (a[1]==sp) {
        if (a[2] ~ /^(X|Z)/) { drop[a[2]]=1; if (a[3]!="") drop[a[3]]=1 }
        else if (a[2] ~ /^(Y|W)/ && a[3]!="" && a[3]!=a[2]) ren[a[3]]=a[2]
      }
      next
    }
    $1 ~ /^(X[0-9]*|Z[0-9]*)$/ {next}
    ($1 in drop) {next}
    { if ($1 in ren) $1=ren[$1]; print }
  ' "$SEX_CHR" "$bed" > "$outbed"

  echo "Wrote $outbed"
done

# 2) Filter peptide FASTAs to match filtered BEDs (keep only genes still present)
#    - assumes peptide file name matches bed base name: SPECIES.bed <-> SPECIES.fa
for bed in "$OUTBED"/*.bed; do
  base="$(basename "$bed" .bed)"
  pep_in="${INPEP}/${base}.fa"
  pep_out="${OUTPEP}/${base}.fa"

  if [[ ! -s "$pep_in" ]]; then
    echo "WARN: missing/empty peptide FASTA: $pep_in (skipping)" >&2
    continue
  fi

  # IDs to keep = col 4 from filtered BED
  keep_ids="$(mktemp)"
  awk -F'\t' '$4!=""{print $4}' "$bed" | sort -u > "$keep_ids"

  # Subset FASTA by header name (requires seqkit)
  seqkit grep -n -f "$keep_ids" "$pep_in" > "$pep_out"

  rm -f "$keep_ids"
  echo "Wrote $pep_out"
done


```
# Make W/Y species lists and check each W/Y is present in its bed
```
REFLISTS=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists
BASE=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/YW_Synteny_Gallus_gallus/genespace_input

# Species, Chr_type, Accession for W/Y rows only
awk -F',' 'NR>1 && $2 ~ /^(Y|W)/ {gsub(/\r/,"",$3); print $1 "\t" $2 "\t" $3}' \
  ${REFLISTS}/sexchrom_accessions.with_gff.csv \
  > ${REFLISTS}/species.sexlimited.all.tsv

# Count genes on each W/Y in the filtered beds
while IFS=$'\t' read -r SP CHR ACC; do
  f="${BASE}/usable_beds_noXZ/${SP}.bed"
  [[ -s "$f" ]] || { printf "%s\t%s\tNO_BED\n" "$SP" "$CHR"; continue; }
  n=$(awk -v c="$CHR" '$1==c' "$f" | wc -l)
  printf "%s\t%s\t%s\n" "$SP" "$CHR" "$n"
done < ${REFLISTS}/species.sexlimited.all.tsv | sort -k3,3n > ${BASE}/WY_gene_counts.tsv

head -40 ${BASE}/WY_gene_counts.tsv   # inspect: 0 / NO_BED = W/Y not in annotated assembly

# Keep only rows with >0 genes on the W/Y -> pair file for plotting (one row per chromosome)
awk -F'\t' '$3 ~ /^[0-9]+$/ && $3 > 0 {print $1 "\t" $2}' ${BASE}/WY_gene_counts.tsv \
  > ${REFLISTS}/species.sexlimited.for_genespace_array.tsv

# Unique species -> list for the GENESPACE runs (one job per species, avoids two jobs in one WORKDIR)
cut -f1 ${REFLISTS}/species.sexlimited.for_genespace_array.tsv | sort -u \
  > ${REFLISTS}/species.sexlimited.unique.txt

wc -l ${REFLISTS}/species.sexlimited.for_genespace_array.tsv ${REFLISTS}/species.sexlimited.unique.txt
```
# Create requisite genespace directory structure
```
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/YW_Synteny_Gallus_gallus/all/genespace

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/YW_Synteny_Gallus_gallus/all/genespace

SEX_CHR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.with_gff.csv"

BED_DIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/YW_Synteny_Gallus_gallus/genespace_input/usable_beds_noXZ"
PEP_DIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/YW_Synteny_Gallus_gallus/genespace_input/usable_peptides_noXZ"

REF_BED="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/Option2/bed/Gallus_gallus.bed"
REF_PEP="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/Option2/peptide/Gallus_gallus.fa"

tail -n +2 "$SEX_CHR" | while IFS=',' read -r SPECIES CHR_TYPE ACC; do
  ACC="${ACC%$'\r'}"

  # Make simple per-species structure
  mkdir -p "${SPECIES}/bed" "${SPECIES}/peptide"

  # Symlink species inputs
  ln -sf "${BED_DIR}/${SPECIES}.bed" "${SPECIES}/bed/${SPECIES}.bed"
  ln -sf "${PEP_DIR}/${SPECIES}.fa"  "${SPECIES}/peptide/${SPECIES}.fa"

  # Symlink references (handy to keep alongside each species)
  ln -sf "$REF_BED" "${SPECIES}/bed/Gallus_gallus_REF.bed"
  ln -sf "$REF_PEP" "${SPECIES}/peptide/Gallus_gallus_REF.fa"
done

```
# Confirm contigs but not scaffolds are kept; confirm X/Z removed
```
awk '{print $1}' Buteo_buteo/bed/Buteo_buteo.bed | sort -u 
awk '{print $1}' Buteo_buteo/bed/Gallus_gallus_REF.bed | sort -u 

grep -cE '^(X|Z)[0-9]*\b' Buteo_buteo/bed/Buteo_buteo.bed   # expect 0

```
# Set up the scripts directory from the GitHub repo
```
YW=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/YW_Synteny_Gallus_gallus
REPO=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/scripts/VGP_Phase1_Neo_Sex

mkdir -p ${YW}/scripts/slurm_output ${YW}/all/plots/slurm_output ${YW}/all/plots_sexchr/slurm_output
mkdir -p "$(dirname ${REPO})"

# clone once; afterwards just pull
if [[ -d ${REPO}/.git ]]; then
  git -C ${REPO} pull --ff-only
else
  git clone https://github.com/dannyjackson/VGP_Phase1_Neo_Sex.git ${REPO}
fi

# record which commit this run used
git -C ${REPO} rev-parse HEAD > ${YW}/scripts/repo_commit.txt
cat ${YW}/scripts/repo_commit.txt

# link the R scripts into the scripts dir (sbatch SCRIPT= paths stay unchanged)
ln -sf ${REPO}/2_Synteny_Analyses/genespace.synteny_to_chicken.wholegenome.r \
       ${YW}/scripts/genespace.synteny_to_chicken.wholegenome.r
ln -sf ${REPO}/2_Synteny_Analyses/genespace.synteny_to_chicken.wholegenome.sexspecific.r \
       ${YW}/scripts/genespace.synteny_to_chicken.wholegenome.sexspecific.r

# confirm the Part A edits are present in the clone
grep -n 'png_name <-' ${YW}/scripts/genespace.synteny_to_chicken.wholegenome.r     # expect a match
grep -n '\.38\.pdf'   ${YW}/scripts/genespace.synteny_to_chicken.wholegenome.sexspecific.r  # expect no match
```
# Submission script
```

#!/bin/bash
#SBATCH --job-name=gs_chicken_wg
#SBATCH --output=slurm_output/%x_%A_%a.out
#SBATCH --error=slurm_output/%x_%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

set -euo pipefail

PAIRFILE=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexlimited.unique.txt
SCRIPT=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/YW_Synteny_Gallus_gallus/scripts/genespace.synteny_to_chicken.wholegenome.r
OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/YW_Synteny_Gallus_gallus/all/plots

# Grab the Nth species (array index) and parse "Species Chromosome"
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PAIRFILE}")
SPECIES=$(awk '{print $1}' <<< "${line}")
WORKDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/YW_Synteny_Gallus_gallus/all/genespace/${SPECIES}

export WORKING_DIR="${WORKDIR}"
export GS_BLKSIZE=5   # set to 3 if W/Y form no blocks

cd "${OUTDIR}"

source myconda
mamba activate genespace

Rscript "${SCRIPT}" "${SPECIES}" "${OUTDIR}/"

# Reduce file number generated by genespace
# make globs that match nothing expand to nothing
shopt -s nullglob

# tar+bgzip each Orthogroup_Sequences dir, then delete it
for r in "${WORKDIR}"/orthofinder/Results*/Orthogroup_Sequences; do
  [ -d "$r" ] || continue
  out="${r%/}.tar.gz"
  tar -C "$(dirname "$r")" -cf - "$(basename "$r")" | bgzip -c > "$out" && rm -rf "$r"
done

# 2) One combined tar+bgzip of ALL Single_Copy_Orthologue_Sequences (across Results*), written to ${OUTDIR}/${SPECIES}/
for r in "${WORKDIR}"/orthofinder/Results*/Single_Copy_Orthologue_Sequences; do
  [ -d "$r" ] || continue
  out="${r%/}.tar.gz"
  tar -C "$(dirname "$r")" -cf - "$(basename "$r")" | bgzip -c > "$out" && rm -rf "$r"
done

```
## Submit batch job
```
# create output directory for plots
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/YW_Synteny_Gallus_gallus/all/plots

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/YW_Synteny_Gallus_gallus/scripts

N=$(wc -l < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexlimited.unique.tsv)

# test
sbatch --array=1-3 genespace_to_chicken_array.wholegenome.sh

# finish
sbatch --array=4-${N} genespace_to_chicken_array.wholegenome.sh

```
## Plot just synteny with sex chrom
```
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/YW_Synteny_Gallus_gallus/all/plots_sexchr

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/YW_Synteny_Gallus_gallus/scripts

N=$(wc -l < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexlimited.for_genespace_array.tsv)
sbatch --array=1-${N} genespace_to_chicken_array.wholegenome.sexspecific.sh

```
### Sbatch script
```
#!/bin/bash
#SBATCH --job-name=gs_chicken_wg
#SBATCH --output=slurm_output/%x_%A_%a.out
#SBATCH --error=slurm_output/%x_%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

set -euo pipefail

PAIRFILE=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexlimited.for_genespace_array.tsv
SCRIPT=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/YW_Synteny_Gallus_gallus/scripts/genespace.synteny_to_chicken.wholegenome.sexspecific.r
OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/YW_Synteny_Gallus_gallus/all/plots_sexchr

# Grab the Nth species (array index) and parse "Species Chromosome"
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PAIRFILE}")
SPECIES=$(awk '{print $1}' <<< "${line}")
SEXCHR=$(awk '{print $2}' <<< "${line}")
WORKDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/YW_Synteny_Gallus_gallus/all/genespace/${SPECIES}

export WORKING_DIR="${WORKDIR}"

cd "${OUTDIR}"

source myconda
mamba activate genespace

Rscript "${SCRIPT}" "${SPECIES}" "${OUTDIR}/" ${SEXCHR}
```


## Create a table of syntenic regions to chicken
### Sex-shared to chicken WG
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/YW_Synteny_Gallus_gallus/all/

```
# Combine all *_phasedBlks.csv files into one CSV (one header), adding provenance columns.

library(readr)
library(dplyr)
library(purrr)
library(stringr)

roots <- c("genespace")

files <- unlist(lapply(roots, function(r) {
  if (dir.exists(r)) list.files(r, pattern = "_phasedBlks\\.csv$", recursive = TRUE, full.names = TRUE) else character()
}), use.names = FALSE) %>% sort()

stopifnot(length(files) > 0)

parse_meta <- function(path) {
  # Example: Z/Accipiter_gentilis/riparian/Accipiter_gentilis_phasedBlks.csv
  parts <- strsplit(path, .Platform$file.sep, fixed = TRUE)[[1]]
  tibble(
    topDir = parts[1] %||% NA_character_,
    species = parts[2] %||% NA_character_,
    subdir = parts[3] %||% NA_character_,
    sourceFile = path
  )
}

# Read all, keep union of columns (fills missing with NA), add provenance cols

df <- map_dfr(files, function(f) {
  dat <- read_csv(
    f,
    col_types = cols(.default = col_character()),
    progress = FALSE,
    show_col_types = FALSE
  )
  bind_cols(parse_meta(f), dat)
})

write_csv(df, "all_phasedBlks_combined.csv")
cat("Wrote all_phasedBlks_combined.csv with", nrow(df), "rows from", length(files), "files\n")

# rename instances of "CP162267.1" to Y
df$chr1[df$chr1 == "CP162267.1"] <- "Y"
df$chr2[df$chr2 == "CP162267.1"] <- "Y"
  
# Filter to just rows where a sex chromosome is in either chr 1 or chr 2
df2 <- df %>%
  filter(grepl("[W-Yw-y]", chr1) | grepl("[W-Yw-y]", chr2))

# remove rows that are chicken Z > other species 
df3 <- df2 %>%
  filter(!(genome1 == "Gallus_gallus_REF" & chr1 == "Z" & !grepl("^(Y|W)", chr2))) %>%
  filter(!(genome2 == "Gallus_gallus_REF" & chr2 == "Z" & !grepl("^(Y|W)", chr1)))

VGP_DATAFRAME <- read_csv("/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/VGP_OrdinalList_Phase1Freeze_v1.1_sex_chroms_seqCov_HalfDeep_SCINKD_May8.25_DwnldDec16.25.csv")
# rename Canis lupus baileyi to Canis lupus
VGP_DATAFRAME$Scientific.Name[VGP_DATAFRAME$Scientific.Name == "Canis lupus baileyi"] <- "Canis lupus"

# rename to Guaruba_guaruba to Guaruba guarouba in df
df3$species[df3$species == "Guaruba_guaruba"] <- "Guaruba_guarouba"

# Join df3 and VGP_DATAFRAME on df3$species and VGP_DATAFRAME$Scientific.Name for columns Lineage, Superorder


df4 <- df3 %>%
  mutate(species_join = gsub("_", " ", species)) %>%
  left_join(
    VGP_DATAFRAME %>% 
      transmute(Scientific.Name, Lineage, Superorder,
                Scientific_join = trimws(Scientific.Name)),
    by = c("species_join" = "Scientific_join")
  ) %>%
  select(-species_join)

unmatched <- df4 %>% filter(is.na(Lineage) & is.na(Superorder)) %>% distinct(species)
unmatched


write_csv(df4, "sexchrs.phasedBlks.csv")
unique_chrs <- unique(df4[c("Lineage", "Superorder", "genome1","genome2","chr1","chr2")])
write_csv(unique_chrs, "unique_chrs.csv")

unique_chrs <- read.csv("unique_chrs.csv")

# condense to just one instance of matching, retaining info if a chr is only matched once not twice
gg_name <- "Gallus_gallus_REF"   # change if your column uses a different exact string

out <- unique_chrs %>%
  # keep rows where chicken is in exactly one of the genome columns
  filter(xor(genome1 == gg_name, genome2 == gg_name)) %>%
  mutate(
    Species = if_else(genome1 == gg_name, genome2, genome1),
    SpChr   = if_else(genome1 == gg_name, chr2,   chr1),

    # chicken chromosome when chicken is genome1 -> chicken chr is chr1
    gg_chr1 = if_else(genome1 == gg_name, chr1, NA_character_),

    # chicken chromosome when chicken is genome2 -> chicken chr is chr2
    gg_chr2 = if_else(genome2 == gg_name, chr2, NA_character_)
  ) %>%
  group_by(Species, SpChr) %>%
  summarise(
    gg_chr1 = list(sort(unique(na.omit(gg_chr1)))),
    gg_chr2 = list(sort(unique(na.omit(gg_chr2)))),

    GgChr_both      = list(intersect(gg_chr1[[1]], gg_chr2[[1]])),
    GgChr_chr1_only = list(setdiff(gg_chr1[[1]], gg_chr2[[1]])),
    GgChr_chr2_only = list(setdiff(gg_chr2[[1]], gg_chr1[[1]])),
    .groups = "drop"
  ) %>%
  select(Species, SpChr, GgChr_both, GgChr_chr1_only, GgChr_chr2_only)


out_csv <- out %>%
  mutate(across(starts_with("GgChr"),
                ~ vapply(.x, \(v) paste(v, collapse = ","), character(1))))

readr::write_csv(out_csv, "syntentic_chromosomes.sex_shared.Gallus_gallus_REF.csv")


# Plot it on a phylogeny


# Plot it on the phylogeny
## Blanks for empty?

library(data.table)
library(ggplot2)
library(patchwork)
library(ape)
library(dplyr)
library(tidyr)
library(purrr)


out <- read.csv("syntentic_chromosomes.sex_shared.Gallus_gallus_REF.csv")

tree_file <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/roadies_v1.1.16b.numbers.scientific.nwk"

tr_full <- ape::read.tree(tree_file)

tr_full$tip.label[tr_full$tip.label == "Accipiter_gentilis"] <- "Astur_gentilis"
tr_full$tip.label <- sub("_$", "", tr_full$tip.label)

canis_subspp <- grep("^Canis_lupus_", tr_full$tip.label, value = TRUE)
if (length(canis_subspp) > 0) {
  keep_tip <- if ("Canis_lupus_baileyi" %in% canis_subspp) "Canis_lupus_baileyi" else canis_subspp[1]
  tr_full <- drop.tip(tr_full, setdiff(canis_subspp, keep_tip))
  tr_full$tip.label[tr_full$tip.label == keep_tip] <- "Canis_lupus"
}
stopifnot(!anyDuplicated(tr_full$tip.label))

# species in the table that are not on the tree (fix names if any print)
setdiff(unique(out$Species), tr_full$tip.label)

tr_ultra <- chronos(tr_full, lambda = 1)


# out: your tibble with Species and list-cols
# tr_ultra: your ultrametric tree

# union the three list-cols

out_long <- out %>%
  mutate(
    gg_all = map2(
      map2(GgChr_both, GgChr_chr1_only, ~ union(.x, .y)),
      GgChr_chr2_only,
      ~ union(.x, .y)
    )
  ) %>%
  select(Species, gg_all) %>%
  unnest(gg_all) %>%
  filter(!is.na(gg_all), gg_all != "") %>%
  distinct(Species, gg_all) %>%
  rename(GgChr = gg_all) %>%
  # --- NEW: split combos like "22,Z" into separate rows "22" and "Z"
  separate_rows(GgChr, sep = ",") %>%
  mutate(GgChr = str_trim(GgChr)) %>%
  filter(GgChr != "") %>%
  distinct(Species, GgChr)

# order chicken chromosomes: numeric first, then others (Z/W/etc)
chr_levels <- out_long %>%
  distinct(GgChr) %>%
  mutate(
    is_num = suppressWarnings(!is.na(as.integer(GgChr))),
    numval = suppressWarnings(as.integer(GgChr))
  ) %>%
  arrange(desc(is_num), numval, GgChr) %>%
  pull(GgChr)

# full grid for all tips x all chrs
grid_df <- tidyr::expand_grid(
  Species = out_long$Species,
  GgChr   = chr_levels
) %>%
  left_join(out_long %>% mutate(present = 1L), by = c("Species", "GgChr")) %>%
  mutate(present = if_else(is.na(present), 0L, present))





# Plot tree

pdf("tree.sexspecific.Gallusgallus.pdf", width = 3, height = 10)

# give more right margin for the dot grid + top margin for labels
par(mar = c(2, 2, 6, 12))

# 1) Compute how much horizontal space you need for the dot grid
#    (dx * nchr) plus some padding.
tmp <- plot(tr_ultra, plot = FALSE)
nchr <- length(chr_levels)

# spacing between chromosome columns
dx <- max(tmp$x.lim) * 0.02

# extra horizontal room to the right of the tree (in tree x-units)
grid_width <- dx * (nchr + 6)

# 2) Plot tree but "squish it left" by expanding xlim to the right
plot(tr_ultra,
     cex = 0.1,
     edge.width = 0.25,
     no.margin = TRUE,
     x.lim = c(0, max(tmp$x.lim) + grid_width))

lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_y <- lp$yy[1:Ntip(tr_ultra)]
tips  <- tr_ultra$tip.label

# start x position for the grid (a bit to the right of the tree)
x0 <- max(tmp$x.lim) + dx * 3

# full tips/y from the plotted tree
lp   <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tips <- tr_ultra$tip.label
tip_y <- lp$yy[1:Ntip(tr_ultra)]

# subset: tips that appear in grid_df
tips_with_data <- intersect(tips, unique(grid_df$Species))
idx_data <- match(tips_with_data, tips)
y_data <- tip_y[idx_data]

# within each chromosome column:
for (j in seq_along(chr_levels)) {
  chrj <- chr_levels[j]
  xj <- x0 + (j - 1) * dx

  present_species <- grid_df %>%
    dplyr::filter(GgChr == chrj, present == 1L) %>%
    dplyr::pull(Species)

  is_present_data <- tips_with_data %in% present_species

  # light-gray circles ONLY for species present in grid_df
  points(rep(xj, length(y_data)), y_data,
         pch = 21, cex = 0.05, col = "lightgray", bg = "white")

  # black filled circles for present==1 among those species
  points(rep(xj, sum(is_present_data)), y_data[is_present_data],
         pch = 21, cex = 0.05, bg = "black")
}

# 2) Column labels at the top
usr <- par("usr")

# Put labels slightly ABOVE the top of the plotting region
y_lab <- 585
text(
  x = x0 + (seq_along(chr_levels) - 1) * dx,
  y = y_lab,
  labels = chr_levels,
  srt = 90,                 # rotate
  adj = c(0, 0.5),          # anchor at bottom of rotated text
  xpd = NA,                 # allow drawing into margin area
  cex = 0.2                 # tweak as needed
)

dev.off()




































df <- df3
gallus <- "Gallus_gallus_REF"

# Helper to safely compute bp length
bp_len <- function(start, end) {
  s <- suppressWarnings(as.double(start))
  e <- suppressWarnings(as.double(end))
  abs(e - s) + 1
}



bp_len <- function(start, end) {
  s <- suppressWarnings(as.double(start))
  e <- suppressWarnings(as.double(end))
  abs(e - s) + 1
}

# VGP lookup table: "Genus species" -> "Genus_species"
vgp_lut <- VGP_DATAFRAME %>%
  transmute(
    TaxonKey = str_replace_all(str_squish(Scientific.Name), " ", "_"),
    Lineage,
    Superorder
  ) %>%
  distinct(TaxonKey, .keep_all = TRUE)


summary_df <- df %>%
  # only Gallus vs non-Gallus comparisons
  filter((genome1 == gallus & genome2 != gallus) | (genome2 == gallus & genome1 != gallus)) %>%
  mutate(
    Species = if_else(genome1 == gallus, genome2, genome1),
    Chromosome = if_else(genome1 == gallus, chr1, chr2),

    BP_Gallus = if_else(genome1 == gallus,
                        bp_len(startBp1, endBp1),
                        bp_len(startBp2, endBp2)),
    BP_Query = if_else(genome1 == gallus,
                       bp_len(startBp2, endBp2),
                       bp_len(startBp1, endBp1)),

    nHits_Gallus = suppressWarnings(as.double(if_else(genome1 == gallus, nHits1, nHits2))),
    nHits_Query  = suppressWarnings(as.double(if_else(genome1 == gallus, nHits2, nHits1))),

    # join key for VGP (underscored binomial)
    TaxonKey = Species
  ) %>%
  # drop reciprocal duplicates (optional; remove if you want both directions counted)
  distinct(Species, Chromosome, blkID, startBp1, endBp1, startBp2, endBp2, nHits1, nHits2, .keep_all = TRUE) %>%
  group_by(Species, Chromosome, TaxonKey) %>%
  summarise(
    TotalBP_Gallus = sum(BP_Gallus, na.rm = TRUE),
    nHits_Gallus   = sum(nHits_Gallus, na.rm = TRUE),
    TotalBP_Query  = sum(BP_Query, na.rm = TRUE),
    nHits_Query    = sum(nHits_Query, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(vgp_lut, by = "TaxonKey") %>%
  relocate(Lineage, Superorder, .after = Species) %>%
  select(Species, Lineage, Superorder, Chromosome,
         TotalBP_Gallus, nHits_Gallus, TotalBP_Query, nHits_Query) %>%
  arrange(Lineage, Superorder, Species, Chromosome)


# Optional: show which Species didn't match VGP after stripping _REF
unmatched <- summary_df %>%
  filter(is.na(Lineage) | is.na(Superorder)) %>%
  distinct(Species)
if (nrow(unmatched) > 0) {
  cat("\nUnmatched Species (after stripping _REF for lookup):\n")
  print(unmatched)
}
# Umatched:
1 Canis_lupus    
2 Guaruba_guaruba

summary_df <- summary_df %>%
  mutate(
    Lineage = case_when(
      Species == "Canis_lupus" ~ "Mammal",
      Species == "Guaruba_guaruba" ~ "Bird",
      TRUE ~ Lineage
    )
  )

summary_df <- summary_df %>%
  mutate(
    Superorder = case_when(
      Species == "Canis_lupus" ~ "Laurasiatheria",
      Species == "Guaruba_guaruba" ~ "Psittaciformes",
      TRUE ~ Superorder
    )
  )

write_csv(summary_df, "Gallus_vs_species_window_sums_with_lineage.sex_limited.csv")
cat("Wrote Gallus_vs_species_window_sums_with_lineage.sex_limited.csv with", nrow(summary_df), "rows\n")

























```
# Make a plot NOT WORKING YET
```
# Plot synteny blocks on Gallus Z (0..86,044,486 bp) with species on Y

library(readr)
library(dplyr)
library(purrr)
library(stringr)

library(ggplot2)

gallus <- "Gallus_gallus_REF"
chr_target <- "Z"
Z_len <- 86044486

to_num <- function(x) suppressWarnings(as.numeric(x))

synteny_plot_df <- df4 %>%
  # keep Gallus vs non-Gallus
  filter((genome1 == gallus & genome2 != gallus) | (genome2 == gallus & genome1 != gallus)) %>%
  mutate(
    Species = if_else(genome1 == gallus, genome2, genome1),

    # Gallus chromosome for this row
    GallusChr = if_else(genome1 == gallus, chr1, chr2),

    # Gallus interval for this row
    g_start = if_else(genome1 == gallus, startBp1, startBp2),
    g_end   = if_else(genome1 == gallus, endBp1,   endBp2),

    g_start = to_num(g_start),
    g_end   = to_num(g_end),

    x_start = pmin(g_start, g_end, na.rm = TRUE),
    x_end   = pmax(g_start, g_end, na.rm = TRUE)
  ) %>%
  # only Gallus Z
  filter(as.character(GallusChr) == chr_target) %>%
  # drop reciprocal duplicates (optional but usually desirable)
  distinct(Species, GallusChr, blkID, x_start, x_end, .keep_all = TRUE) %>%
  filter(!is.na(x_start), !is.na(x_end)) %>%
  mutate(Species = factor(Species, levels = sort(unique(as.character(Species)))))

# Ensure A->Z (top to bottom depends on your preference; this is the common "A at top")
synteny_plot_df <- synteny_plot_df %>%
  mutate(Species = factor(as.character(Species), levels = sort(unique(as.character(Species)))))

p <- ggplot(synteny_plot_df, aes(y = Species)) +
  geom_segment(aes(x = x_start, xend = x_end, yend = Species), linewidth = 0.4, alpha = 0.8) +
  coord_cartesian(xlim = c(0, Z_len), clip = "off") +
  scale_y_discrete(limits = levels(synteny_plot_df$Species)) +  # enforce the order
  scale_x_continuous(labels = scales::comma) +
  labs(
    x = "Gallus_gallus Z position (bp)",
    y = "Species",
    title = "Synteny blocks to Gallus_gallus Z"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 7),
    plot.margin = margin(t = 15, r = 10, b = 15, l = 10)  # <- THIS fixes top/bottom cutoff
  )

p <- p +
  geom_vline(xintercept = 0,  linewidth = 0.2, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = Z_len, linewidth = 0.2, linetype = "dashed", colour = "grey60")

# give the output a touch more vertical room too (optional but helps)
ggsave(
  "synteny_to_gallus_Z.sex_limited.png",
  p,
  width = 10,
  height = max(4, 0.22 * nlevels(synteny_plot_df$Species)), # slightly larger multiplier
  dpi = 300
)


ggsave("synteny_to_gallus_Z.sex_limited.png", p, width = 10, height = max(4, 0.2 * nlevels(synteny_plot_df$Species)), dpi = 300)
```