#!/usr/bin/env bash
# End-to-end dev run on the rat fixture (BRIEF.md, criterion (c)). Meant for a disposable Linux VM:
# it downloads ~0.5 GB of databases and runs DIAMOND and pyhmmer.
#
#   bash scripts/dev_run_rat.sh <workdir> <threads>
#
# Every stage leaves its output on disk and a marker in <workdir>/markers. A re-run skips finished
# stages, so a crash costs only the stage that was running. Progress: <workdir>/run.log.
set -euo pipefail
W=${1:?workdir}; T=${2:-$(nproc)}
mkdir -p "$W"/{data,markers}
LOG="$W/run.log"
log() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a "$LOG"; }
stage() {  # stage <name> <command...>
  local n=$1; shift
  if [ -e "$W/markers/$n" ]; then log "skip $n (done)"; return; fi
  log "start $n"
  if "$@" >>"$LOG" 2>&1; then touch "$W/markers/$n"; log "done $n"; else log "FAILED $n (rc=$?)"; fi
}

GCA=GCA_036323735.1
REL=https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/10116/GCF_036323735.1-RS_2024_02/
D="$W/data"

stage diamond bash -c "cd '$D' && curl -sSfL -o diamond.tgz \
  https://github.com/bbuchfink/diamond/releases/download/v2.1.10/diamond-linux64.tar.gz && tar -xzf diamond.tgz"
stage carbon gene-evidence fetch carbon --gca $GCA --division vertebrate_mammalian --out "$D/carbon"
stage refseq gene-evidence fetch refseq --release-url $REL --gcf GCF_036323735.1 --out "$D/refseq"
stage swissprot gene-evidence fetch swissprot --out "$D/sprot"
stage pfam gene-evidence fetch pfam --parts 16 --out "$D/pfam"

SPREL=$(grep -o 'Swiss-Prot Release [^ ]* of [^ ]*' "$D/sprot/uniprot_reldate.txt" | head -1)
PFREL="Pfam $(zcat "$D/pfam/Pfam.version.gz" | grep -i 'Pfam release' | awk -F: '{gsub(/ /,"",$2); print $2}')"
args=(analyze --gca $GCA --carbon-gff "$D/carbon/annotations.gff3.gz"
      --carbon-proteins "$D/carbon/proteins.faa.gz" --carbon-metadata "$D/carbon/metadata.json"
      --ref-gff "$D"/refseq/*_genomic.gff.gz --ref-proteins "$D"/refseq/*_protein.faa.gz
      --ref-report "$D"/refseq/*_assembly_report.txt --ref-name GCF_036323735.1-RS_2024_02
      --swissprot "$D/sprot/uniprot_sprot.fasta.gz" --swissprot-release "$SPREL"
      --pfam "$D/pfam/Pfam-A.hmm" --pfam-release "$PFREL"
      --diamond "$D/diamond" --threads "$T" --top 100)

# Two independent runs (separate raw tool tables) for the determinism criterion (a).
stage run1 gene-evidence "${args[@]}" --out "$W/run1"
stage run2 gene-evidence "${args[@]}" --out "$W/run2"

stage compare bash -c "cd '$W' && sha256sum run1/graph.json run2/graph.json run1/report.md run2/report.md \
  > determinism.txt && { cmp run1/graph.json run2/graph.json && cmp run1/report.md run2/report.md \
  && echo IDENTICAL || echo DIFFERENT; } >> determinism.txt"
stage check bash -c "gene-evidence check-report '$W/run1/report.md' '$W/run1/graph.json' > '$W/check.txt'; cat '$W/check.txt'"
log "ALL_DONE"
