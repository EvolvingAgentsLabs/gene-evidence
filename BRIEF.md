# BRIEF — gene-evidence, stage 1 (MVP)

Written 2026-10-08, before any code or run. Nothing in this file has run yet. Claims are marked
**[read]** (inferred from documentation or source) or **[ran]** (observed). Sections added after a run
are appended under their own dated heading; the text above them is not edited after the fact.

## The problem

Carbon-A ([`HuggingFaceBio/genbank-annotations`](https://huggingface.co/buckets/HuggingFaceBio/genbank-annotations),
MIT) publishes protein-coding gene predictions for thousands of GenBank assemblies: a GFF3, the CDS and the
protein of every predicted locus, and a per-locus `confidence` and `status` **[read]**. Some of those loci
are absent from the species' reference annotation (NCBI RefSeq). Each of them is either a gene the
reference missed or a prediction artefact. A biologist who wants to follow one up has to assemble, by
hand, what the reference says at that locus, whether the protein looks like a known one, whether it has
domains, and whether it looks like a copy of another gene. Then they have to choose an experiment.

Stage 1 does that assembly mechanically and makes every statement auditable. It does not decide which
candidates are real.

## What the MVP is

`gene-evidence analyze` — one command, no LLM, no learned model.

**Inputs**

- A GenBank assembly accession (GCA) with Carbon-A output: `annotations.gff3.gz` and `proteins.faa.gz`.
- A reference annotation for the same assembly, **or none**: the RefSeq `*_genomic.gff.gz`, the matching
  `*_protein.faa.gz` and `*_assembly_report.txt` of the RefSeq twin (GCF). GenBank sequence names map
  to RefSeq names through the assembly report rows whose relationship is `=` (identical sequence). Loci
  on a sequence without an identical twin are dropped and counted.
- Databases: UniProtKB/Swiss-Prot (FASTA) and Pfam-A (HMM library).
- Optional: RNA evidence as a transcript GTF (see evidence 6).

**Candidates.** Carbon-A protein-coding loci whose CDS (union of CDS segments) overlaps **no**
protein-coding CDS of the reference on the same strand. A reference CDS marked `pseudo=true` is not
protein-coding. Without a reference, every Carbon-A locus is a candidate, and the report says so.

**Outputs** (one directory)

- `graph.json`, the **evidence graph**. Every node is one tool output. It records the tool, the tool
  version, the database name, version and sha256 (if a database was used), the parameters, and the raw
  values. Edges link each candidate to its evidence nodes and to the rule nodes that turned that
  evidence into a tier and an experiment.
- `report.md`: a ranked summary table plus one section per candidate. Each section gives the evidence
  FOR and the evidence AGAINST the candidate being a real, novel protein-coding gene, and a suggested
  validation experiment chosen by fixed rules. **Every sentence and every table row cites graph nodes**
  as `[node-id]`.
- `versions.json`: tools, databases and inputs, with their sha256. `timings.json`: wall time and peak RAM.
  Timings live only here, so the graph and the report contain no clock values.

## Evidence, in this order

The order puts locus context first on purpose. Gene predictors often call processed pseudogenes
(retrocopies) as genes, and the reference usually annotates those loci as pseudogenes, so many candidates
are expected to sit on them. On such a locus a near-identical homolog marks a **retrocopy**, not a new
gene, so homology alone points the wrong way. That expectation is a hypothesis here, not a result. The
dev run below measures the pseudogene fraction on its own fixture **[to run]**.

1. **Locus context.** Which reference `gene`/`pseudogene` features overlap the candidate's span, on
   either strand, and their `gene_biotype`. Each candidate gets one class, checked in this order:
   *pseudogene* (any overlapping biotype containing `pseudo`); *protein-coding gene* (the span overlaps
   a protein-coding gene, which is possible on the other strand or inside an intron, because CDS
   overlap is already excluded); *non-coding gene* (any other biotype); *nothing*. Also recorded: the
   distance in bp to the nearest reference protein-coding gene on the same sequence, either strand,
   0 if they overlap.
2. **Retrocopy signature.** Inputs: the candidate's CDS segment count, and its best DIAMOND blastp hit
   against the **reference's own proteome** (the best paralog; best = highest bitscore, ties broken by
   subject id). **Rule (fixed):** retrocopy signature iff the candidate has exactly **1** CDS segment
   **and** the best paralog hit has identity **≥ 90 %** and query coverage **≥ 80 %**, **and** the
   paralog's parent CDS has **≥ 2** segments, **and** the parent gene does not overlap the candidate.
   Also recorded: a `*` inside the protein (internal stop) and Carbon-A's `status`. Carbon-A reports
   only `complete`/`incomplete` and publishes no frameshift flag **[read]**, so "frameshift" is not
   evidence this stage can have, and the report says so.
3. **Homology.** Best DIAMOND blastp hit against Swiss-Prot (e ≤ 1e-5): identity, query coverage,
   e-value, bitscore, the subject's organism and taxon id. A **homolog** for tiering means e ≤ 1e-5
   **and** query coverage ≥ 50 %.
4. **Protein domains.** pyhmmer `hmmsearch` of Pfam-A against the candidate proteins, with gathering
   thresholds (`--cut_ga`). A domain is an *included* domain hit.
5. **ORF sanity.** Starts with `M`; protein length; Carbon-A `confidence` and `status` as reported. These
   are **labelled "the predictor's own score", never as evidence**. ORF flags (no leading `M`, internal
   stop, `status` ≠ `complete`) count as weak evidence against.
6. **RNA evidence (optional).** If a transcript GTF is supplied, the report records the transcripts
   overlapping the candidate on the same strand, and whether a transcript intron matches a predicted
   splice junction exactly. Without the GTF, this evidence is a recorded no-op node, `rna: not
   supplied`. It is not used for tiering in stage 1. A BAM input is out of scope for stage 1.

## Ranking — a reading aid, not a claim

The tiers below are a fixed, documented rule, written here before any output was seen. Nothing in
them is learned.

| tier | name | rule (first match wins) |
|---|---|---|
| 5 | AGAINST-strong | locus class is *pseudogene*, **or** the retrocopy signature fired. This overrides any FOR evidence. |
| 1 | FOR-strong | Swiss-Prot homolog **and** ≥ 1 Pfam domain |
| 2 | FOR | Swiss-Prot homolog **or** ≥ 1 Pfam domain |
| 4 | AGAINST-weak | any ORF flag, and no FOR evidence |
| 3 | no evidence | everything else |

Order within a tier: number of distinct Pfam families (descending), then best Swiss-Prot bitscore
(descending, 0 if no hit), then locus id (ascending). Carbon-A confidence is **not** used.

**Stated plainly:** no ranking of such candidates has been validated. The tiers sort the reading list so that the evidence against (pseudogene, retrocopy) is seen
first. They do not claim that tier 1 is enriched for real genes.

## Validation experiment — fixed rules (first match wins)

1. Tier 5 with a retrocopy signature or a paralog hit: **paralog-discriminating RT-PCR**, with primers on
   positions where the candidate differs from its closest reference paralog. The report cites the
   paralog and its identity.
2. Tier 5 without a paralog hit: **long-read cDNA sequencing** across the locus, to see whether an
   intact ORF is transcribed at all.
3. ≥ 2 CDS segments: **RT-PCR across the first predicted splice junction**, with primers in the two
   flanking exons. The report cites the junction coordinates and the expected intron size.
4. One CDS segment, protein < 100 aa: **ribosome profiling** (Ribo-seq), for translation of a small ORF.
5. One CDS segment, protein ≥ 100 aa: **long-read cDNA sequencing**.

## Acceptance criteria (fixed now)

| id | criterion | how it is checked |
|---|---|---|
| (a) | **Deterministic:** two runs on the same inputs give byte-identical `graph.json` and `report.md` | the e2e run is done twice into separate directories; sha256 compared. A unit test also does it on synthetic fixtures. |
| (b) | **Traceable:** every report sentence and table row carries ≥ 1 graph-node reference, and every reference resolves | `gene-evidence check-report report.md graph.json` exits non-zero otherwise. A unit test feeds it an uncited sentence and a dangling reference, and **asserts that it fails**. The e2e report must pass it. |
| (c) | **Runs end-to-end** on the dev fixture with **10–100** candidates in the output | the dev run below. Evidence is computed for every candidate, and `--top 100` keeps the 100 highest-ranked in the graph and report. |
| (d) | **Unit tests** with tiny synthetic fixtures, no network: overlap, retrocopy rule, tiering, report-claim checker | `pytest`, in CI |
| (e) | **Every database version and sha256** recorded in the output | `versions.json`, plus database nodes in `graph.json` |
| (f) | **Wall time and peak RAM** of the e2e run reported | `timings.json`, plus a section appended to this file |

**Dev fixture:** *Rattus norvegicus*. Carbon-A on GCA_036323735.1, with the RefSeq twin
GCF_036323735.1, annotation release **RS_2024_02**, as the reference. It is a **development fixture**.
It checks that the pipeline runs and that its output is auditable. It is **not a test of ranking
quality**, and no number from it is quoted as one.

## What stage 1 does not claim

- That the ranking finds real genes. That needs expert review and RNA evidence, which come later.
- That a tier-5 candidate is not a gene. Retrocopies can be expressed, and some pseudogene calls are
  wrong.
- That "no Swiss-Prot hit" means "novel". Swiss-Prot is small and curated.

## Compute and cost

- Unit tests are stdlib plus pytest, with synthetic fixtures. They may run anywhere.
- The e2e run needs Swiss-Prot (~90 MB gz), Pfam-A (~300 MB gz), DIAMOND and pyhmmer. It runs on **one
  Colab session**: a CPU VM, or an L4 VM used only for its CPUs. The run is detached, each stage
  persists its output as it lands, and a log is polled. If the session dies, the next one resumes from
  the persisted stages. **Ceiling: 2 Colab sessions** for this stage. Downloads and tool runs happen
  only on the VM.
- No API spend.

## Stopping condition (written before the run)

If the e2e run fails (a) or (b), the failure is fixed in code and the run is repeated once, within the
2-session ceiling. If it fails again, stage 1 is reported as not passing that criterion, with the
diff. The tier rule and the thresholds above are **not** changed after seeing the dev output. Any
change goes into a dated section below as a redesign, with the reason, and is counted.

Redesign count: 0.

## Result of the dev run (2026-10-09) [ran]

One Colab session, `ge-rat`: a CPU runtime with 2 vCPUs and 12 GB RAM, about 16 minutes from creation to
stop. `scripts/dev_run_rat.sh` ran detached, each stage left a marker, and the log was polled.
Downloads were verified on the VM. The Carbon-A files matched the sha256 values in Carbon-A's own
`metadata.json`. The NCBI files matched the release's `md5checksums.txt`. The Pfam-A file passed an
exact byte count and a full gzip read. The package wheel and the run script were uploaded and checked
by sha256 on the VM before use. Outputs: `runs/2026-10-09-rat-dev/`. The tool tables and DIAMOND
databases stayed on the VM.

This is a re-run of the first dev run of 2026-10-08 (kept in git history) with the current package. Against that run, `graph.json`, `report.md`, `versions.json` and `all_candidates.tsv` differ only in the tool name; the candidate set, tiers, ranking, database versions and input sha256 values are identical. The `graph.json` and `report.md` hashes change with the name; the timings differ by VM.

| id | criterion | result |
|---|---|---|
| (a) | deterministic | **PASS**. Two full independent runs, each with its own DIAMOND and hmmsearch: `graph.json` sha256 `c1d7df17…8a0a7` and `report.md` sha256 `4a959c7b…3055b` are identical (`determinism.txt`). The unit test does the same on synthetic fixtures. |
| (b) | traceable | **PASS**. `gene-evidence check-report` passes on the committed report (`check.txt`). Unit tests show that it **fails** on an uncited sentence, a dangling reference, an uncited table row and an appended uncited claim. |
| (c) | end-to-end, 10–100 candidates | **PASS**. 947 candidates were analysed and the top 100 were kept in the graph and report. A `--top 100000` re-render from the same tool tables also passes the checker. Its compact per-candidate table is `all_candidates.tsv`. |
| (d) | unit tests, no network | **PASS**. 35 tests on overlap, the candidate rule, locus context, the retrocopy rule, tiers, experiments, the checker, determinism, no-reference mode and the RNA no-op. CI runs them on Python 3.9 and 3.12. |
| (e) | versions and sha256 recorded | **PASS**. `versions.json` and the `db:*` nodes: Swiss-Prot Release 2026_03 of 02-Sep-2026 (`uniprot_sprot.fasta.gz` `a9c3496a…4169536`), Pfam 38.2 (`Pfam-A.hmm` `4b0da639…66739e`; the `.gz` is `2d82087b…5dab`, in `run.log`), RefSeq GCF_036323735.1-RS_2024_02 GFF/proteome/report, Carbon-A GFF/proteins (post-processing revision 6ffa159), DIAMOND 2.1.10, pyhmmer 0.12.3. |
| (f) | wall time and peak RAM | **Reported**. `gene-evidence analyze` took 364 s wall on 2 vCPUs: inputs 21 s, DIAMOND vs Swiss-Prot 92 s (including makedb), DIAMOND vs the reference proteome 29 s, hmmsearch vs Pfam-A 211 s, graph and report 2 s. Peak RSS was 480 MB in the Python process and 797 MB in the largest child (DIAMOND). Downloads took about 97 s. The second run took 363 s. |

**Dev fixture numbers.** These describe the fixture, not ranking quality.

- 22,337 Carbon-A loci. 21,390 overlap a RefSeq protein-coding CDS on the same strand, 0 were dropped
  for lack of an identical sequence, and **947 are candidates**.
- Locus context: **pseudogene 578 (61 %)**, nothing 285, protein-coding gene (other strand or intron)
  45, non-coding gene 39. This re-measures, on this fixture, the observation that motivated the
  evidence order.
- The retrocopy signature fired on 261 candidates, 27 of them on loci the reference does not call
  pseudogenes.
- Tiers: 1 FOR-strong 223, 2 FOR 75, 3 no evidence 44, 4 AGAINST-weak **0**, 5 AGAINST-strong 605.
  Tier 4 is empty because every candidate starts with `M`, has `status=complete`, and Carbon-A's
  proteins carry no `*`.
- Experiments over the top 100: long-read cDNA 53, RT-PCR across a junction 47.

**What reading the top 100 showed** [ran: counted from `all_candidates.tsv`]. All of the top 100 are
tier 1, and they are dominated by repeat-rich families: 32 have a best hit to the Smok kinase family
(`SMKY_MOUSE`), 24 hit T-cell receptor or immunoglobulin variable segments, and 2, including
rank 1, hit retroviral Pol polyproteins. That is the transposon / multi-copy failure mode written into
`docs/EVIDENCE.md`: domains and homology reward multi-domain repeat-derived ORFs. Per the stopping
condition, **the tier rule was not changed after seeing this**. It is an input for the next stage's
brief, not a redesign of this one. Redesign count: 0.

**Not done in stage 1:** RNA evidence was not run on the fixture (the interface and the no-op are
tested), BAM input, and any measurement of whether the tiers find real genes.

## Edit for publication (2026-10-09)

Before the repository was made public, two passages above that referred to earlier unpublished work
were rewritten: the paragraph motivating the evidence order now states the retrocopy expectation on
its own, and one sentence about that work was removed from *Ranking*. No rule, threshold, criterion or
result was changed.

## Correction (2026-10-09)

The T-cell receptor / immunoglobulin count in *What reading the top 100 showed* read 23; a recount of
`swissprot_best` in `runs/2026-10-09-rat-dev/all_candidates.tsv` (ranks 1–100) gives 24: 18 T-cell
receptor variable and 6 immunoglobulin variable entries. Corrected above; nothing else changed.
