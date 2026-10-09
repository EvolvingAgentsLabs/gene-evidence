# Images to add

The documentation references three raster images that do not exist yet. Each reference is wrapped in
an HTML comment, so nothing renders broken until the file is added. The two diagrams (the pipeline in
`README.md` and the evidence graph in `docs/EVIDENCE.md`) are Mermaid and need no image.

To add an image: put the file at the exact path below, then make the edit listed under *Un-comment*.

---

## 1. `hero.png`

- **Path:** `docs/img/hero.png`
- **Size / aspect:** 1600 × 640 px (5:2). Keep it under 400 KB.
- **Format:** PNG (or WebP renamed in the reference).
- **What it shows:** an illustrative banner, not data. A horizontal genome track with one predicted gene
  locus highlighted; thin citation lines run from the locus to a few small evidence cards (locus context,
  retrocopy check, homology, domains), each card tagged with a short node id in brackets. It should read
  as "every claim is tied to its evidence", not as a result.
- **Image-model prompt (ready to paste):**

  > Clean flat technical illustration, wide 5:2 banner, white background, muted blue-grey palette with
  > one teal accent. A thin horizontal genome track runs across the lower third, drawn as a line with
  > small exon blocks. One gene locus on the track is highlighted in teal. From that locus, five thin
  > dotted lines rise to five small rounded rectangular cards arranged in an arc above it. Each card
  > has a tiny icon and one short label: "locus context", "retrocopy?", "homology", "domains",
  > "ORF". Under each label a small monospace tag in square brackets, like "[ctx:…]". No other text,
  > no logos, no people, no DNA double helix, no glow, no gradients, no 3D. Generous whitespace,
  > crisp vector lines, suitable for a GitHub README header.

- **Alt text:** A predicted gene locus on a genome track, with evidence cards linked to it by citation
  lines
- **Referenced in:** `README.md`, directly under the title.
- **Un-comment:** in `README.md`, replace

  ```
  <!-- IMAGE PLACEHOLDER: see docs/img/IMAGES.md (hero)
  ![A predicted gene locus on a genome track, with evidence cards linked to it by citation lines](docs/img/hero.png)
  -->
  ```

  with

  ```
  ![A predicted gene locus on a genome track, with evidence cards linked to it by citation lines](docs/img/hero.png)
  ```

---

## 2. `report-example.png`

- **Path:** `docs/img/report-example.png`
- **Size / aspect:** 1400 px wide, height as rendered (about 1400 × 1100). Retina capture scaled to
  1400 px is fine. Keep it under 500 KB.
- **Format:** PNG.
- **What it shows:** a screenshot of one real candidate section of the committed report, as GitHub
  renders it: the rank-1 candidate of the rat development run, with its *Evidence FOR*, *Evidence
  AGAINST*, *Context* and *Suggested validation* blocks, every line ending in `[node-id]` citations.
- **Steps:**
  1. Source: `runs/2026-10-09-rat-dev/run1/report.md`, **lines 117–142** (from
     `### 1. GCA_036323735.1_gene00003511: tier 1 (FOR-strong)` through the line that starts
     `- Tier 1 (FOR-strong) by the fixed tier rule`).
  2. Render it on GitHub: open
     `https://github.com/EvolvingAgentsLabs/gene-evidence/blob/main/runs/2026-10-09-rat-dev/run1/report.md#1-gca_0363237351_gene00003511-tier-1-for-strong`
     in a light-theme browser window 1400 px wide. (Alternative, offline: copy lines 117–142 into a
     scratch `.md` file and render it in any GitHub-flavoured Markdown preview.)
  3. Capture from the `### 1.` heading to the end of *Suggested validation*. Do not crop the citations
     off the right edge; if the long Pfam line wraps, that is fine.
  4. Do not edit the text. The image must match the committed report byte for byte.
- **Alt text:** One candidate section of report.md, with evidence FOR, evidence AGAINST, context and a
  suggested experiment, each line ending in graph-node citations
- **Referenced in:** `README.md`, end of *Inputs and outputs*.
- **Un-comment:** in `README.md`, replace

  ```
  <!-- IMAGE PLACEHOLDER: see docs/img/IMAGES.md (report-example)
  ![One candidate section of report.md, with evidence FOR, evidence AGAINST, context and a suggested experiment, each line ending in graph-node citations](docs/img/report-example.png)
  -->
  ```

  with

  ```
  ![One candidate section of report.md, with evidence FOR, evidence AGAINST, context and a suggested experiment, each line ending in graph-node citations](docs/img/report-example.png)
  ```

---

## 3. `rat-locus-context.png`

- **Path:** `docs/img/rat-locus-context.png`
- **Size / aspect:** 1200 × 600 px (2:1), 150 dpi. Keep it under 200 KB.
- **Format:** PNG.
- **What it shows:** a data chart of the rat development run. Horizontal stacked bars, one per locus
  context class (pseudogene, nothing, protein-coding gene, non-coding gene), length = number of
  candidates, segments = tier. It makes visible why locus context comes first: the pseudogene bar
  (578 of 947) is the largest and is entirely tier 5. Expected counts (from `all_candidates.tsv`):

  | context | tier 1 | tier 2 | tier 3 | tier 5 | total |
  |---|--:|--:|--:|--:|--:|
  | pseudogene | 0 | 0 | 0 | 578 | 578 |
  | nothing | 172 | 60 | 32 | 21 | 285 |
  | protein-coding gene | 26 | 8 | 7 | 4 | 45 |
  | non-coding gene | 25 | 7 | 5 | 2 | 39 |

  (Tier 4 is empty on this run.)
- **Steps:** from the repository root, with matplotlib installed (`pip install matplotlib`), run:

  ```python
  import csv, collections
  import matplotlib.pyplot as plt

  rows = list(csv.DictReader(open("runs/2026-10-09-rat-dev/all_candidates.tsv"), delimiter="\t"))
  ctx = ["pseudogene", "nothing", "protein-coding gene", "non-coding gene"]
  tiers = [("1", "1 FOR-strong", "#2a7f62"), ("2", "2 FOR", "#7fbf9f"),
           ("3", "3 no evidence", "#bdbdbd"), ("5", "5 AGAINST-strong", "#b2453a")]
  n = collections.Counter((r["locus_context"], r["tier"]) for r in rows)
  assert sum(n.values()) == 947

  fig, ax = plt.subplots(figsize=(8, 4), dpi=150)
  left = [0] * len(ctx)
  for t, label, color in tiers:
      w = [n[(c, t)] for c in ctx]
      ax.barh(ctx, w, left=left, color=color, label=label)
      left = [a + b for a, b in zip(left, w)]
  for i, total in enumerate(left):
      ax.text(total + 5, i, str(total), va="center", fontsize=9)
  ax.invert_yaxis()
  ax.set_xlabel("candidates (rat dev run, n = 947)")
  ax.set_title("Locus context of candidates, by tier")
  ax.spines[["top", "right"]].set_visible(False)
  ax.legend(frameon=False, fontsize=8, loc="lower right")
  fig.tight_layout()
  fig.savefig("docs/img/rat-locus-context.png")
  ```

  Check the bar totals against the table above before committing.
- **Alt text:** Rat development run: the 947 candidates by locus context, coloured by tier; the
  pseudogene bar is the largest and is entirely tier 5
- **Referenced in:** `docs/EVIDENCE.md`, end of *The confound to keep in mind*. The path there is
  relative to `docs/` (`img/rat-locus-context.png`).
- **Un-comment:** in `docs/EVIDENCE.md`, replace

  ```
  <!-- IMAGE PLACEHOLDER: see docs/img/IMAGES.md (rat-locus-context)
  ![Rat development run: the 947 candidates by locus context, coloured by tier; the pseudogene bar is the largest and is entirely tier 5](img/rat-locus-context.png)
  -->
  ```

  with

  ```
  ![Rat development run: the 947 candidates by locus context, coloured by tier; the pseudogene bar is the largest and is entirely tier 5](img/rat-locus-context.png)
  ```
