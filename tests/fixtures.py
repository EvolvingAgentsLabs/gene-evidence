"""Tiny synthetic inputs, written into a temporary directory. No network, no external tools: the raw
DIAMOND/hmmsearch tables are pre-seeded so `gene-evidence analyze` resumes from them."""
import json
import os

# locus: (seq, strand, [cds segments], status, confidence, protein)
LOCI = {
    "g1": ("CM1.1", "+", [(1000, 1300)], "complete", "0.5", "M" + "A" * 99),          # on a pseudogene
    "g2": ("CM1.1", "+", [(5000, 5300)], "complete", "0.4", "M" + "C" * 99),          # retrocopy
    "g3": ("CM1.1", "+", [(40000, 40500), (40800, 41000)], "complete", "0.3", "M" + "D" * 230),  # tier 1
    "g4": ("CM1.1", "-", [(50000, 50300)], "complete", "0.2", "M" + "E" * 99),        # Pfam only, lncRNA
    "g5": ("CM1.1", "+", [(60000, 60200)], "complete", "0.9", "A" + "F" * 66),        # no leading M
    "g6": ("CM1.1", "+", [(70000, 70150)], "complete", "0.1", "M" + "G" * 49),        # nothing
    "g7": ("CM1.1", "+", [(2000, 2100)], "complete", "0.8", "M" + "H" * 33),          # overlaps ref CDS
    "g8": ("CM2.1", "+", [(100, 400)], "complete", "0.7", "M" + "I" * 99),            # no identical twin
    "g9": ("CM1.1", "-", [(80000, 80100), (80300, 80400)], "incomplete", "NA", "M" + "K" * 66),
}

REF_GFF = """##gff-version 3
NC_1.1\tRefSeq\tpseudogene\t900\t1400\t.\t+\t.\tID=gene-PS1;Dbxref=GeneID:11;Name=Ps1;gene_biotype=pseudogene;pseudo=true
NC_1.1\tRefSeq\tgene\t1900\t2500\t.\t+\t.\tID=gene-PC1;Dbxref=GeneID:12;Name=Pc1;gene_biotype=protein_coding
NC_1.1\tRefSeq\tCDS\t2000\t2200\t.\t+\t0\tID=cds-XP_A;Parent=rna-1;Dbxref=GeneID:12;protein_id=XP_A.1
NC_1.1\tRefSeq\tgene\t20000\t30000\t.\t+\t.\tID=gene-PAR;Dbxref=GeneID:13;Name=Par;gene_biotype=protein_coding
NC_1.1\tRefSeq\tCDS\t20000\t20100\t.\t+\t0\tID=cds-XP_P;Parent=rna-2;Dbxref=GeneID:13;protein_id=XP_P.1
NC_1.1\tRefSeq\tCDS\t25000\t25100\t.\t+\t0\tID=cds-XP_P;Parent=rna-2;Dbxref=GeneID:13;protein_id=XP_P.1
NC_1.1\tRefSeq\tCDS\t29900\t30000\t.\t+\t0\tID=cds-XP_P;Parent=rna-2;Dbxref=GeneID:13;protein_id=XP_P.1
NC_1.1\tRefSeq\tCDS\t1000\t1300\t.\t+\t0\tID=cds-ps;Parent=gene-PS1;Dbxref=GeneID:11;pseudo=true
NC_1.1\tRefSeq\tgene\t49900\t50400\t.\t-\t.\tID=gene-LNC;Dbxref=GeneID:14;Name=Lnc1;gene_biotype=lncRNA
NC_1.1\tRefSeq\tgene\t79000\t81000\t.\t+\t.\tID=gene-PC2;Dbxref=GeneID:15;Name=Pc2;gene_biotype=protein_coding
NC_1.1\tRefSeq\tCDS\t79500\t80500\t.\t+\t0\tID=cds-XP_B;Parent=rna-3;Dbxref=GeneID:15;protein_id=XP_B.1
"""

REPORT = """# Assembly name: synthetic
# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn
1\tassembled-molecule\t1\tChromosome\tCM1.1\t=\tNC_1.1
2\tassembled-molecule\t2\tChromosome\tCM2.1\t<>\tNC_2.1
"""

SPROT_HITS = (
    "g1\tsp|P1|A_MOUSE\t99.0\t100\t100\t100\t1e-60\t200.0\t100.0\t"
    "sp|P1|A_MOUSE Protein A. Long OS=Mus musculus OX=10090 GN=a PE=1 SV=1\n"
    "g3\tsp|P4|D2_HUMAN\t50.0\t200\t231\t250\t1e-30\t120.0\t90.0\t"
    "sp|P4|D2_HUMAN Protein D2 OS=Homo sapiens OX=9606 PE=1 SV=1\n"
    "g3\tsp|P3|D_HUMAN\t60.5\t200\t231\t250\t1e-40\t150.5\t90.0\t"
    "sp|P3|D_HUMAN Protein D [x] OS=Homo sapiens OX=9606 GN=d PE=1 SV=1\n"
    "g6\tsp|P6|G_YEAST\t40.0\t10\t50\t300\t1e-6\t30.0\t20.0\t"
    "sp|P6|G_YEAST Protein G OS=Saccharomyces cerevisiae OX=4932 PE=1 SV=1\n"
)
PARALOG_HITS = (
    "g2\tXP_P.1\t98.0\t100\t100\t101\t1e-55\t190.0\t95.0\tXP_P.1 parent protein\n"
    "g9\tXP_P.1\t70.0\t50\t67\t101\t1e-10\t60.0\t75.0\tXP_P.1 parent protein\n"
)
PFAM_HITS = (
    "locus\tpfam_acc\tpfam_name\tscore\ti_evalue\tenv_from\tenv_to\n"
    "g3\tPF00001.25\t7tm_1\t120.3\t1e-35\t10\t200\n"
    "g4\tPF00002.30\t7tm_2\t45.0\t1e-12\t5\t90\n"
)


def write(d):
    os.makedirs(d, exist_ok=True)
    rows = []
    for lid, (seq, strand, cds, status, conf, _) in LOCI.items():
        s, e = min(a for a, _ in cds), max(b for _, b in cds)
        rows.append(f"{seq}\tGENERanno\tgene\t{s}\t{e}\t.\t{strand}\t.\tID={lid};status={status};confidence={conf}")
        rows.append(f"{seq}\tGENERanno\tmRNA\t{s}\t{e}\t.\t{strand}\t.\tID={lid}_t;Parent={lid}")
        for a, b in cds:
            rows.append(f"{seq}\tGENERanno\tCDS\t{a}\t{b}\t.\t{strand}\t0\tID={lid}_t.cds;Parent={lid}_t")
    p = {k: os.path.join(d, k) for k in ("carbon.gff3", "proteins.faa", "ref.gff", "report.txt",
                                          "sprot.fasta", "Pfam-A.hmm", "refprot.faa", "rna.gtf")}
    open(p["carbon.gff3"], "w").write("##gff-version 3\n" + "\n".join(rows) + "\n")
    open(p["proteins.faa"], "w").write("".join(f">{k} loc\n{v[5]}\n" for k, v in LOCI.items()))
    open(p["ref.gff"], "w").write(REF_GFF)
    open(p["report.txt"], "w").write(REPORT)
    open(p["sprot.fasta"], "w").write(">sp|P1|A_MOUSE x\nMAAA\n")
    open(p["Pfam-A.hmm"], "w").write("HMMER3/f placeholder\n")
    open(p["refprot.faa"], "w").write(">XP_P.1 parent\nMCCC\n")
    open(p["rna.gtf"], "w").write(
        'NC_1.1\tsrc\texon\t40000\t40500\t.\t+\t.\ttranscript_id "t1";\n'
        'NC_1.1\tsrc\texon\t40800\t41200\t.\t+\t.\ttranscript_id "t1";\n')
    return p


def seed_raw(out, reference=True):
    raw = os.path.join(out, "raw")
    os.makedirs(raw, exist_ok=True)
    open(os.path.join(raw, "sprot_hits.tsv"), "w").write(SPROT_HITS)
    if reference:
        open(os.path.join(raw, "paralog_hits.tsv"), "w").write(PARALOG_HITS)
    open(os.path.join(raw, "pfam_hits.tsv"), "w").write(PFAM_HITS)
    json.dump({"diamond": "0.0.test", "pyhmmer": "0.0.test"}, open(os.path.join(raw, "tools.json"), "w"))


def argv(p, out, reference=True, rna=False, top=100):
    a = ["analyze", "--gca", "GCA_TEST.1", "--carbon-gff", p["carbon.gff3"], "--carbon-proteins",
         p["proteins.faa"], "--swissprot", p["sprot.fasta"], "--swissprot-release", "test",
         "--pfam", p["Pfam-A.hmm"], "--pfam-release", "test", "--threads", "1", "--top", str(top),
         "--out", out]
    if reference:
        a += ["--ref-gff", p["ref.gff"], "--ref-proteins", p["refprot.faa"], "--ref-report",
              p["report.txt"], "--ref-name", "GCF_TEST.1-RS_TEST"]
    if rna:
        a += ["--rna-gtf", p["rna.gtf"]]
    return a
