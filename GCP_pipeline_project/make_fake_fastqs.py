#!/usr/bin/env python3
import random
from pathlib import Path

DNA = "ACGT"
QUAL = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI"

def rand_dna(n: int) -> str:
    return "".join(random.choice(DNA) for _ in range(n))

def rand_qual(n: int) -> str:
    return "".join(random.choice(QUAL) for _ in range(n))

def write_fasta(path: Path, contigs: dict[str, str]) -> None:
    with path.open("w") as f:
        for name, seq in contigs.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

def write_gtf(path: Path, genes: list[dict]) -> None:
    with path.open("w") as f:
        for g in genes:
            attrs = f'gene_id "{g["gene_id"]}"; gene_name "{g["gene_name"]}";'
            f.write(
                f'{g["contig"]}\ttoy\texon\t{g["start"]}\t{g["end"]}\t.\t{g["strand"]}\t.\t{attrs}\n'
            )

def pick_fragment(seq: str, start0: int, end0: int, frag_len: int) -> str:
    """Pick a fragment fully inside [start0,end0) on a contig sequence."""
    if end0 - start0 < frag_len:
        raise ValueError(f"Interval too small for frag_len={frag_len}")
    pos0 = random.randint(start0, end0 - frag_len)
    return seq[pos0:pos0 + frag_len]

def gene_interval_0based(g: dict) -> tuple[int, int]:
    # GTF: 1-based inclusive start/end
    start0 = g["start"] - 1
    end0 = g["end"]          # exclusive for python slicing
    return start0, end0

def make_read_from_gene_chunks(
    contigs: dict[str, str],
    genes: list[dict],
    read_len: int = 120,
    n_gene_chunks: int = 2,
    gene_chunk_len: int = 35,
    min_spacer: int = 5,
    max_spacer: int = 25,
    p_noise_read: float = 0.05,
) -> tuple[str, str]:
    """
    Builds a read by concatenating:
      gene chunk(s) + spacer chunk(s) from the reference background.
    All pieces come from the reference -> BWA can map.
    Returns (read_seq, tag_for_header).
    """
    if random.random() < p_noise_read:
        # Fully random read (may not map). Keep this low.
        return rand_dna(read_len), "NOISE"

    parts: list[str] = []
    used_genes: list[str] = []

    # pick a gene to anchor the read (for better mapping)
    g0 = random.choice(genes)
    used_genes.append(g0["gene_id"])
    contig = g0["contig"]
    contig_seq = contigs[contig]

    g0_start0, g0_end0 = gene_interval_0based(g0)

    # First gene chunk (anchor)
    parts.append(pick_fragment(contig_seq, g0_start0, g0_end0, gene_chunk_len))

    # Alternate spacer + additional gene chunks
    for _ in range(n_gene_chunks - 1):
        spacer_len = random.randint(min_spacer, max_spacer)

        # spacer from anywhere on same contig (background)
        parts.append(pick_fragment(contig_seq, 0, len(contig_seq), spacer_len))

        g = random.choice([x for x in genes if x["contig"] == contig] or genes)
        used_genes.append(g["gene_id"])
        gs0, ge0 = gene_interval_0based(g)
        parts.append(pick_fragment(contig_seq, gs0, ge0, gene_chunk_len))

    seq = "".join(parts)

    # If too short, pad with reference background; if too long, trim.
    if len(seq) < read_len:
        pad_len = read_len - len(seq)
        seq += pick_fragment(contig_seq, 0, len(contig_seq), pad_len)
    else:
        seq = seq[:read_len]

    tag = "+".join(used_genes)
    return seq, tag

def write_fastq_samples(
    out_prefix: str,
    contigs: dict[str, str],
    genes: list[dict],
    n_samples: int = 3,
    n_reads: int = 200,
    read_len: int = 120,
):
    for s in range(1, n_samples + 1):
        out = Path(f"{out_prefix}{s}.fastq")
        with out.open("w") as f:
            for i in range(n_reads):
                seq, tag = make_read_from_gene_chunks(
                    contigs=contigs,
                    genes=genes,
                    read_len=read_len,
                    n_gene_chunks=2,        # tweak: 1,2,3...
                    gene_chunk_len=40,      # tweak
                    min_spacer=5,
                    max_spacer=25,
                    p_noise_read=0.02,      # keep low if you want mapping
                )
                f.write(f"@READ_{s}_{i}_{tag}\n")
                f.write(seq + "\n+\n")
                f.write(rand_qual(len(seq)) + "\n")
        print(f"Wrote {out} ({n_reads} reads, read_len={read_len})")

def main(seed: int = 1337):
    random.seed(seed)

    # Reference contigs
    contigs = {
        "chr1": rand_dna(2500),
        "chr2": rand_dna(2500),
    }

    # GTF genes: exon blocks on those contigs
    genes = [
        {"gene_id": "geneA", "gene_name": "geneA", "contig": "chr1", "start": 200,  "end": 900,  "strand": "+"},
        {"gene_id": "geneB", "gene_name": "geneB", "contig": "chr1", "start": 1200, "end": 1900, "strand": "+"},
        {"gene_id": "geneC", "gene_name": "geneC", "contig": "chr2", "start": 400,  "end": 1400, "strand": "-"},
    ]

    write_fasta(Path("reference.fasta"), contigs)
    write_gtf(Path("annotation.gtf"), genes)

    write_fastq_samples(
        out_prefix="sample_align",
        contigs=contigs,
        genes=genes,
        n_samples=3,
        n_reads=200,
        read_len=120,
    )

    print("Done: reference.fasta, annotation.gtf, sample_align{1..3}.fastq")

if __name__ == "__main__":
    main()
