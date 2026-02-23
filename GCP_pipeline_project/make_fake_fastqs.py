import random

DNA = "ACGT"
QUAL = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI"  # standard-ish ASCII

def rand_seq(n=80):
    return "".join(random.choice(DNA) for _ in range(n))

def rand_qual(n=80):
    return "".join(random.choice(QUAL) for _ in range(n))

def main(out_path="sample.fastq", n_reads=200, read_len=80):
    with open(out_path, "w") as f:
        for i in range(n_reads):
            f.write(f"@FAKE_READ_{i}\n")
            f.write(rand_seq(read_len) + "\n")
            f.write("+\n")
            f.write(rand_qual(read_len) + "\n")
    print(f"Wrote {n_reads} reads to {out_path}")

if __name__ == "__main__":
    main()