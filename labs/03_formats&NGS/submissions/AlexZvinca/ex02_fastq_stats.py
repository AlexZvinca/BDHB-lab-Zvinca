"""
Exercise 04 — FASTQ QC on your own data

TODO:
- Read your FASTQ file from data/work/<handle>/lab03/:
    your_reads.fastq  or  your_reads.fastq.gz
- Compute statistics:
    * total number of reads
    * mean read length
    * proportion of 'N' bases
    * mean Phred score
- Save the report in:
    labs/03_formats&NGS/submissions/<handle>/qc_report_<handle>.txt
"""

import os
import gzip
from pathlib import Path
from Bio import SeqIO

# TODO: replace <handle> with your GitHub username
handle = "AlexZvinca"

in_fastq_plain = Path(f"data/work/{handle}/lab03/your_reads.fastq")
in_fastq_gz = Path(f"data/work/{handle}/lab03/your_reads.fastq.gz")
out_report = Path(f"labs/03_formats&NGS/submissions/{handle}/qc_report_{handle}.txt")
out_report.parent.mkdir(parents=True, exist_ok=True)

# Select the existing file
if in_fastq_plain.exists():
    reader = SeqIO.parse(str(in_fastq_plain), "fastq")
elif in_fastq_gz.exists():
    # Biopython reads from file-like objects; we use gzip.open(..., "rt")
    reader = SeqIO.parse(gzip.open(in_fastq_gz, "rt"), "fastq")
else:
    raise FileNotFoundError(
        f"Could not find either {in_fastq_plain} or {in_fastq_gz}. "
        f"Run ex03_fetch_fastq.py first or copy your own FASTQ."
    )

num_reads = 0
total_length = 0
total_n = 0
total_phred = 0
total_bases = 0

# TODO: complete the aggregation logic
for record in reader:
    # HINT:
    # seq_str = str(record.seq)
    # phred = record.letter_annotations["phred_quality"]

    seq_str = str(record.seq).upper()
    phred = record.letter_annotations["phred_quality"]

    num_reads += 1
    total_length += len(seq_str)
    total_n += seq_str.count("N")

    total_phred += sum(phred)
    total_bases += len(phred)

# TODO: compute final values (watch out for division by zero)
len_mean = total_length / num_reads if num_reads > 0 else 0.0
n_rate = total_n / total_length if total_length > 0 else 0.0
phred_mean = total_phred / total_bases if total_bases > 0 else 0.0

with open(out_report, "w", encoding="utf-8") as out:
    out.write(f"Reads: {num_reads}\n")
    out.write(f"Mean length: {len_mean:.2f}\n")
    out.write(f"N rate: {n_rate:.4f}\n")
    out.write(f"Mean Phred: {phred_mean:.2f}\n")

print(f"[OK] QC report -> {out_report.resolve()}")
