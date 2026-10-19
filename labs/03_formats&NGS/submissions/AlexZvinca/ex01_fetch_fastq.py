"""
Exercise 03 — FASTQ Download (student-owned)

Objective:
- Choose a TP53-related accession (e.g., SRR..., ERR...) and DOWNLOAD a FASTQ file.
- Save it in data/work/<handle>/lab03/your_reads.fastq.gz

Minimum requirements:
- The script must accept an accession (e.g., from command line args).
- The script downloads at least one FASTQ (one file is enough for the exercise).
- The script prints to stdout the path of the downloaded file.

Recommended:
- Support .fastq or .fastq.gz.

NOTE:
- The chosen library does not matter (requests/urllib/etc.), but avoid heavy packages.
"""

import argparse
import csv
import urllib.parse
import urllib.request
from pathlib import Path
from io import StringIO

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("accession")
    parser.add_argument("--handle", required=True)
    return parser.parse_args()


def get_fastq_url(accession):
    base_url = "https://www.ebi.ac.uk/ena/portal/api/filereport"

    params = {
        "accession": accession,
        "result": "read_run",
        "fields": "fastq_ftp",
    }

    url = base_url + "?" + urllib.parse.urlencode(params)

    with urllib.request.urlopen(url) as response:
        text = response.read().decode("utf-8")

    rows = list(csv.DictReader(StringIO(text), delimiter="\t"))

    if not rows or not rows[0]["fastq_ftp"]:
        raise ValueError(f"No FASTQ found for {accession}")

    fastq_url = rows[0]["fastq_ftp"].split(";")[0]

    if not fastq_url.startswith(("http://", "https://", "ftp://")):
        fastq_url = "https://" + fastq_url

    return fastq_url


def get_output_path(handle):
    return Path(f"/workspaces/BDHB-lab-Zvinca/data/work/{handle}/lab03/your_reads.fastq.gz")


def download_fastq(url, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, output_path)

def main():
    # TODO: read the accession (e.g., sys.argv)
    # TODO: query the source (ENA/SRA) for the FASTQ link
    # TODO: download the file to the LOCATION YOU CHOSE
    # TODO: print("Downloaded:", <file_path>)
    
    args = parse_args()

    fastq_url = get_fastq_url(args.accession)
    output_path = get_output_path(args.handle)

    download_fastq(fastq_url, output_path)

    print("Downloaded:", output_path)


if __name__ == "__main__":
    main()
