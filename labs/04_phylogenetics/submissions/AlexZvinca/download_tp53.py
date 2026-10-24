from pathlib import Path
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

Entrez.email = "alexandru.zvinca@student.upt.ro"

accessions = {
    "Human": "NM_000546.6",
    "Mouse": "NM_011640.4",
    "Rat": "NM_030989.4",
    "Cow": "NM_174201.2",
}

output = Path("data/work/AlexZvinca/lab04/your_sequences.fasta")


def download_cds(species, accession):
    with Entrez.efetch(
        db="nuccore",
        id=accession,
        rettype="gb",
        retmode="text",
    ) as handle:
        record = SeqIO.read(handle, "genbank")

    cds_features = [
        feature for feature in record.features
        if feature.type == "CDS"
    ]

    if len(cds_features) != 1:
        raise ValueError(f"Expected one CDS for {accession}")

    cds = cds_features[0].extract(record.seq)

    return SeqRecord(cds, id=species, description=accession)


def main():
    output.parent.mkdir(parents=True, exist_ok=True)

    sequences = []

    for species, accession in accessions.items():
        sequence = download_cds(species, accession)
        sequences.append(sequence)

        print(f"{species}: {len(sequence.seq)} bases")

    SeqIO.write(sequences, output, "fasta")
    print("Saved:", output)


if __name__ == "__main__":
    main()