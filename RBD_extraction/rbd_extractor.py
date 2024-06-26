import io
import lzma
import os
import re
import subprocess
import tarfile

from Bio.Data.IUPACData import protein_letters
import Bio.SeqIO

from IPython.display import display, HTML

import matplotlib.pyplot as plt

import pandas as pd

from plotnine import *
from tqdm.notebook import tqdm

import yaml

with open("config.yaml") as f:
    config = yaml.safe_load(f)

os.makedirs(config["gisaid_mutations_dir"], exist_ok=True)

print(f"Reading GISAID spikes in {config['gisaid_spikes']}")


spikes = []

with open(
    "data/spikeprot0414/spikeprot0414.fasta", "r", encoding="utf-8", errors="ignore"
) as f:
    print("opened")
    fasta_records = Bio.SeqIO.parse(f, "fasta")

    for record in fasta_records:
        spikes.append(record)


print(f"Read {len(spikes)} spike sequences.")

spikes_df = pd.DataFrame({"seqrecord": spikes}).assign(
    description=lambda x: x["seqrecord"].map(lambda rec: rec.description),
    country=lambda x: x["description"].str.split("|").str[-1],
    host=lambda x: x["description"].str.split("|").str[6].str.strip(),
    date=lambda x: x["description"].str.split("|").str[2].str.strip(),
    length=lambda x: x["seqrecord"].map(len),
    n_ambiguous=lambda x: x["seqrecord"].map(
        lambda rec: rec.seq.count("X") + rec.seq.count("x")
    ),
)


spikes_df = spikes_df.query('host == "Human"')
print("Here are the number of human sequences:", len(spikes_df))


min_length, max_length = 1260, 1276
print(f"\nOnly keeping spikes with lengths between {min_length} and {max_length}")
spikes_df = spikes_df.assign(
    valid_length=lambda x: (min_length <= x["length"]) & (x["length"] <= max_length)
)

print("Here are number of sequences with valid lengths:")
print(
    spikes_df.groupby("valid_length").aggregate(
        n_sequences=pd.NamedAgg("seqrecord", "count")
    )
)


spikes_df = spikes_df.query("valid_length")

max_ambiguous = 100
print(f"Filtering sequences with > {max_ambiguous} ambiguous residues")
spikes_df = spikes_df.assign(
    excess_ambiguous=lambda x: x["n_ambiguous"] > max_ambiguous
)
print(
    spikes_df.groupby("excess_ambiguous").aggregate(
        n_sequences=pd.NamedAgg("seqrecord", "count")
    )
)


# ## Align the RBD region of the spikes
# We now align the RBD regions of the spikes.
# We do this **before** we filter sequences with too many ambiguous residues so that we can do that filtering just on the RBD region.
#
# We align with `mafft` using the `--addfragments` and `--keeplength` options (see [here](https://mafft.cbrc.jp/alignment/software/closelyrelatedviralgenomes.html) and [here](https://mafft.cbrc.jp/alignment/software/addsequences.html)) to align relative to a reference that is just the RBD; these options also clip the sequences relative to the reference.
# Note that these options make sense if the following conditions are met:
#   1. Sequences are all very similar.
#   2. We are not worried about insertions.
# For now, both of these appear to be true, but this choice should be kept in mind if there is a lot more divergence or insertions.
#
# We align relative to the reference that is the wildtype sequence used for the experiments:
print(f"Reading reference sequence in data/wildtype_sequence.fasta")
refseq = Bio.SeqIO.read("data/wildtype_sequence.fasta", "fasta")
refprotfile = "data/reference_RBD.fasta"
print(f"Writing the reference RBD to {refprotfile}")
_ = Bio.SeqIO.write(refseq, refprotfile, "fasta")

spikes_file = os.path.join(
    config["gisaid_mutations_dir"], "human_full-length_spikes.fasta"
)
print(f"Writing the spikes to {spikes_file}")
_ = Bio.SeqIO.write(spikes_df["seqrecord"].tolist(), spikes_file, "fasta")

# Now make the alignment.
# Note that we use multiple threads to speed things up, and also align the spikes in chunks.
# The reason that we have to the chunkwise alignment is that some unclear `mafft` error was arising if we tried to align them all at once:

chunksize = 50000

aligned_rbds = []

for i in range(0, len(spikes_df), chunksize):
    spikes_file = os.path.join(
        config["gisaid_mutations_dir"],
        f"human_full-length_spikes_{i + 1}-to-{i + chunksize}.fasta",
    )
    print(f"Writing spikes {i + 1} to {i + chunksize} to {spikes_file}")
    _ = Bio.SeqIO.write(
        spikes_df["seqrecord"].tolist()[i : i + chunksize], spikes_file, "fasta"
    )
    print("Now aligning these sequences...")
    cmds = [
        "mafft",
        "--auto",
        "--thread",
        str(config["max_cpus"]),
        "--keeplength",
        "--addfragments",
        spikes_file,
        refprotfile,
    ]
    res = subprocess.run(cmds, capture_output=True)
    if res.returncode:
        raise RuntimeError(f"Error in alignment:\n{res.stderr}")
    else:
        print("Alignment complete.\n")
        try:
            with io.StringIO(res.stdout.decode("utf-8")) as f:
                iseqs = list(Bio.SeqIO.parse(f, "fasta"))
                # remove reference sequence, which should be first in file
                assert (
                    iseqs[0].seq == refseq.seq
                    and iseqs[0].description == refseq.description
                )
                iseqs = iseqs[1:]
                assert len(iseqs) == min(chunksize, len(spikes_df) - i)
                aligned_rbds += iseqs

        except UnicodeDecodeError:
            # Handle the UnicodeDecodeError here
            # You can print an error message, log the error, or take any other appropriate action
            pass

print("number of sequences with unicode error:", len(spikes_df) - len(aligned_rbds))


# Now put all of the aligned RBDs into a data frame to filter and parse:

rbd_df = pd.DataFrame({"seqrecord": aligned_rbds}).assign(
    description=lambda x: x["seqrecord"].map(lambda rec: rec.description),
    country=lambda x: x["description"].str.split("|").str[-1],
    host=lambda x: x["description"].str.split("|").str[6].str.strip(),
    length=lambda x: x["seqrecord"].map(len),
    date=lambda x: x["description"].str.split("|").str[2],
    n_ambiguous=lambda x: x["seqrecord"].map(
        lambda rec: rec.seq.count("X") + rec.seq.count("x")
    ),
    n_gaps=lambda x: x["seqrecord"].map(lambda rec: rec.seq.count("-")),
    all_valid_aas=lambda x: x["seqrecord"].map(
        lambda rec: re.fullmatch(f"[{protein_letters}]+", str(rec.seq)) is not None
    ),
)

assert all(rbd_df["length"] == len(refseq))

# Based on above plots, we will retain just RBDs with no ambiguous amino acids and no gaps:

rbd_df = rbd_df.query("n_ambiguous == 0").query("n_gaps == 0")
assert rbd_df["all_valid_aas"].all()
print(f"Retained {len(rbd_df)} RBDs with no gap or ambiguous.")

# Now get and plot the number of amino-acid mutations per RBD relative to the reference sequence, plotting on both a linear and log scale.


refseq_str = str(refseq.seq)
rbd_df = rbd_df.assign(
    seq=lambda x: x["seqrecord"].map(lambda rec: str(rec.seq)),
    n_mutations=lambda x: x["seq"].map(
        lambda s: sum(x != y for x, y in zip(s, refseq_str))
    ),
)


print(f"Overall, there are {len(rbd_df)} aligned RBDs that passed filters.")


rbd_df = rbd_df[["seq", "date", "n_mutations"]]

rbd_df.to_csv("rbd_df.csv", index=False)


# Create a function to write a single sequence to the FASTA file
def write_fasta_record(file, header, sequence):
    file.write(f">{header}\n")
    file.write(f"{sequence}\n")


# Open the FASTA file for writing
with open("rbd_df.fasta", "w") as fasta_file:
    # Iterate over each row in the DataFrame
    for index, row in rbd_df.iterrows():
        sequence = row["seq"]
        date = row["date"]
        n_mutations = row["n_mutations"]

        # Create a header for the FASTA record (you can modify this as needed)
        header = f"Sequence_{index} | Date_{date} | NMutations_{n_mutations}"

        # Write the sequence to the FASTA file
        write_fasta_record(fasta_file, header, sequence)
