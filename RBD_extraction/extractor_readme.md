# Protocol to extract RBD from GISAID sequences

## 1. Install repository
https://github.com/jbloomlab/SARS-CoV-2-RBD_MAP_clinical_Abs

## 2. Create conda environment

conda env create -f environment.yml
conda activate SARS-CoV-2-RBD_MAP

## 3. Add folders in repo: 
- add rbd_extractorbis.py in main
- add  reference sequence in ATGC with name wildtype_sequence.fasta in /data (step is not useful?)
- add reference sequence in residues with name reference_RBD.fasta in results/GISAID_mutations
- add GISAID data with name spikeprot0414.fasta in data/spikeprot0414

## 4. Run script
python3 rbd_extractorbis.py 