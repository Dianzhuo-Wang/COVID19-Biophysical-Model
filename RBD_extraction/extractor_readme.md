# Protocol to extract RBD from GISAID sequences

## 1. Install repository
https://github.com/jbloomlab/SARS-CoV-2-RBD_MAP_clinical_Abs

## 2. Create conda environment

conda env create -f environment.yml
conda activate SARS-CoV-2-RBD_MAP

## 3. Add folders in repo: 
- add rbd_extractor.py in main
- add wildtype sequence in residues with name wildtype_sequence.fasta in data
- add GISAID data with name spikeprot0414.fasta in data/spikeprot0414

## 4. Run script
python3 rbd_extractor.py 