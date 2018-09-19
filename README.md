# coarse-grained-DDMP

### Execution:

python coarse-grained-ddmp.py structure_1 chain_1 structure_2 chain_2 helix_indices.txt offset

### Example execution:
python coarse-grained-ddmp.py G45R A G223W A DraNramp_helices.txt 38

where inputs are:

G45R.pdb

G223W.pdb

DraNramp_helices.txt (Helix indices are inputted in new-line-separated text file, with both the beginning and end indices of each helix specified, i.e. there should be an even number of indices in the text file. Indices can be partially determined by DSSP, https://swift.cmbi.umcn.nl/gv/dssp/, but it may be best to curate the cutoffs of helices manually)

Offset is an integer referring to the shift of the first aligned residue, e.g. if the first aligned residue is residue 39, the offset will be 38.

### Outputs:
Aligned FASTA of two sequences (G45R_A-G223W_A_aligned.fasta)

CSV of DDM (G45R_AvsG223W_A.csv)

PNG of DDM (G45R_AvsG223W_A.png)

CSV of coarse-grained DDM (G45R_AvsG223W_A-helices.csv)

PNG of coarse-grained DDM (G45R_AvsG223W_A-helices.png)
