#!/usr/bin/env python
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import Bio.PDB
import Bio.PDB.DSSP as DSSP
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from matplotlib.colors import Normalize

if len(sys.argv) <= 6 or len(sys.argv[2]) != 1 or len(sys.argv[4]) != 1:
    print 'usage: python coarse-grained-ddmp.py 1st_structure_name chain 2nd_structure_name chain helices.txt offset'
    exit()

# input structure names and chains
a_0 = sys.argv[1]
a_chain = sys.argv[2]
b_0 = sys.argv[3]
b_chain = sys.argv[4]
helices_file = sys.argv[5]
offset = int(sys.argv[6])

# input pdb structures
a = a_0+"_"+a_chain
b = b_0+"_"+b_chain

letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
           'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
           'TYR':'Y','VAL':'V'}

# write PDB with just the chain we want
def chain(name, chain):
    with open(name+".pdb", "r") as input_file, open(name+"_"+chain+".pdb", "w") as output_file:
        for line in input_file:
            toks = line.split()
            if len(toks)<5: continue
            if toks[3] not in letters: continue
            if toks[4] == chain:
                output_file.write(line)
        #output_file.write("\n")
    input_file.close()
    output_file.close()

# convert PDB to FASTA
def pdb2fasta(name):
    with open(name+".pdb", "r") as input_file, open(name+".fasta", "w") as output_file:
        output_file.write(">"+name+"\n")
        prev = "-1"
        for line in input_file:
            toks = line.split()
            if len(toks)<1: continue
            if toks[0] != 'ATOM': continue
            if toks[3] not in letters: continue
            if toks[5] != prev:
                output_file.write('%c' % letters[toks[3]])
            prev = toks[5]
        output_file.write("\n")
    input_file.close()
    output_file.close()


# FASTA concatenation and MUSCLE alignment
def muscleAlign(a, b):
    filenames = [a+".fasta", b+".fasta"]
    with open(a+"-"+b+".fasta", 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
    
    cline = MuscleCommandline(input=a+'-'+b+'.fasta', out=a+'-'+b+'_aligned.fasta')
    print cline
    os.system(str(cline))
    #os.system('muscle3.8.31_i86win32.exe -in '+a+'-'+b+'.fasta -out '+a+'-'+b+'_aligned.fasta')

# build DDM and save PNG
def ddmp(a, b, shift):
    model_a = Bio.PDB.PDBParser().get_structure(a, a+'.pdb')[0]
    model_b = Bio.PDB.PDBParser().get_structure(b, b+'.pdb')[0]
    a_chains = list(model_a.get_chains())
    b_chains = list(model_b.get_chains())
    print a_chains
    print b_chains
    # designate which chains to compare
    a_residues = list(a_chains[0])
    b_residues = list(b_chains[0])
    
    # align the residues in the PDB
    alignment = AlignIO.read(a+"-"+b+"_aligned.fasta", "fasta")
    #print alignment[0].id, "\n", alignment[0].seq
    #print alignment[1].id, "\n", alignment[1].seq
    align_len =  alignment.get_alignment_length()
    a_res = -1
    b_res = -1
    # tuple list of matching residues
    res_aligned = []
    for k in range(0, align_len):
        if alignment[0].seq[k] != "-":
            a_res += 1
        if alignment[1].seq[k] != "-":
            b_res += 1
        if alignment[0].seq[k] != "-" and alignment[1].seq[k] != "-":
             res_aligned.append((a_res, b_res))
        else:
            res_aligned.append((-1,-1))
    #print res_aligned
    # number of matching residues
    print "Length of sequence alignment: {}".format(align_len)
    n = sum(1 for i in res_aligned if i != (-1,-1))
    print "Number of aligned residues: {}".format(n)
    
    # DDM & DDM_aligned (without spaces), write DDM .txt file, and plot matrix
    ddm = [[np.nan for x in range(0,align_len)] for y in range(0,align_len)]
    ddm_aligned = [[np.nan for x in range(0,n)] for y in range(0,n)]
    text_file = open(a+"vs"+b+".txt", "w")

    # iterate through all pairwise distances for both a, b, then subtract
    for i in range(0, n):
        for j in range(0, n):
            dist_a = abs(a_residues[res_aligned[i][0]]['CA'] - a_residues[res_aligned[j][0]]['CA'])
            dist_b = abs(b_residues[res_aligned[i][1]]['CA'] - b_residues[res_aligned[j][1]]['CA'])
            ddm_aligned[i][j] = dist_a - dist_b

    # do the same for DDM, except in larger matrix leave NaN for non-aligned
    for i in range(0, align_len):
        for j in range(0, align_len):
            if res_aligned[i] != (-1,-1) and res_aligned[j] != (-1,-1):
                dist_a = abs(a_residues[res_aligned[i][0]]['CA'] - a_residues[res_aligned[j][0]]['CA'])
                dist_b = abs(b_residues[res_aligned[i][1]]['CA'] - b_residues[res_aligned[j][1]]['CA'])
                dd = dist_a - dist_b
                ddm[i][j] = dd
                text_file.write("{}\t".format(dd))
            else:
                text_file.write("{}\t".format("N"))
        text_file.write("\n")
    
    text_file.close()
    
    #print ddm

    # plot and save DDMP
    ddm_shifted = np.empty((align_len+shift+1, align_len+shift+1))
    ddm_shifted[:] = np.nan
    for i in range(0, align_len):
        for j in range(0, align_len):
            ddm_shifted[i+shift+1][j+shift+1] = ddm[i][j]
    np.savetxt(a+'vs'+b+".csv", ddm_shifted, delimiter=",")
    print a+"vs"+b+".csv has been saved."
    return ddm_shifted

def plotddm(ddm_shifted, helix_indices):
    first_aligned = np.where(~np.isnan(ddm_shifted))[0][0]
    print "First aligned residue: " + str(first_aligned)
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    ax1.set_aspect('equal')
    ax1.set_xlim([first_aligned, len(ddm_shifted)])
    ax1.set_ylim([len(ddm_shifted), first_aligned])
    ax1.set_xticks(range(first_aligned, len(ddm_shifted), 50))
    ax1.set_yticks(helix_indices)
    ddm_masked = np.ma.masked_where(np.isnan(ddm_shifted), ddm_shifted)
    palette = plt.cm.seismic
    palette.set_bad(color="lightgray")
    #print ddm
    plt.imshow(ddm_shifted, vmin=-10, vmax=10, cmap=palette)

    #plt.imshow(ddm, norm=norm, cmap=palette)
    plt.colorbar()
    plt.savefig(a+"vs"+b+".png")
    
    print a+"vs"+b+".png has been saved."

def rmsdddm(ddm_shifted, helix_indices):
    rmsd_ddm = [[0 for x in range(0, len(helix_indices)/2)] for y in range(0, len(helix_indices)/2)]
    #print len(helix_indices)
    for i in range(0, len(helix_indices)/2):
        for j in range(0, len(helix_indices)/2):
            helix1a = helix_indices[2*i]
            helix1b = helix_indices[2*i+1]
            helix2a = helix_indices[2*j]
            helix2b = helix_indices[2*j+1]
            rmsd = 0
            numnotnan = 0
            for k in range(helix1a, helix1b+1):
                for l in range(helix2a, helix2b+1):
                    if not (np.isnan(ddm_shifted[k][l])):
                        numnotnan += 1
                        rmsd += (ddm_shifted[k][l])**2
            rmsd_ddm[i][j] = np.sqrt(rmsd/(numnotnan+1))
    np.savetxt(a+'vs'+b+"-helices.csv", rmsd_ddm, delimiter=",")
    print a+"vs"+b+"-helices.csv has been saved."
    #print rmsd_ddm
    fig = plt.figure()
    plt.matshow(rmsd_ddm, cmap=plt.cm.afmhot, vmin=0, vmax=5)
    plt.colorbar()
    plt.savefig(a+"vs"+b+"-helices.png")
    print a+"vs"+b+"-helices.png has been saved."

if __name__ == '__main__':
    chain(a_0, a_chain)
    chain(b_0, b_chain)
    # FASTA concatenation and MUSCLE alignment
    pdb2fasta(a)
    pdb2fasta(b)
    muscleAlign(a, b)
    f = open(helices_file, 'r')
    helices = []
    vals = f.read().split('\n')
    vals = vals[:-1]
    for val in vals:
        helices.append(int(val))
    f.close()
    helix_starts = helices[::2]
    DDM = ddmp(a,b,offset)
    plotddm(DDM, helix_starts)
    rmsdddm(DDM, helices)
