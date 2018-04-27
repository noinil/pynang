#!/usr/bin/env python3

import glob
import numpy as np

def get_conpl_seq(seq):
    bp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    s = ''
    for b in seq:
        c = bp[b]
        s += c
    return s[::-1]


def main():
    # =======================
    # Prepare data structures
    # =======================
    # make list of all possible single base, base-step, trinucleotide base-steps...
    base_type = ['A', 'C', 'G', 'T']
    dinuc_list = []
    trinuc_list = []
    for b1 in base_type:
        for b2 in base_type:
            dinuc_list.append(b1+b2)
            for b3 in base_type:
                trinuc_list.append(b1+b2+b3)

    # create dictionaries for different bond types...
    bond_SB = {}
    for b in base_type:
        bond_SB[b] = []
    bond_SP = {}
    bond_PS = {}
    for bs in dinuc_list:
        bond_SP[bs] = []
        bond_PS[bs] = []
    # create dictionaries for different angle types...
    angle_PSP = {}
    angle_PSB = {}
    angle_SPS = {}
    angle_BSP = {}
    for bs in dinuc_list:
        angle_PSB[bs] = []
        angle_SPS[bs] = []
        angle_BSP[bs] = []
    # create dictionaries for different dihedral types...
    dihedral_PSPS = {}
    dihedral_SPSP = {}
    dihedral_SPSB = {}
    dihedral_BSPS = {}
    for bs in dinuc_list:
        dihedral_BSPS[bs] = []
        dihedral_SPSB[bs] = []
    for bt in trinuc_list:
        angle_PSP[bt] = []
        dihedral_PSPS[bt] = []
        dihedral_SPSP[bt] = []

    # ============================================================================

    def analyze_dna_structure(ninfo_name, seq):
        # create a list for every CG DNA particle
        full_seq = []
        for i, b in enumerate(seq):
            if i > 0:
                full_seq.append('P')
            full_seq.append('S')
            full_seq.append(b)

        # read in ninfo file
        with open(ninfo_name, 'r') as fin:
            for line in fin:
                words = line.split()
                if line.startswith('bond'):
                    imp, jmp, dist = int(words[6]) - 1, int(words[7]) - 1, float(words[8])
                    i_b, j_b = full_seq[imp], full_seq[jmp]
                    if i_b == 'P':
                        b1 = full_seq[imp - 1]
                        b2 = full_seq[jmp + 1]
                        bond_PS[b1+b2].append(dist)
                    elif i_b == 'S' and j_b == 'P':
                        b1 = full_seq[imp + 1]
                        b2 = full_seq[jmp + 2]
                        bond_SP[b1+b2].append(dist)
                    elif i_b == 'S' and j_b in base_type:
                        bond_SB[j_b].append(dist)
                if line.startswith('angl'):
                    imp, jmp, kmp, angl = int(words[7]) - 1, int(words[8]) - 1, int(words[9]) - 1, float(words[10])
                    i_b, j_b, k_b = full_seq[imp], full_seq[jmp], full_seq[kmp]
                    if i_b == 'S':
                        b1 = full_seq[imp + 1]
                        b2 = full_seq[kmp + 1]
                        angle_SPS[b1+b2].append(angl)
                    elif i_b == 'P':
                        if k_b == 'P':
                            b1 = full_seq[imp - 1]
                            b2 = full_seq[jmp + 1]
                            b3 = full_seq[kmp + 2]
                            angle_PSP[b1 + b2 + b3].append(angl)
                        else:
                            b1 = full_seq[imp - 1]
                            b2 = k_b
                            angle_PSB[b1+b2].append(angl)
                    else:
                        b1 = i_b
                        b2 = full_seq[kmp + 2]
                        angle_BSP[b1+b2].append(angl)
                if line.startswith('dihd'):
                    imp, jmp, kmp, lmp, dih = int(words[8]) - 1, int(words[9]) - 1, int(words[10]) - 1, int(words[11]) - 1, float(words[12])
                    i_b, j_b, k_b, l_b = full_seq[imp], full_seq[jmp], full_seq[kmp], full_seq[lmp]
                    if i_b == 'P':
                        b1 = full_seq[imp - 1]
                        b2 = full_seq[jmp + 1]
                        b3 = full_seq[lmp + 1]
                        trib = b1 + b2 + b3
                        dihedral_PSPS[trib].append(dih)
                    elif i_b == 'S':
                        if l_b == 'P':
                            b1 = full_seq[imp + 1]
                            b2 = full_seq[kmp + 1]
                            b3 = full_seq[lmp + 2]
                            trib = b1 + b2 + b3
                            dihedral_SPSP[trib].append(dih)
                        else:
                            b1 = full_seq[imp + 1]
                            b2 = l_b
                            dihedral_SPSB[b1 + b2].append(dih)
                    else:
                        b1 = i_b
                        b2 = full_seq[lmp + 1]
                        dihedral_BSPS[b1 + b2].append(dih)
    # ============================================================================

    # ========================
    # Read in data and analyze
    # ========================
    for seq_file in glob.glob('./*.fasta'):
        # print(seq_file)
        with open(seq_file, 'r') as fin:
            seq = fin.readline().strip()
            ceq = get_conpl_seq(seq)
        s1_ninfo_file = seq_file[:-6] + '_s1.ninfo'
        s2_ninfo_file = seq_file[:-6] + '_s2.ninfo'
        analyze_dna_structure(s1_ninfo_file, seq)
        analyze_dna_structure(s2_ninfo_file, ceq)

    fout_final = open("DNA_3SPN2C_structure_parameters.dat", "w")
    fout_str = "{0}    {1}     Average = {2:>10.3f}  +-  {3:>10.3f} \n"

    fout_final.write(' ======================================== \n')
    for k, v in bond_SB.items():
        varr = np.array(v)
        fout_final.write(fout_str.format('Bond SB:', k, varr.mean(), varr.std() ))
    fout_final.write(' ========== \n')
    for k, v in bond_SP.items():
        varr = np.array(v)
        fout_final.write(fout_str.format('Bond SP:', k, varr.mean(), varr.std() ))
    fout_final.write(' ========== \n')
    for k, v in bond_PS.items():
        varr = np.array(v)
        fout_final.write(fout_str.format('Bond PS:', k, varr.mean(), varr.std() ))
    fout_final.write(' ======================================== \n \n')

    fout_final.write(' ======================================== \n')
    for k, v in angle_PSP.items():
        varr = np.array(v)
        fout_final.write(fout_str.format('Angle PSP', k, varr.mean(), varr.std() ))
    fout_final.write(' ========== \n')
    for k, v in angle_SPS.items():
        varr = np.array(v)
        fout_final.write(fout_str.format('Angle SPS', k, varr.mean(), varr.std() ))
    fout_final.write(' ========== \n')
    for k, v in angle_PSB.items():
        varr = np.array(v)
        fout_final.write(fout_str.format('Angle PSB', k, varr.mean(), varr.std() ))
    fout_final.write(' ========== \n')
    for k, v in angle_BSP.items():
        varr = np.array(v)
        fout_final.write(fout_str.format('Angle BSP', k, varr.mean(), varr.std() ))
    fout_final.write(' ======================================== \n \n')

    fout_final.write(' ======================================== \n')
    for k, v in dihedral_PSPS.items():
        varr = np.array(v)
        fout_final.write(fout_str.format('Dihedral PSPS', k, varr.mean(), varr.std() ))
    fout_final.write(' ========== \n')
    for k, v in dihedral_SPSP.items():
        varr = np.array(v)
        fout_final.write(fout_str.format('Dihedral SPSP', k, varr.mean(), varr.std() ))
    fout_final.write(' ========== \n')
    for k, v in dihedral_SPSB.items():
        varr = np.array(v)
        fout_final.write(fout_str.format('Dihedral SPSB', k, varr.mean(), varr.std() ))
    fout_final.write(' ========== \n')
    for k, v in dihedral_BSPS.items():
        varr = np.array(v)
        fout_final.write(fout_str.format('Dihedral BSSP', k, varr.mean(), varr.std() ))
    fout_final.write(' ======================================== \n \n')

    fout_final.close()

if __name__ == '__main__':
    main()
