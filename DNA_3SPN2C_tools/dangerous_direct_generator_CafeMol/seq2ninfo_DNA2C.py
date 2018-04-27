#!/usr/bin/env python3

def read_DNA_sequence(file_name):
    """Read in DNA sequence from .fasta file.
    """
    with open(file_name, 'r') as fin:
        for line in fin:
            if line.startswith('>'):
                continue
            seq = line.strip()
            for b in seq:
                if b not in ('A', 'C', 'G', 'T'):
                    break
            else:            
                return seq

def get_complementary_seq(seq):
    bp_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    s = ''
    for b in seq:
        s += bp_dict[b]
    return s[::-1]

def main(DNA_sequence_file_name, flag_phosphate_5, flag_infinite_DNA):
    ###########################################################################
    #                                Parameters                               #
    ###########################################################################
    # ===========
    # bond length
    # ===========
    bond_bdna = {
        "SB" : { "A" : 4.864, "C" : 4.300, "G" : 4.973, "T" : 4.379},
        "SP" : {
            "AA" : 3.688, "AC" : 3.018, "AG" : 3.836, "AT" : 3.287,
            "CA" : 4.386, "CC" : 3.538, "CG" : 4.676, "CT" : 3.999,
            "GA" : 3.736, "GC" : 3.256, "GG" : 3.633, "GT" : 3.285,
            "TA" : 4.191, "TC" : 3.707, "TG" : 4.391, "TT" : 3.868,},
        "PS" : {
            "XA" : 3.747, "XC" : 3.725, "XG" : 3.723, "XT" : 3.758,
            "AA" : 3.745, "AC" : 3.704, "AG" : 3.725, "AT" : 3.729,
            "CA" : 3.753, "CC" : 3.786, "CG" : 3.686, "CT" : 3.784,
            "GA" : 3.740, "GC" : 3.700, "GG" : 3.766, "GT" : 3.760,
            "TA" : 3.751, "TC" : 3.710, "TG" : 3.716, "TT" : 3.759,}}

    # ==========
    # bond angle
    # ==========
    angle_bdna = {
        "PSP" : {
            "AAA" : 120.685, "AAC" : 112.882, "AAG" : 113.827, "AAT" : 117.435,
            "ACA" : 119.061, "ACC" : 120.353, "ACG" : 113.240, "ACT" : 121.103,
            "AGA" : 122.182, "AGC" : 118.658, "AGG" : 120.489, "AGT" : 122.928,
            "ATA" : 117.235, "ATC" : 112.084, "ATG" : 111.714, "ATT" : 119.324,
            "CAA" : 122.866, "CAC" : 115.083, "CAG" : 116.036, "CAT" : 119.640,
            "CCA" : 120.442, "CCC" : 121.712, "CCG" : 114.602, "CCT" : 122.446,
            "CGA" : 124.721, "CGC" : 121.204, "CGG" : 122.937, "CGT" : 125.429,
            "CTA" : 119.317, "CTC" : 114.156, "CTG" : 113.756, "CTT" : 121.413,
            "GAA" : 120.809, "GAC" : 112.897, "GAG" : 113.816, "GAT" : 117.461,
            "GCA" : 119.550, "GCC" : 120.788, "GCG" : 113.687, "GCT" : 121.506,
            "GGA" : 121.512, "GGC" : 118.019, "GGG" : 119.634, "GGT" : 122.157,
            "GTA" : 117.087, "GTC" : 111.922, "GTG" : 111.501, "GTT" : 119.185,
            "TAA" : 122.361, "TAC" : 114.671, "TAG" : 115.653, "TAT" : 119.219,
            "TCA" : 121.235, "TCC" : 122.532, "TCG" : 115.417, "TCT" : 123.284,
            "TGA" : 123.936, "TGC" : 120.395, "TGG" : 122.319, "TGT" : 124.730,
            "TTA" : 119.004, "TTC" : 113.847, "TTG" : 113.465, "TTT" : 121.093},
        "SPS" : {
            "AA" : 94.805, "AC" : 94.462, "AG" : 95.308, "AT" : 95.232,
            "CA" : 95.110, "CC" : 98.906, "CG" : 92.244, "CT" : 97.476,
            "GA" : 94.973, "GC" : 92.666, "GG" : 97.929, "GT" : 97.640,
            "TA" : 94.886, "TC" : 93.066, "TG" : 93.999, "TT" : 95.122},
        "PSB" : {
            "XA" : 108.200, "XC" : 103.850, "XG" : 111.750, "XT" : 98.523,
            "AA" : 108.826, "AC" : 105.066, "AG" : 112.796, "AT" : 99.442,
            "CA" : 107.531, "CC" : 103.509, "CG" : 110.594, "CT" : 97.807,
            "GA" : 108.064, "GC" : 103.135, "GG" : 112.654, "GT" : 98.577,
            "TA" : 108.414, "TC" : 103.853, "TG" : 111.732, "TT" : 98.271},
        "BSP" : {
            "AA" : 113.855, "AC" : 114.226, "AG" : 112.201, "AT" : 111.931,
            "CA" : 113.822, "CC" : 112.056, "CG" : 116.081, "CT" : 111.008,
            "GA" : 114.665, "GC" : 118.269, "GG" : 110.102, "GT" : 111.146,
            "TA" : 113.984, "TC" : 115.457, "TG" : 113.397, "TT" : 113.606}}

    # ==============
    # dihedral angle
    # ==============
    dihedral_bdna = {
        "PSPS" : {
            "AAA" : -335.622, "AAC" : -332.885, "AAG" : -331.259, "AAT" : -336.185,
            "ACA" : -336.388, "ACC" : -335.577, "ACG" : -336.063, "ACT" : -337.660,
            "AGA" : -339.083, "AGC" : -339.751, "AGG" : -334.497, "AGT" : -339.668,
            "ATA" : -332.487, "ATC" : -331.938, "ATG" : -330.672, "ATT" : -335.597,
            "CAA" : -336.021, "CAC" : -332.981, "CAG" : -331.273, "CAT" : -336.309,
            "CCA" : -335.364, "CCC" : -334.499, "CCG" : -335.058, "CCT" : -336.547,
            "CGA" : -338.746, "CGC" : -339.509, "CGG" : -333.638, "CGT" : -339.033,
            "CTA" : -331.817, "CTC" : -331.269, "CTG" : -329.902, "CTT" : -334.955,
            "GAA" : -334.534, "GAC" : -331.854, "GAG" : -330.223, "GAT" : -335.116,
            "GCA" : -334.009, "GCC" : -333.155, "GCG" : -333.791, "GCT" : -335.211,
            "GGA" : -337.783, "GGC" : -338.478, "GGG" : -333.379, "GGT" : -338.439,
            "GTA" : -331.220, "GTC" : -330.726, "GTG" : -329.471, "GTT" : -334.288,
            "TAA" : -336.903, "TAC" : -333.864, "TAG" : -332.178, "TAT" : -337.225,
            "TCA" : -336.627, "TCC" : -335.754, "TCG" : -336.236, "TCT" : -337.799,
            "TGA" : -339.780, "TGC" : -340.478, "TGG" : -334.803, "TGT" : -340.164,
            "TTA" : -332.217, "TTC" : -331.655, "TTG" : -330.303, "TTT" : -335.342},
        "SPSP" : {
            "AAA" :   -0.215, "AAC" :   -6.669, "AAG" :   -8.623, "AAT" :   -6.140,
            "ACA" : -356.300, "ACC" : -357.745, "ACG" : -357.543, "ACT" : -358.626,
            "AGA" : -349.949, "AGC" : -348.414, "AGG" :   -0.166, "AGT" : -355.422,
            "ATA" : -359.491, "ATC" :   -0.267, "ATG" :   -2.823, "ATT" : -358.801,
            "CAA" : -359.648, "CAC" :   -6.270, "CAG" :   -8.270, "CAT" :   -5.727,
            "CCA" :   -1.694, "CCC" :   -3.186, "CCG" :   -2.836, "CCT" :   -4.099,
            "CGA" : -352.058, "CGC" : -350.458, "CGG" :   -2.541, "CGT" : -357.699,
            "CTA" :   -3.434, "CTC" :   -4.154, "CTG" :   -6.746, "CTT" :   -2.748,
            "GAA" :   -5.294, "GAC" :  -11.589, "GAG" :  -13.574, "GAT" :  -11.159,
            "GCA" :   -6.965, "GCC" :   -8.477, "GCG" :   -7.947, "GCT" :   -9.399,
            "GGA" : -354.234, "GGC" : -352.619, "GGG" :   -4.326, "GGT" : -359.682,
            "GTA" :   -5.833, "GTC" :   -6.486, "GTG" :   -9.031, "GTT" :   -5.212,
            "TAA" : -357.232, "TAC" :   -3.956, "TAG" :   -5.933, "TAT" :   -3.335,
            "TCA" : -357.663, "TCC" : -359.135, "TCG" : -358.953, "TCT" :   -0.035,
            "TGA" : -349.881, "TGC" : -348.369, "TGG" :   -0.320, "TGT" : -355.458,
            "TTA" :   -2.210, "TTC" :   -2.963, "TTG" :   -5.556, "TTT" :   -1.521},
        "SPSB" : {
            "AA" : -134.575, "AC" : -125.211, "AG" : -133.016, "AT" : -122.792,
            "CA" : -134.805, "CC" : -130.229, "CG" : -135.453, "CT" : -126.633,
            "GA" : -138.911, "GC" : -134.485, "GG" : -136.077, "GT" : -128.440,
            "TA" : -132.922, "TC" : -127.162, "TG" : -133.947, "TT" : -125.592},
        "BSPS" : {
            "AA" : -203.347, "AC" : -207.858, "AG" : -207.117, "AT" : -209.246,
            "CA" : -211.608, "CC" : -211.364, "CG" : -214.383, "CT" : -213.819,
            "GA" : -196.641, "GC" : -197.077, "GG" : -200.529, "GT" : -201.472,
            "TA" : -216.960, "TC" : -219.034, "TG" : -219.283, "TT" : -218.799}}

    # ===================
    # bond force constant
    # ===================
    # unit: kJ/mol/A^2
    bond_force_constant = 0.6   

    # ====================
    # angle force constant
    # ====================
    # unit: kJ/mol/rad^2
    angle_force_constant = {
        "BSP" : {
            "AA":460, "AT":370, "AC":442, "AG":358,
            "TA":120, "TT":460, "TC":383, "TG":206,
            "CA":206, "CT":358, "CC":278, "CG":278,
            "GA":383, "GT":442, "GC":336, "GG":278},
        "PSB" : {
            "XA":292, "XT":407, "XC":359, "XG":280,
            "AA":460, "TA":120, "CA":206, "GA":383,
            "AT":370, "TT":460, "CT":358, "GT":442,
            "AC":442, "TC":383, "CC":278, "GC":336,
            "AG":358, "TG":206, "CG":278, "GG":278},
        "PSP" : { "all" : 300},
        "SPS" : {
            "AA":355, "AT":147, "AC":464, "AG":368,
            "TA":230, "TT":355, "TC":442, "TG":273,
            "CA":273, "CT":368, "CC":165, "CG":478,
            "GA":442, "GT":464, "GC":228, "GG":165}}

    # =======================
    # dihedral force constant
    # =======================
    # unit: kJ/mol/rad^2

    # ==========
    # base types
    # ==========
    base_type = ['A', 'T', 'C', 'G']

    ###########################################################################


    ###########################################################################
    #                   Read in sequence and make ninfo file                  #
    ###########################################################################
    # read in DNA sequence:
    seq_DNA_a = read_DNA_sequence(DNA_sequence_file_name)
    seq_DNA_b = get_complementary_seq(seq_DNA_a)
    len_DNA = len(seq_DNA_a)
    print(" DNA lenth: ", len_DNA, " bp.")
    print(" ---------- ")
    print(" Strand 1:")
    print(" 5'-", seq_DNA_a, "-3'")
    print(" Strand 2:")
    print(" 3'-", seq_DNA_b[::-1], "-5'")
    print(" ---------- ")

    # Construct local topology information
    def CG_particle_names(seq, p5):
        dna_pname = []
        for i, b in enumerate(seq):
            if i > 0 or p5:
                dna_pname.append('P')
            dna_pname.append('S')
            dna_pname.append(b)
        return dna_pname
    strand_a_particles = CG_particle_names(seq_DNA_a, flag_phosphate_5)
    strand_b_particles = CG_particle_names(seq_DNA_b, flag_phosphate_5)

    # ===============================
    # functions generating parameters
    # ===============================
    def gen_bond_info(dna_particles):
        bonds = []
        for i, p in enumerate(dna_particles):
            if p == "P":
                # P-S
                if i == 0 and not flag_infinite_DNA:
                    b1 = "X"
                else:
                    b1 = dna_particles[i - 1]
                b2 = dna_particles[i + 2]
                dist = bond_bdna["PS"][b1 + b2]
                j = i + 1
                bonds.append((j, j + 1, dist, bond_force_constant))
            elif p == "S":
                # S-B
                b1 = dna_particles[i + 1]
                dist = bond_bdna["SB"][b1]
                j = i + 1
                bonds.append((j, j + 1, dist, bond_force_constant))
                # S-P
                if i < len(dna_particles) - 2:
                    b1 = dna_particles[i + 1]
                    b2 = dna_particles[i + 4]
                    dist = bond_bdna["SP"][b1 + b2]
                    j = i + 1
                    bonds.append((j, j + 2, dist, bond_force_constant))
                elif flag_infinite_DNA:
                    b1 = dna_particles[i + 1]
                    b2 = dna_particles[2]
                    dist = bond_bdna["SP"][b1 + b2]
                    j = i + 1
                    bonds.append((j, 3, dist, bond_force_constant))
        return bonds[:]
    def gen_angle_info(dna_particles):
        angles = []
        J_2_cal = 4.184
        for i, p in enumerate(dna_particles):
            if p == "P":
                # P-S-B
                if i == 0 and not flag_infinite_DNA:
                    b1 = "X"
                else:
                    b1 = dna_particles[i - 1]
                b2 = dna_particles[i + 2]
                a = angle_bdna["PSB"][b1 + b2]
                k = angle_force_constant["PSB"][b1 + b2] / J_2_cal
                j = i + 1
                angles.append((j, j + 1, j + 2, a, k))
                # P-S-P
                b1 = dna_particles[i - 1]
                if i < len(dna_particles) - 3:
                    b3 = dna_particles[i + 5]
                    a = angle_bdna["PSP"][b1 + b2 + b3]
                    k = angle_force_constant["PSP"]["all"] / J_2_cal
                    j = i + 1
                    angles.append((j, j + 1, j + 3, a, k))
                elif flag_infinite_DNA:
                    b3 = dna_particles[2]
                    a = angle_bdna["PSP"][b1 + b2 + b3]
                    k = angle_force_constant["PSP"]["all"] / J_2_cal
                    j = i + 1
                    angles.append((j, j + 1, 1, a, k))
            elif p == "S":
                # S-P-S
                b1 = dna_particles[i + 1]
                if i < len(dna_particles) - 2:
                    b2 = dna_particles[i + 4]
                    a = angle_bdna["SPS"][b1 + b2]
                    k = angle_force_constant["SPS"][b1 + b2] / J_2_cal
                    j = i + 1
                    angles.append((j, j + 2, j + 3, a, k))
                elif flag_infinite_DNA:
                    b2 = dna_particles[2]
                    a = angle_bdna["SPS"][b1 + b2]
                    k = angle_force_constant["SPS"][b1 + b2] / J_2_cal
                    j = i + 1
                    angles.append((j, 1, 2, a, k))
            else:               # p == "B"
                # B-S-P
                b1 = dna_particles[i]
                if i < len(dna_particles) - 1:
                    b2 = dna_particles[i + 3]
                    a = angle_bdna["BSP"][b1 + b2]
                    k = angle_force_constant["BSP"][b1 + b2] / J_2_cal
                    j = i + 1
                    angles.append((j, j - 1, j + 1, a, k))
                elif flag_infinite_DNA:
                    b2 = dna_particles[2]
                    a = angle_bdna["BSP"][b1 + b2]
                    k = angle_force_constant["BSP"][b1 + b2] / J_2_cal
                    j = i + 1
                    angles.append((j, j-1, 1, a, k))
        return angles[:]
    def gen_dihedral_info(dna_particles):
        dihedrals = []
        for i, p in enumerate(dna_particles):
            if p == "P":
                # P-S-P-S
                b1 = dna_particles[i - 1]
                b2 = dna_particles[i + 2]
                if i < len(dna_particles) - 3:
                    b3 = dna_particles[i + 5]
                    d = dihedral_bdna["PSPS"][b1 + b2 + b3]
                    k = "N2P2"
                    j = i + 1
                    dihedrals.append((j, j + 1, j + 3, j + 4, d, k))
                elif flag_infinite_DNA:
                    b3 = dna_particles[2]
                    d = dihedral_bdna["PSPS"][b1 + b2 + b3]
                    k = "N2P2"
                    j = i + 1
                    dihedrals.append((j, j + 1, 1, 2, d, k))
            elif p == "S":
                # S-P-S-B
                b1 = dna_particles[i + 1]
                if i < len(dna_particles) - 2:
                    b2 = dna_particles[i + 4]
                    d = dihedral_bdna["SPSB"][b1 + b2]
                    k = "N2P1"
                    j = i + 1
                    dihedrals.append((j, j + 2, j + 3, j + 4, d, k))
                elif flag_infinite_DNA:
                    b2 = dna_particles[2]
                    d = dihedral_bdna["SPSB"][b1 + b2]
                    k = "N2P1"
                    j = i + 1
                    dihedrals.append((j, 1, 2, 3, d, k))
                # S-P-S-P
                if i < len(dna_particles) - 5:
                    b2 = dna_particles[i + 4]
                    b3 = dna_particles[i + 7]
                    d = dihedral_bdna["SPSP"][b1 + b2 + b3]
                    k = "N2P2"
                    j = i + 1
                    dihedrals.append((j, j + 2, j + 3, j + 5, d, k))
                elif flag_infinite_DNA:
                    if i < len(dna_particles) - 2:
                        b2 = dna_particles[i + 4]
                        b3 = dna_particles[2]
                        d = dihedral_bdna["SPSP"][b1 + b2 + b3]
                        k = "N2P2"
                        j = i + 1
                        dihedrals.append((j, j + 2, j + 3, 1, d, k))
                    else:
                        b2 = dna_particles[2]
                        b3 = dna_particles[5]
                        d = dihedral_bdna["SPSP"][b1 + b2 + b3]
                        k = "N2P2"
                        j = i + 1
                        dihedrals.append((j, 1, 2, 4, d, k))
            else:               # p == "B"
                # B-S-P-S
                b1 = dna_particles[i]
                if i < len(dna_particles) - 1:
                    b2 = dna_particles[i + 3]
                    d = dihedral_bdna["BSPS"][b1 + b2]
                    k = "N2P1"
                    j = i + 1
                    dihedrals.append((j, j - 1, j + 1, j + 2, d, k))
                elif flag_infinite_DNA:
                    b2 = dna_particles[2]
                    d = dihedral_bdna["BSPS"][b1 + b2]
                    k = "N2P1"
                    j = i + 1
                    dihedrals.append((j, j - 1, 1, 2, d, k))
        return dihedrals

    # =====================
    # Output to ninfo files
    # =====================
    strand_a_bonds = gen_bond_info(strand_a_particles)
    strand_b_bonds = gen_bond_info(strand_b_particles)
    strand_a_angls = gen_angle_info(strand_a_particles)
    strand_b_angls = gen_angle_info(strand_b_particles)
    strand_a_dihes = gen_dihedral_info(strand_a_particles)
    strand_b_dihes = gen_dihedral_info(strand_b_particles)
    dna_bnd_list = [strand_a_bonds, strand_b_bonds]
    dna_ang_list = [strand_a_angls, strand_b_angls]
    dna_dih_list = [strand_a_dihes, strand_b_dihes]

    def write_ninfo():
        ninfo_bnd_head = "<<<< native bond length\n"
        ninfo_bnd_tail = ">>>>\n"
        ninfo_bnd_line = "bond  {serial:>5d}   {iunit:>4d}   {iunit:>4d} {bbond[0]:>6d} {bbond[1]:>6d} {bbond[0]:>6d} {bbond[1]:>6d}   {bbond[2]:>10.4f}       1.0000       1.0000       0.1434\n"

        ninfo_ang_head = "<<<< native bond angles \n"
        ninfo_ang_tail = ">>>> \n"
        ninfo_ang_line = "angl  {serial:>5d}   {iunit:>4d}   {iunit:>4d} {bangle[0]:>6d} {bangle[1]:>6d} {bangle[2]:>6d} {bangle[0]:>6d} {bangle[1]:>6d} {bangle[2]:>6d}   {bangle[3]:>10.4f}       1.0000       1.0000   {bangle[4]:>10.4f}\n"

        ninfo_dih_head = "<<<< native dihedral angles \n"
        ninfo_dih_tail = ">>>> \n"
        ninfo_dih_line = "dihd  {serial:>5d}   {iunit:>4d}   {iunit:>4d} {dihedral[0]:>6d} {dihedral[1]:>6d} {dihedral[2]:>6d} {dihedral[3]:>6d} {dihedral[0]:>6d} {dihedral[1]:>6d} {dihedral[2]:>6d} {dihedral[3]:>6d}   {dihedral[4]:>10.4f}       1.0000       1.0000       1.6730       0.3000 {dihedral[5]}\n"


        for j in range(2):
            ninfo_name = "strand{0}.ninfo".format(j + 1)
            ninfo_file = open(ninfo_name, 'w')
            strnd_bnd_list = dna_bnd_list[j]
            strnd_ang_list = dna_ang_list[j]
            strnd_dih_list = dna_dih_list[j]
            # write bond information
            ninfo_file.write(ninfo_bnd_head)
            for i, b in enumerate(strnd_bnd_list):
                ninfo_file.write(ninfo_bnd_line.format(serial=i + 1, iunit=j + 1, bbond=b))
            ninfo_file.write(ninfo_bnd_tail)
            ninfo_file.write("\n")
            # write angle information
            ninfo_file.write(ninfo_ang_head)
            for i, a in enumerate(strnd_ang_list):
                ninfo_file.write(ninfo_ang_line.format(serial=i + 1, iunit=j + 1, bangle=a))
            ninfo_file.write(ninfo_ang_tail)
            ninfo_file.write("\n")
            # write angle information
            ninfo_file.write(ninfo_dih_head)
            for i, d in enumerate(strnd_dih_list):
                ninfo_file.write(ninfo_dih_line.format(serial=i + 1, iunit=j + 1, dihedral=d))
            ninfo_file.write(ninfo_dih_tail)

            ninfo_file.close()

    write_ninfo()


if __name__ == '__main__':
    import argparse
    def parse_arguments():
        parser = argparse.ArgumentParser(description='Generate DNA curvature parameter file from DNA sequence.')
        parser.add_argument('fasta', type=str, help="DNA sequence file name.")
        parser.add_argument('-p', '--headPHOS', help="Include 5' phosphate group or not.", action="store_true")
        parser.add_argument('-i', '--infinite', help="Infinite DNA.", action="store_true")
        return parser.parse_args()
    args = parse_arguments()
    if args.infinite and not args.headPHOS:
        print(" 5'-phosphate must be included while setting infinite DNA. ")
        exit()
    main(args.fasta, args.headPHOS, args.infinite)
