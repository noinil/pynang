#!/usr/bin/env python3

import numpy as np
import MDAnalysis 

def main(PDB_name, flag_head_phos, flag_psf_output):
    # ==========================
    # Variables can be modified:
    # ==========================
    chainID_DNA_a = 'A'
    chainID_DNA_b = 'B'

    # ============================================
    # Constants for CG particle masses and charges
    # ============================================
    std_base_mass = {'A': 134.1, 'G': 150.1, 'C': 110.1, 'T': 125.1}

    # ============================================
    # Read in structural informatino from PDB file
    # ============================================
    u = MDAnalysis.Universe(PDB_name)

    selstr_DNA_a = "nucleic and segid {0}".format(chainID_DNA_a)
    selstr_DNA_b = "nucleic and segid {0}".format(chainID_DNA_b)

    selstr_P = "(resid {0} and (name P or name OP* or name O5')) or (resid {1} and name O3')"
    selstr_S = "resid {0} and (name C5' or name C4' or name C3' or name C2' or name C1' or name O4' or name O2')"
    selstr_B = "resid {0} and not (name C5' or name C4' or name C3' or name C2' or name C1' or name O4' or name O2' or name O3' or name O5' or name OP* or name P or name H*)"
    selstr_P0 = "resid {0} and (name P)"

    sel_dna_a = u.select_atoms(selstr_DNA_a)
    sel_dna_b = u.select_atoms(selstr_DNA_b)

    def cg_dna_top(dna_atom_group):
        """Translate atomistic DNA structure into coarse-grained topology.
        Input: dna_atom_group - Atomistic DNA atom/residue group.
        Output: Coarse-grained DNA particles.
        """
        if flag_head_phos == 0:
            cg_dna_particle_num = len(dna_atom_group.residues) * 3 - 1
        else:
            cg_dna_particle_num = len(dna_atom_group.residues) * 3 
        cg_dna_coors = np.empty([cg_dna_particle_num, 3])
        cg_dna_particle_resname = []
        cg_dna_particle_name = []
        cg_dna_particle_charge = []
        cg_dna_particle_mass = []
        resid_list = list(dna_atom_group.residues.resids)
        for i, j in enumerate(resid_list):
            tmp_resname = dna_atom_group.residues[i].resname
            # Phosphate
            if i > 0 or flag_head_phos == 1:
                if i == 0:
                    res_P = dna_atom_group.select_atoms(selstr_P0.format(j))
                else:
                    res_P = dna_atom_group.select_atoms(selstr_P.format(j, j-1))
                cg_dna_coors[i * 3 - 1 + flag_head_phos] = res_P.center_of_mass()
                cg_dna_particle_name.append('DP')
                cg_dna_particle_charge.append(-0.6)
                cg_dna_particle_mass.append(94.97)
                cg_dna_particle_resname.append(tmp_resname)
            # Sugar
            res_S = dna_atom_group.select_atoms(selstr_S.format(j))
            cg_dna_coors[i * 3 + flag_head_phos] = res_S.center_of_mass()
            cg_dna_particle_name.append('DS')
            cg_dna_particle_charge.append(0.0)
            cg_dna_particle_mass.append(83.11)
            cg_dna_particle_resname.append(tmp_resname)
            # Base
            res_B = dna_atom_group.select_atoms(selstr_B.format(j))
            cg_dna_coors[i * 3 + 1 + flag_head_phos] = res_B.center_of_mass()
            cg_dna_particle_name.append('DB')
            cg_dna_particle_charge.append(0.0)
            cg_dna_particle_mass.append(std_base_mass[tmp_resname[-1]])
            cg_dna_particle_resname.append(tmp_resname)
        return (cg_dna_coors, cg_dna_particle_name, cg_dna_particle_charge, cg_dna_particle_mass, cg_dna_particle_resname)

    cg_dna_a_coors, cg_dna_a_p_name, cg_dna_a_p_charge, cg_dna_a_p_mass, cg_dna_a_p_resname = cg_dna_top(sel_dna_a)
    cg_dna_b_coors, cg_dna_b_p_name, cg_dna_b_p_charge, cg_dna_b_p_mass, cg_dna_b_p_resname = cg_dna_top(sel_dna_b)

    # Assign CG particle ID and residue ID
    cg_dna_a_p_ID = []
    cg_dna_b_p_ID = []

    cg_dna_a_r_ID = []
    cg_dna_b_r_ID = []

    cg_dna_a_p_num = cg_dna_a_coors.shape[0]
    cg_dna_b_p_num = cg_dna_b_coors.shape[0]

    cg_dna_a_r_num = len(sel_dna_a.residues)
    cg_dna_b_r_num = len(sel_dna_b.residues)

    for i in range(cg_dna_a_p_num):
        cg_dna_a_p_ID.append(i + 1)
        if flag_head_phos == 1:
            cg_dna_a_r_ID.append(i // 3 + 1)
        else:
            cg_dna_a_r_ID.append((i + 1) // 3 + 1)
    for i in range(cg_dna_b_p_num):
        cg_dna_b_p_ID.append(i + 1 + cg_dna_a_p_num)
        if flag_head_phos == 1:
            cg_dna_b_r_ID.append(i // 3 + 1 + cg_dna_a_r_num)
        else:
            cg_dna_b_r_ID.append((i + 1) // 3 + 1 + cg_dna_a_r_num)

    def output_psf():
        """Output psf file for protein-DNA complex.
        """
        PSF_HEAD_STR = "PSF CMAP \n\n"
        PSF_TITLE_STR0 = "      3 !NTITLE \n"
        PSF_TITLE_STR1 = "REMARKS PSF file created with PWMcos tools. \n"
        PSF_TITLE_STR2 = "REMARKS DNA: {0:>5d} bp. \n"
        PSF_TITLE_STR5 = "REMARKS ======================================== \n"
        PSF_TITLE_STR6 = "       \n"
        PSF_TITLE_STR = PSF_TITLE_STR0 + PSF_TITLE_STR1 + PSF_TITLE_STR2 + PSF_TITLE_STR5 + PSF_TITLE_STR6
        PSF_ATOM_TITLE = " {atom_num:>6d} !NATOM \n"
        PSF_ATOM_LINE = " {atom_ser:>6d} {seg_id:>3} {res_ser:>5d} {res_name:>3} {atom_name:>3} {atom_type:>5}  {charge:>10.6f}  {mass:>10.6f}          0 \n"

        psf_file_name = "dna_cg.psf"
        psf_file = open(psf_file_name, 'w')
        psf_file.write(PSF_HEAD_STR)
        psf_file.write(PSF_TITLE_STR.format(cg_dna_a_r_num))
        psf_file.write(PSF_ATOM_TITLE.format(atom_num = cg_dna_a_p_num + cg_dna_b_p_num))

        for i in range(cg_dna_a_p_num):
            psf_file.write(PSF_ATOM_LINE.format(atom_ser = i + 1,
                                                seg_id = 'a',
                                                res_ser = cg_dna_a_r_ID[i],
                                                res_name = cg_dna_a_p_resname[i],
                                                atom_name = cg_dna_a_p_name[i],
                                                atom_type = cg_dna_a_p_name[i][-1],
                                                charge = cg_dna_a_p_charge[i],
                                                mass = cg_dna_a_p_mass[i]))
        for i in range(cg_dna_b_p_num):
            psf_file.write(PSF_ATOM_LINE.format(atom_ser = i + 1 + cg_dna_a_p_num,
                                                seg_id = 'b',
                                                res_ser = cg_dna_b_r_ID[i], 
                                                res_name = cg_dna_b_p_resname[i],
                                                atom_name = cg_dna_b_p_name[i],
                                                atom_type = cg_dna_b_p_name[i][-1],
                                                charge = cg_dna_b_p_charge[i],
                                                mass = cg_dna_b_p_mass[i]))
        psf_file.close()

    if flag_psf_output:
        output_psf()


    ###########################################################################
    #                       Determine .ninfo parameters                       #
    ###########################################################################

    def compute_bond(coor1, coor2):
        vec = coor1 - coor2
        return np.linalg.norm(vec)
    def compute_angle(coor1, coor2, coor3):
        vec1 = coor1 - coor2
        vec2 = coor3 - coor2
        n_v1 = np.linalg.norm(vec1)
        n_v2 = np.linalg.norm(vec2)
        return np.arccos(np.clip(np.dot(vec1, vec2) / n_v1 / n_v2, -1.0, 1.0)) / np.pi * 180
    def compute_dihedral(coor1, coor2, coor3, coor4):
        v_12 = coor2 - coor1
        v_23 = coor3 - coor2
        v_34 = coor4 - coor3
        n123 = np.cross(v_12, v_23)
        n234 = np.cross(v_23, v_34)
        norm_n123 = np.linalg.norm(n123)
        norm_n234 = np.linalg.norm(n234)
        dih = np.arccos(np.clip(np.dot(n123, n234) / norm_n123 / norm_n234, -1.0, 1.0)) / np.pi * 180.0 
        # determine sign of dih
        n1234 = np.cross(n123, n234)
        zajiao = np.dot(n1234, v_23)
        if zajiao < 0:
            dih = - dih
        return dih - 180.0

    def get_angle_param(angle_type, base_step):
        # get angle parameters based on angle type and base step (sequence)
        # Angle force constants (kJ/mol/rad^2)

        angle_params = {
            # Base-Sugar-Phosphate
            "BSP" : {
                "AA":460, "AT":370, "AC":442, "AG":358,
                "TA":120, "TT":460, "TC":383, "TG":206,
                "CA":206, "CT":358, "CC":278, "CG":278,
                "GA":383, "GT":442, "GC":336, "GG":278
            },
            # Phosphate-Sugar-Base
            "PSB" : {
                "AA":460, "TA":120, "CA":206, "GA":383,
                "AT":370, "TT":460, "CT":358, "GT":442,
                "AC":442, "TC":383, "CC":278, "GC":336,
                "AG":358, "TG":206, "CG":278, "GG":278
            },
            # Phosphate-Sugar-Phosphate
            "PSP" : { "all" : 300
            },
            # Sugar-Phosphate-Sugar
            "SPS" : {
                "AA":355, "AT":147, "AC":464, "AG":368,
                "TA":230, "TT":355, "TC":442, "TG":273,
                "CA":273, "CT":368, "CC":165, "CG":478,
                "GA":442, "GT":464, "GC":228, "GG":165
            }
        }
        return angle_params[angle_type][base_step] / 4.184

    dna_bnd_list = [[], []]
    dna_ang_list = [[], []]
    dna_dih_list = [[], []]
    for j in range(2):          # j =                         = 0: chain A, j = = 1: chain B
        cg_dna_p_num     = cg_dna_a_p_num     if j == 0 else cg_dna_b_p_num
        cg_dna_coors     = cg_dna_a_coors     if j == 0 else cg_dna_b_coors
        cg_dna_p_ID      = cg_dna_a_p_ID      if j == 0 else cg_dna_b_p_ID
        cg_dna_p_name    = cg_dna_a_p_name    if j == 0 else cg_dna_b_p_name
        cg_dna_p_resname = cg_dna_a_p_resname if j == 0 else cg_dna_b_p_resname
        for i_dna in range(cg_dna_p_num):
            if cg_dna_p_name[i_dna] == "DS":
                # bond S--B
                coor_s = cg_dna_coors[i_dna]
                coor_b = cg_dna_coors[i_dna + 1]
                r_sb = compute_bond(coor_s, coor_b)
                dna_bnd_list[j].append((cg_dna_p_ID[i_dna], cg_dna_p_ID[i_dna + 1], i_dna + 1, i_dna + 1 + 1, r_sb))
                if i_dna + 4 < cg_dna_p_num:
                    # bond S--P+1
                    coor_p3 = cg_dna_coors[i_dna + 2]
                    r_sp3 = compute_bond(coor_s, coor_p3)
                    dna_bnd_list[j].append((cg_dna_p_ID[i_dna], cg_dna_p_ID[i_dna + 2], i_dna + 1, i_dna + 2 + 1, r_sp3))
                    # Angle S--P+1--S+1
                    resname5 = cg_dna_p_resname[i_dna][-1]
                    resname3 = cg_dna_p_resname[i_dna + 3][-1]
                    coor_s3 = cg_dna_coors[i_dna + 3]
                    ang_sp3s3 = compute_angle(coor_s, coor_p3, coor_s3)
                    k = get_angle_param("SPS", resname5 + resname3)
                    dna_ang_list[j].append((cg_dna_p_ID[i_dna], cg_dna_p_ID[i_dna + 2], cg_dna_p_ID[i_dna + 3], i_dna + 1, i_dna + 2 + 1, i_dna + 3 + 1, ang_sp3s3, k))
                    # Dihedral S--P+1--S+1--B+1
                    coor_b3 = cg_dna_coors[i_dna + 4]
                    dih_sp3s3b3 = compute_dihedral(coor_s, coor_p3, coor_s3, coor_b3)
                    dna_dih_list[j].append((cg_dna_p_ID[i_dna], cg_dna_p_ID[i_dna + 2], cg_dna_p_ID[i_dna + 3], cg_dna_p_ID[i_dna + 4], i_dna + 1, i_dna + 2 + 1, i_dna + 3 + 1, i_dna + 4 + 1, dih_sp3s3b3, "N2P1"))
                    # Dihedral S--P+1--S+1--P+2
                    if i_dna + 6 < cg_dna_p_num:
                        coor_p33 = cg_dna_coors[i_dna + 5]
                        dih_sp3s3p33 = compute_dihedral(coor_s, coor_p3, coor_s3, coor_p33)
                        dna_dih_list[j].append((cg_dna_p_ID[i_dna], cg_dna_p_ID[i_dna + 2], cg_dna_p_ID[i_dna + 3], cg_dna_p_ID[i_dna + 5], i_dna + 1, i_dna + 2 + 1, i_dna + 3 + 1, i_dna + 5 + 1, dih_sp3s3p33, "N2P2"))
            elif cg_dna_p_name[i_dna] == "DP":
                # bond P--S
                coor_p = cg_dna_coors[i_dna]
                coor_s = cg_dna_coors[i_dna + 1]
                r_ps = compute_bond(coor_p, coor_s)
                dna_bnd_list[j].append((cg_dna_p_ID[i_dna], cg_dna_p_ID[i_dna + 1], i_dna + 1, i_dna + 1 + 1, r_ps))
                # angle P--S--B
                resname5 = cg_dna_p_resname[i_dna - 1][-1]
                resname3 = cg_dna_p_resname[i_dna + 2][-1]
                coor_b = cg_dna_coors[i_dna + 2]
                ang_psb = compute_angle(coor_p, coor_s, coor_b)
                k = get_angle_param("PSB", resname5 + resname3)
                dna_ang_list[j].append((cg_dna_p_ID[i_dna], cg_dna_p_ID[i_dna + 1], cg_dna_p_ID[i_dna + 2], i_dna + 1, i_dna + 1 + 1, i_dna + 2 + 1, ang_psb, k))
                if i_dna + 4 < cg_dna_p_num:
                    # angle P--S--P+1
                    coor_p3 = cg_dna_coors[i_dna + 3]
                    ang_psp3 = compute_angle(coor_p, coor_s, coor_p3)
                    k = get_angle_param("PSP", "all")
                    dna_ang_list[j].append((cg_dna_p_ID[i_dna], cg_dna_p_ID[i_dna + 1], cg_dna_p_ID[i_dna + 3], i_dna + 1, i_dna + 1 + 1, i_dna + 3 + 1, ang_psp3, k))
                    # Dihedral P--S--P+1--S+1
                    coor_s3 = cg_dna_coors[i_dna + 4]
                    dih_psp3s3 = compute_dihedral(coor_p, coor_s, coor_p3, coor_s3)
                    dna_dih_list[j].append((cg_dna_p_ID[i_dna], cg_dna_p_ID[i_dna + 1], cg_dna_p_ID[i_dna + 3], cg_dna_p_ID[i_dna + 4], i_dna + 1, i_dna + 1 + 1, i_dna + 3 + 1, i_dna + 4 + 1, dih_psp3s3, "N2P2"))
            elif cg_dna_p_name[i_dna] == "DB":
                if i_dna + 3 < cg_dna_p_num:
                    # angle B--S--P+1
                    resname5 = cg_dna_p_resname[i_dna][-1]
                    resname3 = cg_dna_p_resname[i_dna + 1][-1]
                    coor_b = cg_dna_coors[i_dna]
                    coor_s = cg_dna_coors[i_dna - 1]
                    coor_p3 = cg_dna_coors[i_dna + 1]
                    ang_bsp3 = compute_angle(coor_b, coor_s, coor_p3)
                    k = get_angle_param("BSP", resname5 + resname3)
                    dna_ang_list[j].append((cg_dna_p_ID[i_dna], cg_dna_p_ID[i_dna - 1], cg_dna_p_ID[i_dna + 1], i_dna + 1, i_dna - 1 + 1, i_dna + 1 + 1, ang_bsp3, k))
                    # Dihedral B--S--P+1--S+1
                    coor_s3 = cg_dna_coors[i_dna + 2]
                    dih_bsp3s3 = compute_dihedral(coor_b, coor_s, coor_p3, coor_s3)
                    dna_dih_list[j].append((cg_dna_p_ID[i_dna], cg_dna_p_ID[i_dna - 1], cg_dna_p_ID[i_dna + 1], cg_dna_p_ID[i_dna + 2], i_dna + 1, i_dna - 1 + 1, i_dna + 1 + 1, i_dna + 2 + 1, dih_bsp3s3, "N2P1"))
            else:
                print("Error! Wrong DNA type...")
                

    # =========================
    # Output CafeMol ninfo file
    # =========================
    def write_ninfo():
        ninfo_bnd_head = "<<<< native bond length\n"
        ninfo_bnd_tail = ">>>>\n"
        ninfo_bnd_line = "bond  {serial:>5d}   {iunit:>4d}   {iunit:>4d} {bbond[0]:>6d} {bbond[1]:>6d} {bbond[2]:>6d} {bbond[3]:>6d}   {bbond[4]:>10.4f}       1.0000       1.0000       0.1434\n"

        ninfo_ang_head = "<<<< native bond angles \n"
        ninfo_ang_tail = ">>>> \n"
        ninfo_ang_line = "angl  {serial:>5d}   {iunit:>4d}   {iunit:>4d} {bangle[0]:>6d} {bangle[1]:>6d} {bangle[2]:>6d} {bangle[3]:>6d} {bangle[4]:>6d} {bangle[5]:>6d}   {bangle[6]:>10.4f}       1.0000       1.0000   {bangle[7]:>10.4f}\n"

        ninfo_dih_head = "<<<< native dihedral angles \n"
        ninfo_dih_tail = ">>>> \n"
        ninfo_dih_line = "dihd  {serial:>5d}   {iunit:>4d}   {iunit:>4d} {dihedral[0]:>6d} {dihedral[1]:>6d} {dihedral[2]:>6d} {dihedral[3]:>6d} {dihedral[4]:>6d} {dihedral[5]:>6d} {dihedral[6]:>6d} {dihedral[7]:>6d}   {dihedral[8]:>10.4f}       1.0000       1.0000       1.6730       0.3000 {dihedral[9]}\n"


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
        parser = argparse.ArgumentParser(description='Generate CafeMol ninfo file from DNA PDB.')
        parser.add_argument('pdb', type=str, help="PDB file name.")
        parser.add_argument('-t', '--psf', help="Output psf file.", action="store_true")
        parser.add_argument('-p', '--headPHOS', help="Specify whether the 5' phosphate group will be used in CafeMol or not.", action="store_true")
        return parser.parse_args()
    args = parse_arguments()
    flag_head_phos = 1 if args.headPHOS else 0
    main(args.pdb, flag_head_phos, args.psf)
