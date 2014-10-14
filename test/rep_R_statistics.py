#!/usr/bin/env python

huge_dict = {}
N = 400

def anaf(fin_name):
    coors = []
    global N
    with open(fin_name, 'r') as fin:
        local_coors = []
        for lines in fin:
            if len(lines.split()) == 0 or lines.split()[0] == '#':
                continue
            if "Chain" in lines.split():
                if  len(local_coors) != 0:
                    coors.append(local_coors[:])
                    local_coors.clear()
                continue
            else:
                line = lines.split()
                res_name = line[1]
                if line[3] == "nan" or line[4] == "nan" or line[5] == "nan":
                    continue
                x, y, z = float(line[3]), float(line[4]), float(line[5])
                local_coors.append( ( (x, y, z), res_name) )
        coors.append(local_coors[:])

    # print("Reading", fin_name, "OK!")
    # for chains in coors:
    #     for reses in chains:
    #         print("Res", reses[1], "  Coors:", reses[0])

    n_chain = len(coors)
    for i in range(0, n_chain):
        n_res0 = len(coors[i])
        for r0 in range(0, n_res0):
            res0 = coors[i][r0]
            x0, y0, z0 = res0[0][0], res0[0][1], res0[0][2]
            rname0 = res0[1]
            for j in range(i, n_chain):
                n_res1 = len(coors[j])
                for r1 in range(0, n_res1):
                    res1 = coors[j][r1]
                    if j == i and r1 <= r0:
                        continue
                    x1, y1, z1 = res1[0][0], res1[0][1], res1[0][2]
                    rname1 = res1[1]
                    dx, dy, dz = x0 - x1, y0 - y1, z0 - z1
                    dist01 = (dx * dx + dy * dy + dz * dz) ** 0.5
                    tmp_n = round(dist01 * 2)
                    if rname0 > rname1 or len(rname0) > len(rname1):
                        rn0, rn1 = rname1, rname0
                    else:
                        rn0, rn1 = rname0, rname1
                    sn0 = "INTER" if i < j else "INTRA"
                    sn1 = "TOTAL"
                    # -------------------- Add to huge_dict --------------------
                    thiskey, thiskey1 = (rn0, rn1, sn0), (rn0, rn1, sn1)
                    if thiskey not in huge_dict:
                        huge_dict[thiskey] = [0 for i in range(0, N)]
                    huge_dict[thiskey][tmp_n] += 1
                    if thiskey1 not in huge_dict:
                        huge_dict[thiskey1] = [0 for i in range(0, N)]
                    huge_dict[thiskey1][tmp_n] += 1
                    if len(rn0) == 1:
                        tkey = (rn0, "TOTAL")
                        if tkey not in huge_dict:
                            huge_dict[tkey] = [0 for i in range(0, N)]
                        huge_dict[tkey][tmp_n] += 1
                    else:
                        tkey = ("CA - CA")
                        if tkey not in huge_dict:
                            huge_dict[tkey] = [0 for i in range(0, N)]
                        huge_dict[tkey][tmp_n] += 1
                        tkey1, tkey2 = (rn0, "TOTAL"), (rn1, "TOTAL")
                        if tkey1 not in huge_dict:
                            huge_dict[tkey1] = [0 for i in range(0, N)]
                        huge_dict[tkey1][tmp_n] += 1
                        if tkey2 not in huge_dict:
                            huge_dict[tkey2] = [0 for i in range(0, N)]
                        huge_dict[tkey2][tmp_n] += 1


if __name__ == '__main__':
    # import sys
    # fin_name = sys.argv[1]
    # main(fin_name)
    import numpy as np
    import matplotlib.pyplot as plt
    import glob
    for filename in glob.glob(r'/Users/noinil/Learningspace/*.pos'):
        print(filename)
        # pass
        anaf(filename)

    for i in huge_dict.keys():
        i_distro = huge_dict[i]
        total_amount = sum(i_distro)
        for j in range(0, N):
            huge_dict[i][j] /= total_amount

    pro_name = ['ALA', 'ARG', 'GLU', 'GLY', 'PHE', 'ASP', 'PRO', 'THR', 'SER',\
                'LYS', 'HIS', 'ILE', 'CYS', 'VAL', 'LEU', 'MET', 'TYR', 'TRP',\
                'GLN', 'ASN']
    na_name = ['A', 'T', 'G', 'C', 'P', 'S']

    # OK, output
    # hds = huge_dict.keys()
    # for i in sorted(hds):
    #     print(i, huge_dict[i])

    X = [i/2 for i in range(0, N)]
    XLIM = 10
    # for i in range(0, 20):
    #     fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(6,6))
    #     # out_name = "near_" + pro_name[i] + "_na.png"
    #     out_name = pro_name[i] + "_na.png"
    #     # plt.xlim(0, 10)
    #     for j in range(0, 6):
    #         pn, nn = pro_name[i], na_name[j]
    #         if pn > nn or len(pn) > len(nn):
    #             pn, nn = nn, pn
    #         plotkey1, plotkey2, plotkey0 = (pn, nn, 'INTER'), (pn, nn, 'INTRA'), (pn, nn, 'TOTAL')
    #         Y1 = huge_dict[plotkey1]
    #         m, n = j // 3, j % 3
    #         axes[m, n].set_xlim(0, XLIM)
    #         # axes[m, n].set_ylim(0, 500)
    #         axes[m, n].set_title(nn+'~'+pn, fontsize=10)
    #         axes[m, n].plot(X, Y1)
    #     fig.subplots_adjust(hspace=0.4)
    #     fig.savefig(out_name, dpi=72)

    # fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(6,6))
    # # out_name = "near_nucleic_acid_total.png"
    # out_name = "nucleic_acid_total.png"
    # for i in range(0, 6):
    #     nn = na_name[i]
    #     plotkey = (nn, "TOTAL")
    #     Y1 = huge_dict[plotkey]
    #     m, n = i // 3, i % 3
    #     axes[m, n].set_xlim(0, XLIM)
    #     axes[m, n].set_title(nn+' total', fontsize = 10)
    #     axes[m, n].plot(X, Y1)
    # Y1 = huge_dict[("CA - CA")]
    # axes[2, 0].set_xlim(0, XLIM)
    # axes[2, 0].set_title('CA-CA total', fontsize = 10)
    # axes[2, 0].plot(X, Y1)
    # fig.subplots_adjust(hspace=0.4)
    # fig.savefig(out_name, dpi = 72)

    for i in range(0, 20):
        fig, axes = plt.subplots(nrows=5, ncols=4, figsize=(10,10))
        out_name = "near_" + pro_name[i] + "_pro.png"
        # out_name = pro_name[i] + "_pro.png"
        # plt.xlim(0, 10)
        for j in range(0, 20):
            pn, nn = pro_name[i], pro_name[j]
            if pn > nn:
                pn, nn = nn, pn
            plotkey1, plotkey2, plotkey0 = (pn, nn, 'INTER'), (pn, nn, 'INTRA'), (pn, nn, 'TOTAL')
            Y1 = huge_dict[plotkey1]
            Y2 = huge_dict[plotkey2]
            Y0 = huge_dict[plotkey0]
            m, n = j // 4, j % 4
            axes[m, n].set_xlim(0, XLIM)
            # axes[m, n].set_ylim(0, 500)
            axes[m, n].set_title(nn+'~'+pn, fontsize=10)
            axes[m, n].plot(X, Y1, Y2, Y0)
        fig.subplots_adjust(hspace=0.4)
        fig.savefig(out_name, dpi=72)

    fig, axes = plt.subplots(nrows=5, ncols=4, figsize=(10,10))
    out_name = "near_amino_acid_total.png"
    # out_name = "amino_acid_total.png"
    for i in range(0, 20):
        nn = pro_name[i]
        plotkey = (nn, "TOTAL")
        Y1 = huge_dict[plotkey]
        m, n = i // 4, i % 4
        axes[m, n].set_xlim(0, XLIM)
        axes[m, n].set_title(nn+' total', fontsize = 10)
        axes[m, n].plot(X, Y1)
    fig.subplots_adjust(hspace=0.5)
    fig.savefig(out_name, dpi = 72)
