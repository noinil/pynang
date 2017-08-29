#!/usr/bin/env python3

import numpy as np

class Cafe_ts:
    def __init__(self, ts_file_name):
        self.energy_tags = {}
        self.unit_tags = {}

        col_tag = []
        unit_logs = {}
        with open(ts_file_name, 'r') as ts_fin:
            for line in ts_fin:
                if line.startswith('#unit-unit'):
                    words = line.split()
                    col_tag = words[1:]
                elif line.startswith('#all') or (line.startswith('#') and line[1].isdigit()):
                    row_tag = line[:10]
                    data_str = line[11:].split()
                    local_data = [int(data_str[0])]
                    for d in data_str[1:]:
                        local_data.append(float(d))
                    if row_tag in unit_logs.keys():
                        unit_logs[row_tag].append(local_data[:])
                    else:
                        unit_logs[row_tag] = [local_data[:]]
        for i, t in enumerate(col_tag):
            self.energy_tags[t] = i
        unit_tag_list = list(unit_logs.keys())
        arr_tmp = []
        for i, t in enumerate(unit_tag_list):
            self.unit_tags[t] = i
            arr_unit = np.array(unit_logs[t]).T
            arr_tmp.append(arr_unit)
        self.data = np.array(arr_tmp)

        # ------ Read in summary ------
        print(" ==============================================")
        print(" Read in CafeMol ts file:", ts_file_name)
        nstep = self.data.shape[2]
        for rowt, rlogs in self.unit_tags.items():
            print("  > ", nstep, " steps of ", rowt, " interactions.")
        print("  Each entry has the following columns:")
        for ctags, coli in self.energy_tags.items():
            print("  - {0} : {1}".format(coli, ctags))
        print(" Dimensions of 3d data in this ts file:", self.data.shape)
        print(" ==============================================")

    def get_t_series(self, unit_tag, energy_term):
        unit_index = self.unit_tags[unit_tag.ljust(10, ' ')]
        colm_index = self.energy_tags[energy_term]
        return np.copy(self.data[unit_index, colm_index])


def main():
    pass

if __name__ == '__main__':
    main()

