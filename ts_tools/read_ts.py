#!/usr/bin/env python3

class Cafe_ts:
    def __init__(self, ts_file_name):
        self.column_tags = {}
        self.unit_logs = {}

        with open(ts_file_name, 'r') as ts_fin:
            col_tag = []
            col_num = 0
            for line in ts_fin:
                if line.startswith('#unit-unit'):
                    words = line.split()
                    col_tag = words[1:]
                    col_num = len(col_tag)
                    for i, t in enumerate(col_tag):
                        self.column_tags[t] = i
                elif line.startswith('#all') or (line.startswith('#') and line[1].isdigit()):
                    row_tag = line[:10]
                    data_str = line[11:].split()
                    local_data = [int(data_str[0])]
                    for d in data_str[1:]:
                        local_data.append(float(d))
                    if row_tag in self.unit_logs.keys():
                        self.unit_logs[row_tag].append(local_data[:])
                    else:
                        self.unit_logs[row_tag] = [local_data[:]]

        # ------ Read in summary ------
        print(" ==============================================")
        print(" Read in CafeMol ts file:", ts_file_name)
        for rowt, rlogs in self.unit_logs.items():
            nstep = len(rlogs)
            print("  > ", nstep, " steps of ", rowt, " interactions.")
        print("  Each entry has the following columns:")
        for ctags, coli in self.column_tags.items():
            print("  - {0} : {1}".format(coli, ctags))
        print("\n ==============================================")


def main():
    # r00 = Cafe_ts("r00.ts")
    pass

if __name__ == '__main__':
    main()

