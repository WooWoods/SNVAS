"""
    utils module
    ~~~~~~~~~~~~

    A set of reusable instruments.
"""

import os
import re



def dir_check(dirname):
    if os.path.exists(dirname):
        return
    os.makedirs(dirname)

def file_check(filename):
    return os.path.isfile(filename)

def parse_column(val):
    cols = []
    val = str(val)
    val = re.sub(r'\(|\)', '', val)
    for i in val.split(','):
        if '-' in i:
            tmp = i.split('-')
            cols.extend(range(int(tmp[0]) - 1, int(tmp[1])))
        else:
            cols.append(int(i) - 1)
    return cols

def formater_type(values, positions, formater):
    fmt = [formater.normal] * len(values)
    for p in positions:
        try:
            if float(values[p]) <= 0.05:
                fmt[p] = formater.remarkable
        except:
            pass
    return fmt

def print_readme(sheet, readmefile, formater):
    sheet.set_column(0, 0, 30)
    sheet.set_column(1, 1, 60)
    row = 0
    with open(readmefile, 'rt') as fh:
        for line in fh:
            arr = line.strip().split('\t')
            for i, j in enumerate(arr):
                sheet.write(row, i, j, formater.normal)
            row += 1




