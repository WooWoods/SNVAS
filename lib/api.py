import re
import numpy as np
import pandas as pd

from .mathematics import LogitRegression, ChiSquare


def LRanalysis(filename, formula):
    LR = LogitRegression(filename=filename, formula=formula)
    result = LR.gofit()
    print(result.summary2())

    conf = pd.DataFrame(result.conf_int())
    conf['OR'] = result.params
    conf.columns = ['L95', 'U95', 'OR']
    # reindex conf, first element of formula is `y`, so index from 1
    tmpcol = re.split(r'[~+]', formula)
    conf.index = ['Intercept'] + tmpcol[1:]
    res = np.exp(conf)
    res['pvalue'] = result.pvalues
    return res

def Chi_test(filename, y, x):
    chi = ChiSquare(filename=filename, items=(y, x))
    result = chi.calculator()
    return result

