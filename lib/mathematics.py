"""
    mathematics module
    ~~~~~~~~~~~~~~~~~~

    Implements frequently used algorithms.
"""

import os
import re

import pandas as pd
import numpy as np
import patsy
import statsmodels.api as sm
from scipy.stats import chi2_contingency, ttest_ind


class LogitRegression:
    """Logistic regression model.

    :param y: an numpy.array like object, refering to the response variable. Must be submitted along with `X`.
    :param X: similar to `y`, refering to indepent variables.
    :param filename: a <TAB> separated text file, submitted along with `formula`.
    :param formula: a pasty readable string. e.g.
                    'y ~ a + b + a:b'
    """

    def __init__(self, y=None, X=None, filename=None, formula=None):
        self.y = y
        self.X = X
        self.filename = filename
        self.formula = formula

        # to make sure `y` and `X` are None or available simultaneously
        if self.y is not None or self.X is not None:
            if  self.y is None or self.X is None:
                raise('y and X must be provided simultaneously.')

        if any([self.filename, self.formula]) and not all([self.filename, self.formula]):
            raise('file and formula must be provided simultaneously.')

    def data_prepare(self):
        if self.filename is not None:
            dataset = pd.read_table(self.filename, header=0, index_col=None, sep='\t')
            self.y, self.X = patsy.dmatrices(self.formula, dataset)

    def gofit(self):
        self.data_prepare()
        try:
            logit = sm.Logit(self.y, self.X)
            result = logit.fit()
            return result
        except Exception:
            pass


class ChiSquare:
    """Chi square calculator using scipy.stats.chi2_contingency.

    :param filename: a <TAB> delimited txt file.
    :param items: must be submitted with filename to tell what group and
                  variable to be parsed.
                  e.g. ['class', 'sex'] or [0, 1].
    :param dataset: a contingency like dataset. This could be a numpy
                  ndarray object or others. It can be provided solely.
                  e.g. [[100, 97],
                        [84, 127]]
                  mostly this dataset will be a staced pd.series object,
                  e.g.
                  group    gender
                  case     1.0       168
                           2.0        88
                  control  1.0        28
                           2.0        26
                  dtype: int64
    """

    def __init__(self, filename=None, dataset=None, items=None):
        self.filename = filename
        self.dataset = dataset

        if self.filename and self.dataset:
            raise('Accept only one kind of data input.')

        if self.filename is not None:
            if items is None:
                raise('Loss varibles to be analysised. Please refer to the __doc__.')
            table = pd.read_table(self.filename, header=0, index_col=0, sep='\t')
            header = table.columns
            if not contain_item(header, items) and re.search(r'\d', str(items)):
                try:
                    items = list(map(lambda k: header[k], items))
                except IndexError:
                    raise('Items provided not found in the table.')

            self.dataset = table.groupby(items).size().unstack()

    def calculator(self):
        from collections import namedtuple
        Result = namedtuple('Result', 'chi p OR L95 U95')
        chi, p, dof, expected = chi2_contingency(self.dataset)
        OR, L95, U95 = odd_ratio(self.dataset)
        result = Result(chi=chi, p=p, OR=OR, L95=L95, U95=U95)
        # return chi, p, OR, L95, U95
        return result

    def put_down(self, path, var):
        filename = os.path.join(path, '%s_chi.txt' % var)
        result = self.calculator()
        with open(filename ,'wt') as fh:
            fh.write('Chi-score\tp-value\tOR\tL95\tU95\n')
            fh.write('\t'.join(map(lambda x: str(x), result)) + '\n')
            fh.write(str(self.dataset))
        return result

def contain_item(header, items):
    h = set(header)
    i = set(items)
    return i & h == i

def odd_ratio(dataset):
    try:
        darray = np.array(dataset).reshape(2,2)
        OR = darray[0][0] * darray[1][1] / (darray[0][1] * darray[1][0])
        L = np.log(OR)
        SE = np.sqrt(sum(1/darray.reshape(4,1)))[0]
        L95 = np.exp(L - 1.96 * SE)
        U95 = np.exp(L + 1.96 * SE)
        return OR, L95, U95
    except ZeroDivisionError:
        pass


class Ttest:
    """Implement two samples t-test.

    :param vector_one: sample one data.
    :param vector_two: sample two data.
    :filename: a <TAB> delimited txt file, it must be provided along with
          `groupby` and `var`.
    :param groupby: by what column to separate samples into groups.
          could be column number or column name.
    :param var: the variable to be analysised. Same with `groupby`.
    """

    def __init__(self, vector_one=None, vector_two=None, filename=None, groupby=None, var=None):
        self.vector_one = vector_one
        self.vector_two = vector_two
        self.filename = filename
        self.groupby = groupby
        self.var = var

        if self.filename is not None:
            if self.groupby is None or self.var is None:
                raise('`group` and `variable` info not clear,\
                        please refer to the __doc__ for details.')
            table = pd.read_table(self.filename, header=0, index_col=0, sep='\t')
            header = table.columns
            if contain_item(header, [self.groupby, self.var]):
                self.groupby = header.index(self.groupby)
                self.var = header.index(self.var)
            else:
                if not isinstance(self.groupby, int) or not isinstance(self.var, int):
                    raise('`groupby` and `var` incorrect.')
            groupinfo = list(set(table.iloc[:, self.groupby].values))
            self.groupinfo = groupinfo
            if len(groupinfo) > 2:
                raise('`group` not binomial variable.')

            self.vector_one = table.iloc[:, self.var][table.iloc[:, self.var].notnull()][table.iloc[:, self.groupby] == groupinfo[0]].values
            self.vector_two = table.iloc[:, self.var][table.iloc[:, self.var].notnull()][table.iloc[:, self.groupby] == groupinfo[1]].values
        else:
            if not all([self.vector_one, self.vector_two]):
                raise('No effective data input.')

    def calculator(self):
        return ttest_ind(self.vector_one, self.vector_two)

    def summary(self, vector):
        vec = np.array(vector)
        count = len(vec)
        mean_ = vec.mean()
        median_ = median(vec)
        stderr = vec.std()
        min_ = vec.min()
        max_ = vec.max()
        return {'count':count, 'mean': mean_, 'median':median_, 'stde':stderr, 'min':min_, 'max':max_}

    def put_down(self, path, var):
        output = os.path.join(path, '%s_ttest.txt' % var)
        fmt = '{name}\t{count}\t{mean}\t{median}\t{stde}\t{min}\t{max}\n'
        t_score, p = self.calculator()
        with open(output, 'wt') as fh:
            summary_one = self.summary(self.vector_one)
            summary_two = self.summary(self.vector_two)
            if self.groupinfo:
                summary_one['name'] = self.groupinfo[0]
                summary_two['name'] = self.groupinfo[1]
            else:
                summary_one['name'] = 'group 1'
                summary_two['name'] = 'group 2'
            fh.write(fmt.format_map(summary_one))
            fh.write(fmt.format_map(summary_two))
            fh.write('t\t%s' % t_score)
            fh.write('p\t%s' % p)
        return t_score, p, summary_one, summary_two


def median(array):
    length = len(array)
    if length == 1:
        return array[0]
    half = length // 2
    if length % 2 == 0:
        return (array[half - 1] + array[half]) / 2
    else:
        return array[half]


