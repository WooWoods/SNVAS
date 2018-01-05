"""
    pheno independ test module
    ~~~~~~~~~~~~~~~~~~~~~~~~~~

    Implements independ test among groups for pheno.
"""

import os
import re

import xlsxwriter
import numpy as np
import pandas as pd

from . import ChiSquare, Ttest
from .utils import dir_check, parse_column, formater_type
from .xlsx_formater import Formater




class PhenoIndepTest:
    def __init__(self, asso_inst):
        self.config = asso_inst.config
        self.info_file = self.config.get('INFOFILE', None)
        self.pheno_chi = self.config.get('CHI_TEST', None)
        self.pheno_ttest = self.config.get('TTEST', None)
        self.path = self.config.get('ROUTINE', None)

        self.resultdir = os.path.join(self.path, 'result/pheno_test')
        self.reportdir = os.path.join(self.path, 'report')
        dir_check(self.resultdir)
        dir_check(self.reportdir)

        self.t_result_container = []
        self.chi_result_container = []

    def go(self):
        self.iter_test()
        self.to_excel()

    def iter_test(self):
        header = open(self.info_file).readline().strip().split('\t')
        if self.pheno_chi is not None and re.search(r'\d', str(self.pheno_chi)):
            cols = parse_column(self.pheno_chi)
            for col in cols:
                var_name = header[col + 1]
                result, dataset, cata = self.var_chisq([0, col], var_name)
                chi, p, OR, L95, U95 = result
                chi_result = ChitestHandler(var_name, dataset, chi, p, cata)
                self.chi_result_container.append(chi_result)

        if self.pheno_ttest is not None and re.search(r'\d', str(self.pheno_ttest)):
            cols = parse_column(self.pheno_ttest)
            for col in cols:
                var_name = header[col + 1]
                tvalue, p, sum_one, sum_two = self.var_ttest(col, var_name)
                t_result = TtestHandler(var_name, tvalue, p, sum_one, sum_two)
                self.t_result_container.append(t_result)

    def var_chisq(self, var, var_name):
        table = pd.read_table(self.info_file, header=0, index_col=0, sep='\t')
        table.replace('-9', np.nan, inplace=True)
        header = table.columns
        if not contain_item(header, var) and re.search(r'\d', str(var)):
            try:
                var = list(map(lambda k: header[k], var))
            except IndexError:
                raise Exception('Items provided not found in the table.')
        dataset = table.groupby(var).size().unstack()

        chi = ChiSquare(dataset=np.array(dataset))
        result = chi.put_down(self.resultdir, var_name)
        return result, chi.dataset, list(dataset.columns)

    def var_ttest(self, var, var_name, group=0):
        t_test = Ttest(filename=self.info_file, groupby=group, var=var)
        tvalue, p, sum_one, sum_two = t_test.put_down(self.resultdir, var_name)
        return tvalue, p, sum_one, sum_two

    def to_excel(self):
        workbook = xlsxwriter.Workbook(os.path.join(self.reportdir, 'PhenoTest.xlsx'))
        formater = Formater(workbook)
        sheet_chi = workbook.add_worksheet('Chi-test')
        sheet_t = workbook.add_worksheet('T-test')
        row_chi = 0
        row_t = 0
        for result in self.t_result_container:
            row_t = self.t_result_printer(result, sheet_t, row_t, formater)
        for result in self.chi_result_container:
            row_chi = self.chi_result_printer(result, sheet_chi, row_chi, formater)
        workbook.close()

    def chi_result_printer(self, result, sheet, row, formater):
        sheet.write(row, 0, result.item_name, formater.header)
        row += 1

        for i, j in enumerate(['项目'] + result.cata):
            sheet.write(row, i, str(j), formater.header)
        row += 1
        case_data = ['case'] + list(np.array(result.dataset).reshape(2,-1)[0])
        ctl_data = ['control'] + list(np.array(result.dataset).reshape(2,-1)[1])
        for i, j in enumerate(case_data):
            sheet.write(row, i, str(j), formater.normal)
        row += 1
        for i, j in enumerate(ctl_data):
            sheet.write(row, i, str(j), formater.normal)
        row += 1
        sheet.write(row, 0, 'case--control', formater.normal)
        row += 1
        sheet.write(row, 0, 'chi-score', formater.normal)
        sheet.write(row, 1, str(result.chi), formater.normal)
        row += 1
        sheet.write(row, 0, 'p', formater.normal)
        p = result.p
        fmt = formater.normal
        if p <= 0.05:
            fmt = formater.remarkable
        sheet.write(row, 1, str(result.p), fmt)
        row += 1
        return row

    @staticmethod
    def parse_nan(datalist):
        return ['NA' if i is np.nan else i for i in datalist]

    def t_result_printer(self, result, sheet, row, formater):
        header = ['项目', '例数', '均值', '中位数', '标准差', '最小值', '最大值']
        keys = ['count', 'mean', 'median', 'stde', 'min', 'max']
        sheet.write(row, 0, result.item_name, formater.header)
        row += 1
        for i, j in enumerate(header):
            sheet.write(row, i, j, formater.header)
        row += 1
        case_data = ['case'] + [result.summary_one.get(key, '') for key in keys]
        ctl_data = ['control'] + [result.summary_two.get(key, '') for key in keys]
        for i, j in enumerate(case_data):
            sheet.write(row, i, j, formater.normal)
        row += 1
        for i, j in enumerate(ctl_data):
            sheet.write(row, i, j, formater.normal)
        row += 1
        sheet.write(row, 0, 'case--control', formater.normal)
        row += 1
        sheet.write(row, 0, 't', formater.normal)
        sheet.write(row, 1, result.tvalue, formater.normal)
        row += 1
        sheet.write(row, 0, 'p', formater.normal)
        p = result.p
        fmt = formater.normal
        if p <= 0.05:
            fmt = formater.remarkable
        sheet.write(row, 1, result.p, fmt)
        row += 1
        return row


class TtestHandler:
    def __init__(self, item_name, tvalue, p, summary_one, summary_two):
        self.item_name = item_name
        self.tvalue = tvalue
        self.p = p
        self.summary_one = summary_one
        self.summary_two = summary_two

class ChitestHandler:
    def __init__(self, item_name, dataset, chi, p, cata):
        self.item_name = item_name
        self.dataset = dataset
        self.chi = chi
        self.p = p
        self.cata = cata


def contain_item(header, items):
    h = set(header)
    i = set(items)
    return i & h == i

