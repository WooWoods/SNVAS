"""
    pheno independ test module
    ~~~~~~~~~~~~~~~~~~~~~~~~~~

    Implements independ test among groups for pheno.
"""

import os
import re

import xlsxwriter

from . import ChiSquare, Ttest
from .utils import dir_check, parse_column



class PhenoIndepTest:
    def __init__(self, asso_inst):
        self.config = asso_inst
        self.info_file = self.config.get('INFOFILE', None)
        self.pheno_chi = self.config.get('Chi-test', None)
        self.pheno_ttest = self.config.get('T-test', None)
        self.path = self.config.get('ROUTINE', None)

        self.resultdir = os.path.join(self.path, 'pheno_test')
        self.reportdir = os.path.join(self.path, 'report')
        dir_check(self.resultdir)
        dir_check(self.reportdir)

    def iter_test(self):
        header = open(self.info_file).readline().strip().split()
        if self.pheno_chi is not None and re.search(r'\d', str(self.pheno_chi)):
            cols = parse_column(self.pheno_chi)
            for col in cols:
                var_name = header[col]
                var_chisq([1, col], var_name)

        if self.pheno_ttest is not None and re.search(r'\d', str(self.pheno_ttest)):
            cols = parse_column(self.pheno_ttest)
            for col in cols:
                var_name = header[col]
                var_ttest(col, var_name)

    def var_chisq(self, var_name):
        chi = ChiSquare(filename=self.info_file, items=var)
        chi.put_down(self.resultdir, var_name)

    def var_ttest(self, group=1, var_name):
        t_test = Ttest(filename=self.info_file, group=group, var=var)
        t_test.put_down(self.resultdir, var_name)

    def to_excel(self):
        workbook = xlsxwriter.add_workbook(os.path.join(self.reportdir, 'PhenoTest.xlsx'))
        worksheet1 = workbook.add_sheet('Chi-test')
        worksheet2 = workbook.add_sheet('T-test')






