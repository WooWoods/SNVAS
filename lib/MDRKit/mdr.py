"""
    mdr module
    ~~~~~~~~~~

    Implements mdr and permutation analysis.
"""

import os
import re
import subprocess

import xlsxwriter
from collections import namedtuple

from ..utils import dir_check
from ..xlsx_formater import Formater

class MdrOperate:
    """Gene-gene interaction analysis with Multi Dimensional Reduction method."""
    def __init__(self, asso_inst):
        self.config = asso_inst.config
        self.mdrfile = self.config.get('MDR', None)
        self.path = self.config.get('ROUTINE', None)
        self.tmpdir = self.config.get('TMPDIR', None)
        self.reportdir = os.path.join(self.path, 'report')
        self.resultdir = os.path.join(self.path, 'result/mdr')
        dir_check(self.reportdir)
        dir_check(self.resultdir)

    def go(self):
        result = self.run_mdr()
        models = self.read_mdr_result(result)
        self.to_excel(models)

    def run_mdr(self):
        fmdr = self.mdrfile
        head = open(fmdr).readline()
        max_model = min(3, len(head.split()) - 1)

        mdr_jar = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mdr.jar')
        output = os.path.join(self.resultdir, "mdroutput.txt")

        mdr_command = 'java -jar {0} -min=1 -max={1} -cv=10 -table_data=true -minimal_output=true\
                {2} > {3}'.format(mdr_jar, max_model, fmdr, output)
        commandfile = os.path.join(self.tmpdir, 'mdrun.sh')
        with open(commandfile, 'wt') as fh:
            fh.write(mdr_command)
        subprocess.run(["sh", commandfile])
        return output

    def read_mdr_result(self, result):
        models = []
        wanted_line = False
        with open(result, 'rt') as fh:
            for line in fh:
                if re.match(r'^\s+$', line):
                    continue
                if re.search(r'Finished level', line):
                    wanted_line = True
                    continue
                if not wanted_line:
                    continue
                model, cvc, actrain, actest, *rest = line.strip().split()
                models.append([model, actrain, actest, cvc])
                wanted_line = False
        return models

    def to_excel(self, models):
        workbook = xlsxwriter.Workbook(os.path.join(self.reportdir,'mdr_result.xlsx'))
        formater = Formater(workbook)
        worksheet= workbook.add_worksheet('MDR')
        worksheet.set_row(0, 30)

        header = ('Model', 'bal. acc. CV traning', 'bal. acc. CV testing', 'CV Consistency')
        row = 0
        for i, title in enumerate(header):
            worksheet.write(row, i, title, formater.header)
        row += 1

        for m in models:
            for i, j in enumerate(m):
                worksheet.write(row, i, j, formater.normal)
            row += 1






