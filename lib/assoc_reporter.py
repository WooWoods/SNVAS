"""
    assoc_reporter module
    ~~~~~~~~~~~~~~~~~~~~~

    Stand alone reporter for only association analysis by plink.
"""

import os
import re
import shutil

import xlsxwriter
from collections import defaultdict
import numpy as np
from .utils import file_check, dir_check, formater_type, print_readme
from .xlsx_formater import Formater


def reporter(assoc_inst):
    hwe_reporter = HweReporter(assoc_inst)
    hwe_reporter.report()
    chi_reporter = ChiReporter(assoc_inst)
    chisq_info_container = chi_reporter.report()
    logit_reporter = LogitReporter(assoc_inst)
    logit_info_container = logit_reporter.report()
    logit_covar_info_container = None
    covar = False
    if assoc_inst.config.get('CORRECTION', None):
        covar = True
        logit_covar_reporter = LogitReporter(assoc_inst, covar)
        logit_covar_info_container = logit_covar_reporter.report()
        if assoc_inst.config.get('PHENO', None):
            logit_pheno_covar_reporter = PhenoLogitReporter(assoc_inst, covar)
            logit_pheno_covar_info_container = logit_pheno_covar_reporter.report()
    if assoc_inst.config.get('PHENO', None):
        logit_pheno_reporter = PhenoLogitReporter(assoc_inst)
        logit_pheno_info_container = logit_pheno_reporter.report()

    all_reporter = AllReport(assoc_inst, chisq_info_container, logit_info_container, covar, logit_covar_info_container)
    all_reporter.report()


class AllReport:
    def __init__(self, assoc_inst, chisq_info_container, logit_info_container, covar=False, *args):
        self.reportdir = os.path.join(assoc_inst.config.get('ROUTINE'), 'report')
        raw_datadir = os.path.join(self.reportdir, 'Raw_data')
        dir_check(raw_datadir)
        tmpdir = os.path.join(assoc_inst.config.get('ROUTINE'), 'tmp')
        shutil.copy(os.path.join(tmpdir, 'sample.map'), raw_datadir)
        shutil.copy(os.path.join(tmpdir, 'sample.ped'), raw_datadir)

        self.report_cutoff = assoc_inst.config.get('REPORT_CUTOFF', None) or 1
        self.chisq = chisq_info_container
        self.logit = logit_info_container
        self.report_covar = False
        if covar:
            self.report_covar = True
            if args[0] is None:
                self.report_covar = False
            self.logit_covar = args[0]

        self.modelname = {
                'dom': 'Dominant',
                'rec': 'Recessive',
                'allele': 'Allele',
                }

    def report(self):
        workbook = xlsxwriter.Workbook(os.path.join(self.reportdir, 'Report.xlsx'))
        formater = Formater(workbook)
        sheet_chi = workbook.add_worksheet('卡方检验')
        sheet_logit = workbook.add_worksheet('逻辑回归')
        if self.report_covar:
            sheet_logit_covar = workbook.add_worksheet('逻辑回归校正')
        header_chi = 'SNP,Class,Model,Genotype,Case,Control,ChiScore,OR(95%CI),P-value,FDR_BH adjusted'.split(',')
        header_logit = 'SNP,Class,Model,Genotype,Case,Control,OR(95%CI),P-value,FDR_BH adjusted'.split(',')

        row_chi = 0
        row_logit = 0
        row_logit_covar = 0
        for snp in self.chisq:
            chi_handler = self.chisq.get(snp)
            if not self.judge_retain(chi_handler):
                continue
            row_chi = self.chi_block_performer(sheet_chi, row_chi, chi_handler, header_chi, formater)
            row_logit = self.logit_block_performer(sheet_logit, row_logit, chi_handler, self.logit, header_logit, formater)
            if self.report_covar:
                row_logit_covar = self.logit_block_performer(sheet_logit_covar, row_logit_covar, self.chisq.get(snp), self.logit_covar, header_logit, formater)
        workbook.close()

    def judge_retain(self, handler):
        try:
            return any(list(map(lambda p: p < self.report_cutoff, handler.allp)))
        except:
            return False

    def logit_block_performer(self, sheet, row, handler, res_container, header, formater):
        orci_fmt = '{0}({1}-{2})'

        for i, j in enumerate(header):
            sheet.write(row, i, j, formater.normal)
        row += 1
        merge_row = row

        blockhandler = ChiBlockHandler(handler)
        logit_handler = res_container.get(handler.snp)
        logitkeys = ['OR', 'L95', 'U95', 'p', 'fdr']
        sheet.merge_range(merge_row, 0, merge_row + 7, 0, handler.snp, formater.normal)
        sheet.merge_range(merge_row, 1, merge_row + 7, 1, 'ALL', formater.normal)
        for key in ['00', '01', '11']:
            arr = blockhandler.codom.get(key)
            if key == '00':
                tmp_arr = arr[0:6] + ['-'] * 3
            elif key == '01':
                tmplogit = [logit_handler.data.get('HET').get(k, '') for k in logitkeys]
                tmp_arr = arr[0:6] + [orci_fmt.format(tmplogit[0], tmplogit[1], tmplogit[2]), tmplogit[3], tmplogit[4]]
            else:
                tmplogit = [logit_handler.data.get('HOM').get(k, '') for k in logitkeys]
                tmp_arr = arr[0:6] + [orci_fmt.format(tmplogit[0], tmplogit[1], tmplogit[2]), tmplogit[3], tmplogit[4]]

            fmt = formater_type(tmp_arr, [7, 8], formater)
            for i, j in enumerate(tmp_arr):
                sheet.write(row, i, j, fmt[i])
            row += 1
        sheet.merge_range(merge_row, 2, merge_row + 2, 2, 'Codominant', formater.normal)
        merge_row += 3

        for model in ['dom', 'rec']:
            tmplogit = [logit_handler.data.get(model.upper()).get(k, '') for k in logitkeys]
            for key in ['0', '1']:
                arr = blockhandler.__dict__.get(model).get(key)
                tmp_arr = arr[0:6] + [orci_fmt.format(tmplogit[0], tmplogit[1], tmplogit[2]), tmplogit[3], tmplogit[4]]
                fmt = formater_type(tmp_arr, [7,8], formater)
                for i, j in enumerate(tmp_arr[3:]):
                    sheet.write(row, i + 3, j, fmt[i+3])
                row += 1
            sheet.merge_range(merge_row, 2, merge_row + 1, 2, self.modelname.get(model), fmt[2])
            sheet.merge_range(merge_row, 6, merge_row + 1, 6, tmp_arr[6], fmt[6])
            sheet.merge_range(merge_row, 7, merge_row + 1, 7, tmp_arr[7], fmt[7])
            sheet.merge_range(merge_row, 8, merge_row + 1, 8, tmp_arr[8], fmt[8])
            merge_row += 2

        tmplogit = [logit_handler.data.get('ADD').get(k, '') for k in logitkeys]
        add_arr = [handler.snp, 'ALL', 'Additive', '-', '-', '-'] + [orci_fmt.format(tmplogit[0], tmplogit[1], tmplogit[2]), tmplogit[3], tmplogit[4]]
        fmt = formater_type(add_arr, [7,8], formater)
        for i, j in enumerate(add_arr):
            sheet.write(row, i, j, fmt[i])
        row += 2
        return row

    def chi_block_performer(self, sheet, row, handler, header, formater):
        for i, j in enumerate(header):
            sheet.write(row, i, j, formater.normal)
        row += 1
        merge_row = row

        blockhandler = ChiBlockHandler(handler)
        sheet.merge_range(merge_row, 0, merge_row + 8, 0, handler.snp, formater.normal)
        sheet.merge_range(merge_row, 1, merge_row + 8, 1, 'ALL', formater.normal)
        for key in ['00', '01', '11']:
            arr = blockhandler.codom.get(key)
            fmt = formater_type(arr, [8,9], formater)
            for i, j in enumerate(arr[3:]):
                sheet.write(row, i + 3, j, fmt[i + 3])
            row += 1
        sheet.merge_range(merge_row, 2, merge_row + 2, 2, 'Codominant', fmt[2])
        sheet.merge_range(merge_row, 6, merge_row + 2, 6, arr[6], fmt[6])
        sheet.merge_range(merge_row, 7, merge_row + 2, 7, arr[7], fmt[7])
        sheet.merge_range(merge_row, 8, merge_row + 2, 8, arr[8], fmt[8])
        sheet.merge_range(merge_row, 9, merge_row + 2, 9, arr[9], fmt[9])
        merge_row += 3

        for model in ['dom', 'rec', 'allele']:
            for key in ['0', '1']:
                arr = blockhandler.__dict__.get(model).get(key)
                fmt = formater_type(arr, [8,9], formater)
                for i, j in enumerate(arr):
                    sheet.write(row, i, j, fmt[i])
                row += 1
            sheet.merge_range(merge_row, 2, merge_row + 1, 2, self.modelname.get(model), fmt[2])
            sheet.merge_range(merge_row, 6, merge_row + 1, 6, arr[6], fmt[6])
            sheet.merge_range(merge_row, 7, merge_row + 1, 7, arr[7], fmt[7])
            sheet.merge_range(merge_row, 8, merge_row + 1, 8, arr[8], fmt[8])
            sheet.merge_range(merge_row, 9, merge_row + 1, 9, arr[9], fmt[9])
            merge_row += 2
        row += 2
        return row

class ChiBlockHandler:
    def __init__(self, handler):
        genofmt = '{0}/{1}'
        self.snp = handler.snp
        self.ref = handler.Majorallele
        self.alt = handler.Minorallele
        self.homr = genofmt.format(self.ref, self.ref)
        self.het = genofmt.format(self.ref, self.alt)
        self.homa = genofmt.format(self.alt, self.alt)
        self.codom = {}
        self.dom = {}
        self.rec = {}
        self.allele = {}

        self.parse_info_codom(handler)
        self.parse_info_dom(handler)
        self.parse_info_rec(handler)
        self.parse_info_allele(handler)

    def parse_info_codom(self, handler):
        genoaff = handler.data['GENO'].get('AFF').split('/')
        genounaff = handler.data['GENO'].get('UNAFF').split('/')
        chi = handler.data['GENO'].get('chi')
        p = handler.data['GENO'].get('p')
        self.codom['00'] = [self.snp, 'ALL', 'Codiminant', self.homr, genoaff[2], genounaff[2], chi, '', p, '']
        self.codom['01'] = [self.snp, 'ALL', 'Codiminant', self.het, genoaff[1], genounaff[1], chi, '', p, '']
        self.codom['11'] = [self.snp, 'ALL', 'Codiminant', self.homa, genoaff[0], genounaff[0], chi, '', p, '']

    def parse_info_dom(self, handler):
        genoaff = handler.data['DOM'].get('AFF').split('/')
        genounaff = handler.data['DOM'].get('UNAFF').split('/')
        chi = handler.data['DOM'].get('chi')
        p = handler.data['DOM'].get('p')
        self.dom['0'] = [self.snp, 'ALL', 'Diminant', self.homr, genoaff[1], genounaff[1], chi, '', p, '']
        self.dom['1'] = [self.snp, 'ALL', 'Diminant', '-'.join([self.het, self.homa]), genoaff[0], genounaff[0], chi, '', p, '']

    def parse_info_rec(self, handler):
        genoaff = handler.data['REC'].get('AFF').split('/')
        genounaff = handler.data['REC'].get('UNAFF').split('/')
        chi = handler.data['REC'].get('chi')
        p = handler.data['REC'].get('p')
        self.rec['0'] = [self.snp, 'ALL', 'Recessive', '-'.join([self.homr,self.het]), genoaff[1], genounaff[1], chi, '', p, '']
        self.rec['1'] = [self.snp, 'ALL', 'Recessive', self.homa, genoaff[0], genounaff[0], chi, '', p, '']

    def parse_info_allele(self, handler):
        genoaff = handler.data['ALLELIC'].get('AFF').split('/')
        genounaff = handler.data['ALLELIC'].get('UNAFF').split('/')
        chi = handler.data['ALLELIC'].get('chi')
        p = handler.data['ALLELIC'].get('p')
        OR = handler.data['ALLELIC'].get('ORCI')
        FDR = handler.data['ALLELIC'].get('FDR')
        self.allele['0'] = [self.snp, 'ALL', 'Allele', self.ref, genoaff[1], genounaff[1], chi, OR, p, FDR]
        self.allele['1'] = [self.snp, 'ALL', 'Allele', self.alt, genoaff[0], genounaff[0], chi, OR, p, FDR]



class HweReporter:
    def __init__(self, assoc_inst):
        self.basepath = assoc_inst.config.get('basepath')
        self.snpinfo = assoc_inst.config.get('SNPFILE')
        self.reportdir = os.path.join(assoc_inst.config.get('ROUTINE'), 'report')
        dir_check(self.reportdir)
        self.resultdir = os.path.join(assoc_inst.config.get('ROUTINE'), 'result')
        self.info_container = {}

    def report(self):
        workbook = xlsxwriter.Workbook(os.path.join(self.reportdir, 'HWE.xlsx'))
        formater = Formater(workbook)
        sheet = workbook.add_worksheet('HWE')
        sheet.set_row(0, 30)
        sheet_readme = workbook.add_worksheet('ReadMe')
        readmefile = os.path.join(self.basepath, 'ReadMetxt/readme_hwe.txt')
        print_readme(sheet_readme, readmefile, formater)

        self.record_hwe_result()
        self.record_maf()
        self.parse_annotation()
        header = 'SNP,CHR,Position(hg19),Minor allele,Major allele,GeneName,Mrna,Region,\
                CHBS_1000g,Total(11/01/00),Total MAF,HWE,Case(11/01/00),\
                Case_majorallele_number,Case_minorallele_number,Case MAF,HWE_Case,\
                Control(11/01/00),Control_majorallele_number,Control_minorallele_number,\
                Control MAF,HWE_Control'.split(',')
        row = 0
        for i, j in enumerate(map(lambda s: s.strip(), header)):
            sheet.write(row, i, j, formater.header)
        row += 1

        for snp in self.info_container:
            handler = self.info_container.get(snp)
            line = handler.output()
            fmt = formater_type(line, [11, 16, -1], formater)
            for i, v in enumerate(line):
                sheet.write(row, i, v, fmt[i])
            row += 1
        workbook.close()

    def record_hwe_result(self):
        snpinfo = self.get_snpinfo()
        hwefile = os.path.join(self.resultdir, 'hwe/hwe.hwe')
        count = 0
        with open(hwefile, 'rt') as fh:
            for line in fh:
                count += 1
                if count == 1:
                    continue
                arr = line.split()
                if arr[1] in self.info_container:
                    handler = self.info_container.get(arr[1])
                else:
                    handler = HweHandler(arr[1])
                    self.info_container[arr[1]] = handler
                handler.Chr = arr[0]
                handler.pos = snpinfo[arr[1]][0]
                handler.Minorallele = arr[3]
                handler.Majorallele = arr[4]
                handler.add_info(arr)
    def get_snpinfo(self):
        snp = {}
        with open(self.snpinfo, 'rt') as fh:
            for line in fh:
                if re.match(r'^\s$', line):
                    continue
                arr = line.split()
                snp[arr[0]] = arr[2:5]
        return snp

    def record_maf(self):
        maffile = os.path.join(self.resultdir, 'hwe/freq.frq')
        ccmaffile = os.path.join(self.resultdir, 'hwe/freq.frq.cc')

        count = 0
        with open(maffile, 'rt') as fh:
            for line in fh:
                count += 1
                if count == 1:
                    continue
                arr = line.split()
                if arr[1] in self.info_container:
                    handler = self.info_container.get(arr[1])
                else:
                    continue
                handler.data.get('ALL')['maf'] = arr[4]

        count = 0
        with open(ccmaffile, 'rt') as fh:
            for line in fh:
                count += 1
                if count == 1:
                    continue
                arr = line.split()
                if arr[1] in self.info_container:
                    handler = self.info_container.get(arr[1])
                else:
                    continue
                handler.data.get('AFF')['maf'] = arr[4]
                handler.data.get('UNAFF')['maf'] = arr[5]

    def parse_annotation(self):
        f1000g = os.path.join(self.resultdir, 'hwe/library.hg19_ALL.sites.2012_02_dropped')
        fgeneanno = os.path.join(self.resultdir, 'hwe/library.variant_function')
        fmrnaanno = os.path.join(self.resultdir, 'hwe/library.exonic_variant_function')

        try:
            with open(f1000g, 'rt') as fh:
                for line in fh:
                    arr = line.strip().split()
                    if arr[-1] in self.info_container:
                        handler = self.info_container.get(arr[-1])
                    else:
                        continue
                    handler.g1000 = arr[1]
        except FileNotFoundError:
            pass
        try:
            with open(fgeneanno, 'rt') as fh:
                for line in fh:
                    arr = line.strip().split()
                    if arr[-1] in self.info_container:
                        handler = self.info_container.get(arr[-1])
                    else:
                        continue
                    handler.region = arr[0]
                    gene = re.match(r'^(\w+)', arr[1]).group(1)
                    handler.gene = gene
        except FileNotFoundError:
            pass
        try:
            with open(fmrnaanno, 'rt') as fh:
                for line in fh:
                    arr = line.strip().split('\t')
                    if arr[-1] in self.info_container:
                        handler = self.info_container.get(arr[-1])
                    else:
                        continue
                    mrnainfo = re.split(r'[:,]', arr[2])
                    mrnas = list(filter(lambda x: re.match(r'^NM', x), mrnainfo))
                    handler.mrna = ','.join(mrnas)
        except FileNotFoundError:
            pass


class HweHandler:
    """Hwe result handler, to save result for each snv."""
    def __init__(self, snp):
        self.snp = snp
        self.Chr = None
        self.Minorallele = None
        self.Majorallele = None
        self.data = {}

    def add_info(self, arr):
        test, minor, major, geno, *rest, p = arr[2:]
        self.data.setdefault(test, {}).__setitem__('geno', geno)
        self.data.setdefault(test, {}).__setitem__('p', p)

    def output(self):
        line_arr = [self.snp, self.Chr, self.pos, self.Minorallele, self.Majorallele]
        tmp = [self.__dict__.get(key, '') for key in ['gene', 'mrna', 'region', 'g1000']]
        line_arr.extend(tmp)
        alltmp = self.dic_info_list(self.data['ALL'], 'ALL')
        line_arr.extend(alltmp)
        afftmp = self.dic_info_list(self.data['AFF'], 'AFF')
        line_arr.extend(afftmp)
        unafftmp = self.dic_info_list(self.data['UNAFF'], 'UNAFF')
        line_arr.extend(unafftmp)
        return line_arr

    @staticmethod
    def dic_info_list(dic, key):
        if key == 'ALL':
            return [dic.get(key, '') for key in ['geno', 'maf', 'p']]
        try:
            geno = list(map(lambda s: int(s), dic.get('geno', '').split('/')))
            refnum = geno[-1] * 2 + geno[1]
            altnum = geno[0] * 2 + geno[1]
        except:
            refnum = 0
            altnum = 0
        tmp = [dic.get(key, '') for key in ['geno', 'maf', 'p']]
        tmp.insert(1, altnum)
        tmp.insert(1, refnum)
        return tmp


class ChiReporter:
    """Put chi-square analysis result into xlsx files."""
    def __init__(self, assoc_inst):
        self.basepath = assoc_inst.config.get('basepath')
        self.reportdir = os.path.join(assoc_inst.config.get('ROUTINE'), 'report')
        dir_check(self.reportdir)
        self.resultdir = os.path.join(assoc_inst.config.get('ROUTINE'), 'result')
        self.info_container = {}

    def report(self):
        workbook = xlsxwriter.Workbook(os.path.join(self.reportdir, 'ChiScore.xlsx'))
        formater = Formater(workbook)
        sheet = workbook.add_worksheet('ALL')
        sheet.set_row(0, 30)
        sheet_readme = workbook.add_worksheet('ReadMe')
        readmefile = os.path.join(self.basepath, 'ReadMetxt/readme_chi.txt')
        print_readme(sheet_readme, readmefile, formater)

        self.record_model_result()
        self.record_assoc_result()

        header = 'SNP,CHR,Major allele,Minor allele,Model,AFF(11|10|00),\
                UNAFF(11|10|00),ChiScore,OR(95%CI),P-value,FDR_BH adjusted'.split(',')
        row = 0
        for i, j in enumerate(map(lambda s: s.strip(),header)):
            sheet.write(row, i, j, formater.header)
        row += 1

        for snp in self.info_container:
            handler = self.info_container.get(snp)
            lines = handler.output()
            for line in lines:
                fmt = formater_type(line, [9, 10], formater)
                for i, j in enumerate(line):
                    sheet.write(row, i, j, fmt[i])
                row += 1
        workbook.close()
        return self.info_container

    def record_model_result(self):
        modelfile = os.path.join(self.resultdir, 'chi-test/model_chi.model')
        count = 0
        with open(modelfile, 'rt') as fh:
            for line in fh:
                count += 1
                if count == 1:
                    continue
                arr = line.split()
                if arr[4] == 'TREND':
                    continue
                snp = arr[1]
                if snp in self.info_container:
                    handler = self.info_container.get(snp)
                else:
                    handler = ChiHandler(snp)
                    self.info_container[snp] = handler
                handler.Chr = arr[0]
                handler.Minorallele = arr[2]
                handler.Majorallele = arr[3]
                handler.add_info(arr)

    def record_assoc_result(self):
        assocfile = os.path.join(self.resultdir, 'chi-test/chi.assoc')
        adjustfile = os.path.join(self.resultdir, 'chi-test/chi.assoc.adjusted')
        count = 0
        with open(assocfile, 'rt') as fh:
            for line in fh:
                count += 1
                if count == 1:
                    continue
                arr = line.split()
                snp = arr[1]
                if snp in self.info_container:
                    handler = self.info_container.get(snp)
                else:
                    continue
                handler.data['ALLELIC']['ORCI'] = '%s(%s-%s)' %(arr[9], arr[11], arr[12])

        count = 0
        with open(adjustfile, 'rt') as fh:
            for line in fh:
                count += 1
                if count == 1:
                    continue
                arr = line.split()
                snp = arr[1]
                if snp in self.info_container:
                    handler = self.info_container.get(snp)
                else:
                    continue
                handler.data['ALLELIC']['FDR'] = arr[-2]


class ChiHandler:
    """To save result of chi-square analysis."""
    def __init__(self, snp):
        self.snp = snp
        self.Chr = None
        self.Minorallele = None
        self.Majorallele = None
        self.data = {}
        self.allp = []

    def output(self):
        base_arr = [self.snp, self.Chr, self.Majorallele, self.Minorallele]
        keys = ['AFF', 'UNAFF', 'chi', 'ORCI', 'p', 'FDR']
        codom = base_arr + ['Codominant'] + [self.data['GENO'].get(key, '') for key in keys]
        dom = base_arr + ['Dominant'] + [self.data['DOM'].get(key, '') for key in keys]
        rec = base_arr + ['Recessive'] + [self.data['REC'].get(key, '') for key in keys]
        allele = base_arr + ['Allele'] + [self.data['ALLELIC'].get(key, '') for key in keys]
        return (codom, dom, rec, allele)

    def add_info(self, arr):
        test, affgeno, unaffgeno, chisq, df, p = arr[4:]
        self.data.setdefault(test, {}).__setitem__('AFF', affgeno)
        self.data.setdefault(test, {}).__setitem__('UNAFF', unaffgeno)
        self.data.setdefault(test, {}).__setitem__('chi', chisq)
        self.data.setdefault(test, {}).__setitem__('p', p)
        if re.match(r'\d', p):
            self.allp.append(float(p))


class LogitReporter:
    """Put result of logistic analysis into a xlsx."""
    def __init__(self, assoc_inst, covar=False):
        self.basepath = assoc_inst.config.get('basepath')
        self.reportdir = os.path.join(assoc_inst.config.get('ROUTINE'), 'report')
        dir_check(self.reportdir)
        self.resultdir = os.path.join(assoc_inst.config.get('ROUTINE'), 'result/logistic-test')
        self.info_container = {}
        self.report_covar = covar
        if self.report_covar:
            self.resultdir = os.path.join(assoc_inst.config.get('ROUTINE'), 'result/logistic-test/logit_covar')

    def report(self):
        self.iter_models()
        if self.report_covar:
            workbook = xlsxwriter.Workbook(os.path.join(self.reportdir, 'Logistic_CORRECT.xlsx'))
        else:
            workbook = xlsxwriter.Workbook(os.path.join(self.reportdir, 'Logistic.xlsx'))
        formater = Formater(workbook)
        sheet = workbook.add_worksheet('ALL')
        sheet.set_row(0, 30)
        sheet_readme = workbook.add_worksheet('ReadMe')
        readmefile = os.path.join(self.basepath, 'ReadMetxt/readme_logit.txt')
        print_readme(sheet_readme, readmefile, formater)

        header = 'SNP,CHR,BP,Alt Allele,Model,NMISS,OR,SE,L95,U95,STAT,P-value,FDR_BH adjusted'.split(',')
        row = 0
        for i, j in enumerate(header):
            sheet.write(row, i, j, formater.header)
        row += 1

        for snp in self.info_container:
            handler = self.info_container.get(snp)
            lines = handler.output()
            for line in lines:
                fmt = formater_type(line, [11,12], formater)
                for i, j in enumerate(line):
                    sheet.write(row, i, j, fmt[i])
                row += 1
        workbook.close()
        return self.info_container

    def iter_models(self):
        models = [('dominant', 'DOM'), ('recessive', 'REC'), ('', 'ADD'), ('hethom', 'HETHOM')]
        for model, mark in models:
            logitfile = os.path.join(self.resultdir, 'logistic%s.assoc.logistic' % model)
            self.record_logit_result(logitfile, mark)

    def record_logit_result(self, filename, mark):
        """extracting result  from logistic analysis result files
        of different genetic models.
        """
        adjusted = filename + '.adjusted'
        try:
            with open(filename, 'rt') as fh:
                for line in fh:
                    arr = line.split()
                    snp = arr[1]
                    if arr[4] in ('ADD', 'DOM', 'REC', 'HOM', 'HET'):
                        if snp in self.info_container:
                            handler = self.info_container.get(snp)
                        else:
                            handler = LogitHandler(snp)
                            self.info_container[snp] = handler
                        handler.Chr = arr[0]
                        handler.pos = arr[2]
                        handler.Minorallele = arr[3]
                        handler.add_info(arr)
        except FileNotFoundError:
            pass

        count = 0
        try:
            with open(adjusted, 'rt') as fh:
                for line in fh:
                    count += 1
                    if count == 1:
                        continue
                    arr = line.split()
                    snp = arr[1]
                    handler = self.info_container.get(snp)
                    if mark == 'HETHOM':
                        handler.data.get('HET').__setitem__('fdr', arr[-2])
                        handler.data.get('HOM').__setitem__('fdr', arr[-2])
                    else:
                        handler.data.get(mark).__setitem__('fdr', arr[-2])
        except FileNotFoundError:
            pass


class LogitHandler:
    """To save logistic analysis result for each snv."""
    def __init__(self, snp):
        self.snp = snp
        self.Chr = None
        self.pos = None
        self.Minorallele = None
        self.Majorallele = None
        self.data = {}

    def add_info(self, arr):
        test, nmiss, OR, SE, L95, U95, stat, p = arr[4:]
        self.data.setdefault(test, {}).__setitem__('nmiss', nmiss)
        self.data.setdefault(test, {}).__setitem__('OR', OR)
        self.data.setdefault(test, {}).__setitem__('SE', SE)
        self.data.setdefault(test, {}).__setitem__('L95', L95)
        self.data.setdefault(test, {}).__setitem__('U95', U95)
        self.data.setdefault(test, {}).__setitem__('stat', stat)
        self.data.setdefault(test, {}).__setitem__('p', p)

    def output(self):
        base_arr = [self.snp, self.Chr, self.pos, self.Minorallele]
        keys = ['nmiss', 'OR', 'SE', 'L95', 'U95', 'stat', 'p', 'fdr']
        dom = base_arr + ['Dominant'] + [self.data['DOM'].get(key, '') for key in keys]
        rec = base_arr + ['Recessive'] + [self.data['REC'].get(key, '') for key in keys]
        add = base_arr + ['Additive'] + [self.data['ADD'].get(key, '') for key in keys]
        hom = base_arr + ['HOM'] + [self.data['HOM'].get(key, '') for key in keys]
        het = base_arr + ['HET'] + [self.data['HET'].get(key, '') for key in keys]
        return (dom, rec, add, hom, het)


class PhenoLogitReporter(LogitReporter):
    """Put result of logistic analysis for phenotypes and genotypes into a xlsx."""
    def __init__(self, assoc_inst, covar=False):
        self.basepath = assoc_inst.config.get('basepath')
        self.reportdir = os.path.join(assoc_inst.config.get('ROUTINE'), 'report')
        dir_check(self.reportdir)
        self.resultdir = os.path.join(assoc_inst.config.get('ROUTINE'), 'result/logistic-test/phenoassoc')
        self.info_container = {}
        self.report_covar = covar
        if self.report_covar:
            self.resultdir = os.path.join(assoc_inst.config.get('ROUTINE'), 'result/logistic-test/phenoassoc_covar')

    def report(self):
        self.iter_models()
        readmefile = os.path.join(self.basepath, 'ReadMetxt/readme_phenologit.txt')
        if self.report_covar:
            workbook = xlsxwriter.Workbook(os.path.join(self.reportdir, 'PhenoLogistic_CORRECT.xlsx'))
        else:
            workbook = xlsxwriter.Workbook(os.path.join(self.reportdir, 'PhenoLogistic.xlsx'))
        formater = Formater(workbook)
        sheet = workbook.add_worksheet('ALL')
        sheet.set_row(0, 30)
        sheet_readme = workbook.add_worksheet('ReadMe')
        print_readme(sheet_readme, readmefile, formater)
        header = 'PhenoName,SNP,CHR,BP,Alt Allele,Model,NMISS,Beta,SE,L95,U95,STAT,P-value,FDR_BH adjusted'.split(',')
        row = 0
        for i, j in enumerate(header):
            sheet.write(row, i, j, formater.header)
        row += 1

        for pheno_snp in self.info_container:
            handler = self.info_container.get(pheno_snp)
            line_head = pheno_snp.split('-') # turn joined snp-pheno into list: [pheno, snp]
            lines = handler.output()
            for line in lines:
                arr = line_head + line[1:]
                fmt = formater_type(arr, [12,13], formater)
                for i, j in enumerate(arr):
                    sheet.write(row, i, j, fmt[i])
                row += 1
        workbook.close()
        return self.info_container

    def iter_models(self):
        import glob
        models_mark = dict([('dominant', 'DOM'), ('recessive', 'REC'), ('add', 'ADD'), ('hethom', 'HETHOM')])
        for filename in glob.glob('%s/*linear' % self.resultdir):
            self.record_logit_result(filename, models_mark)

    def record_logit_result(self, filename, models_mark):
        adjusted = filename + '.adjusted'
        basename_info = os.path.basename(filename).split('.')
        pheno_name = basename_info[1]
        model = basename_info[0].split('_')[1]
        mark = models_mark.get(model)

        try:
            with open(filename, 'rt') as fh:
                for line in fh:
                    arr = line.split()
                    pheno_snp = '-'.join([pheno_name, arr[1]]) # to accommodate this specific situation to `LogitHandler`
                    if arr[4] in ('ADD', 'DOM', 'REC', 'HOM', 'HET'):
                        if pheno_snp in self.info_container:
                            handler = self.info_container.get(pheno_snp)
                        else:
                            handler = LogitHandler(pheno_snp)
                            self.info_container[pheno_snp] = handler
                        handler.Chr = arr[0]
                        handler.pos = arr[2]
                        handler.Minorallele = arr[3]
                        handler.add_info(arr)
        except FileNotFoundError:
            pass

        count = 0
        try:
            with open(adjusted, 'rt') as fh:
                for line in fh:
                    count += 1
                    if count == 1:
                        continue
                    arr = line.split()
                    pheno_snp = '-'.join([pheno_name, arr[1]])
                    handler = self.info_container.get(pheno_snp)
                    if mark == 'HETHOM':
                        handler.data.get('HET').__setitem__('fdr', arr[-2])
                        handler.data.get('HOM').__setitem__('fdr', arr[-2])
                    else:
                        handler.data.get(mark).__setitem__('fdr', arr[-2])
        except FileNotFoundError:
            pass










