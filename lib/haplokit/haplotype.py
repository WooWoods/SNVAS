"""
    haplokit module
    ~~~~~~~~~~~~~~~~

    Implements haplotype analysis.
"""

import os
import re
import subprocess
from copy import deepcopy

import pandas as pd
import numpy as np
import xlsxwriter
import patsy
from collections import defaultdict, UserDict

from ..utils import dir_check, file_check, parse_column, formater_type, print_readme
from ..mathematics import LogitRegression
from ..xlsx_formater import Formater
from .block_read import BlockIdentifier

def hap_analysis(assoc_inst):
    split_ped = SplitPed(assoc_inst)
    genes = split_ped.go()
    haploview = Haploview(assoc_inst, genes)
    haploview.go()
    hapassoc = HapAssocAnalysis(assoc_inst)
    hapassoc.hap_go()


class SplitPed:
    """Split sample.ped into pieces according to genes, i.e. one
    gene, one ped file. These ped files will be submmited to haloview.
    """
    def __init__(self, assoc_inst):
        self.config = assoc_inst.config
        self.path = self.config.get('ROUTINE', None)
        self.basepath = self.config.get('basepath', None)
        self.hapfile = self.config.get('RAW_HAP', None) or \
                os.path.join(self.path, 'data/raw_hap.txt')
        self.tmpdir = self.config.get('TMPDIR', None)\
                or os.path.join(self.path, 'tmp')
        self.pedfile = os.path.join(self.tmpdir, 'sample.ped')
        self.mapfile = os.path.join(self.tmpdir, 'sample.map')

        self.header = ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO']

    def go(self):
        ped, map_ = self.load_table()
        hap = self.parse_hap()
        self.gene_ped(hap, ped, map_)
        genes = hap.keys()
        return genes

    def load_table(self):
        ped = pd.read_table(self.pedfile, header=None, index_col=None, sep='\t')
        map_= pd.read_table(self.mapfile, header=None, index_col=None, sep='\t')

        snps = map_.iloc[:, 1].values
        header = self.header + list(snps)
        ped.columns = header
        return ped, map_.iloc[:, [1,3]]

    def parse_hap(self):
        hap = {}
        with open(self.hapfile, 'rt') as fh:
            for line in fh:
                arr = line.strip().split()
                hap[arr[1]] = arr[2:]
        return hap

    def gene_ped(self, hap, ped, map_):
        """select corresponding columns from sample.ped
        file by snvs of a gene.
        """
        dirname = os.path.join(self.path, 'report/Raw_data')
        dir_check(dirname)
        for gene in hap:
            tmpped = os.path.join(dirname, '{0}.ped'.format(gene))
            tmpinfo = os.path.join(dirname, '{0}.info'.format(gene))
            snps = hap[gene]
            header = self.header + snps
            geneped = ped[header]
            genemap = map_.loc[map_.iloc[:, 0].isin(snps)]
            genemap.to_csv(tmpinfo, header=False, index=False, sep='\t')

            for n in geneped.columns[6:]:
                col = geneped[n]
                self.replace_indel(col)
            geneped.to_csv(tmpped, header=False, index=False, sep='\t')

    @staticmethod
    def replace_zero(dataframe):
        """missing genotypes are expressed by 'N N' in haploview."""
        _ = dataframe.replace('0 0', 'N N', inplace=True)

    @staticmethod
    def replace_indel(col):
        """Since 'haploview' software recognise noly ['A', 'C', 'C', 'G'] single bases,
        indels expressed with '-' and 'AAA' like oligos need to be replace by single
        base in case the program breaks.
        """
        bases = set(['A', 'T', 'C', 'G'])
        if '-' not in ''.join(col.values):
            return
        col_bases = set([i for pair in col.values for i in pair.split() if re.match(r'[A-Za-z]', i, re.I)])
        col_bases_index0 = set([i[0] for i in col_bases])
        replacement = random.choice(list(bases - col_bases_index0))
        base = list(col_bases)[0]
        homa = '- -'
        het = '%s -' % base
        het2 = '- %s' % base
        homr = '%s %s' %(base, base)
        _ = col.replace(homa, '%s %s' %(replacement, replacement), inplace=True)
        _ = col.replace(het, '%s %s' %(base[0], replacement), inplace=True)
        _ = col.replace(het2, '%s %s' %(replacement, base[0]), inplace=True)
        _ = col.replace(homr, '%s %s' %(base[0], base[0]), inplace=True)


class Haploview:
    """Haploview plot and LD calculation."""
    def __init__(self, assoc_inst, genes):
        self.config = assoc_inst.config
        self.path = self.config.get('ROUTINE', None)
        self.reportdir = os.path.join(self.path, 'report')
        self.genes = genes

    def go(self):
        files = self.read_dir()
        self.haploview(files)
        self.LD_block_xlsx()

    def read_dir(self):
        import glob
        peddir = os.path.join(self.reportdir, 'Raw_data')
        return glob.glob('%s/*ped' % peddir)

    def outpath(self):
        D = os.path.join(self.reportdir, 'haploview/D_Prime')
        R2 = os.path.join(self.reportdir, 'haploview/R_Square')
        dir_check(D)
        dir_check(R2)
        return D, R2

    def haploview(self, pedlist):
        hap_jar = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Haploview.jar')
        D, R2 = self.outpath()
        for fped in pedlist:
            finfo = re.sub(r'ped', 'info', fped)
            gene = re.sub(r'\.ped', '', os.path.basename(fped))
            for to_dir, ldvalues in ((R2, 'RSQ'), (D, 'DPRIME')):
                out_item = os.path.join(to_dir, gene + '_' + ldvalues)
                subprocess.run(["java", "-jar", hap_jar, "-n", "-out", out_item,
                                "-pedfile", fped, "-info", finfo,
                                "-png", "-ldcolorscheme", "GOLD", "-ldvalues",
                                ldvalues, "-blockoutput", "GAB"])
            subprocess.run(["java", "-jar", hap_jar, "-n", "-out",
                            os.path.join(self.reportdir, 'haploview/%s' %gene),
                            "-pedfile", fped, "-info", finfo,
                            "-dprime", "-blockoutput", "GAB"])

    def LD_block_xlsx(self):
        workbook = xlsxwriter.Workbook(os.path.join(self.reportdir, 'LD_block.xlsx'))
        formater = Formater(workbook)
        sheet = workbook.add_worksheet('LD_block')
        header = "Gene,L1,L2,D',LOD,r^2,CIlow,CIhi,Dist,T-int".split(',')
        sheet.set_row(0, 20)
        row = 0
        for i, j in enumerate(header):
            sheet.write(row, i, j, formater.header)
        row += 1
        for gene in self.genes:
            fblock = os.path.join(self.reportdir, 'haploview', gene + '.LD')
            try:
                count = 0
                with open(fblock, 'rt') as fh:
                    for line in fh:
                        count += 1
                        if count == 1:
                            continue
                        arr = line.split()
                        sheet.write(row, 0, gene, formater.normal)
                        for i, r in enumerate(arr):
                            sheet.write(row, i+1, r, formater.normal)
                        row += 1
            except Exception:
                print("file %s did not exist." % fblock)


class HapAssocAnalysis:
    def __init__(self, assoc_inst):
        self.config = assoc_inst.config
        self.plink = '/home/wuj/.local/bin/plink'
        self.path = self.config.get('ROUTINE', None)
        self.basepath = self.config.get('basepath')
        self.resultdir = os.path.join(self.path, 'result/haplotype')
        self.reportdir = os.path.join(self.path, 'report')
        dir_check(self.resultdir)
        self.info_file = self.config.get('INFOFILE', None)
        self.pedfile = self.config.get('PEDFILE', None) or \
                os.path.join(self.path, 'tmp/sample.ped')
        self.cov_num = self.config.get('HAP_COV', None) or self.config.get('CORRECTION', None)
        self.bedfile = self.config.get('BED', None) or \
                os.path.join(self.path, 'tmp/bsample')
        self.snp_file = self.config.get('SNPFILE', None)
        self.chrinfo = self.load_chrinfo()
        self.sampleshaps = pd.DataFrame()
        self.result_wrapper = []

        if self.config.get('HAPFILE', None) is not None and\
                file_check(self.config.get('HAPFILE')):
            self.hapfile = self.config.get('HAPFILE')
        else:
            block_identi = BlockIdentifier(assoc_inst)
            self.hapfile = block_identi.block_go()

        if self.cov_num is not None:
            self.cov_num = parse_column(self.cov_num)
            self.covar_result_wrapper = []

    def hap_go(self):
        phasebase = self.hap_phase()
        freqfile = self.hap_freq()

        sample_pheno = self.load_sampleinfo()

        hap_for_handle = self.filter_low_freq_hap(freqfile)
        self.block_sites = self.record_hapfile()

        for block in hap_for_handle:
            res_wrapper = HapResultWrapper(block)
            cov_res_wrapper = HapResultWrapper(block)
            phasefile = phasebase + block
            phase_tab = self.load_phase(block, phasefile)
            for hap in hap_for_handle[block]:
                ready_table, freq = self.table_for_LR(phase_tab, sample_pheno, hap)
                result = self.perform_LR(block, hap, ready_table)
                result.update(freq)
                res_wrapper[hap] = result
                self.result_wrapper.append(res_wrapper)
                if self.cov_num:
                    result = self.perform_LR(block, hap, ready_table, covar=True)
                    result.update(freq)
                    cov_res_wrapper[hap] = result
                    self.covar_result_wrapper.append(cov_res_wrapper)
        self.put_to_excel()

    def put_to_excel(self):
        workbook = xlsxwriter.Workbook(os.path.join(self.reportdir, 'haplotype.xlsx'))
        formater = Formater(workbook)
        sheet = workbook.add_worksheet('单倍型分析')
        sheet.set_row(0, 30)
        sheet_hap = workbook.add_worksheet('haplotype')
        sheet_hap.set_row(0, 30)
        self.sampleshap_to_excel(sheet_hap, formater)
        sheet_readme = workbook.add_worksheet('ReadMe')
        readmefile = os.path.join(self.basepath, 'ReadMetxt/readme_hap.txt')
        print_readme(sheet_readme, readmefile, formater)

        header = ('Hap', 'CHR', 'SNPS', 'HAPLOTYPE', 'case_F', 'control_F', 'OR', '95%CI', 'P-value')
        row = 0
        for i, j in enumerate(header):
            sheet.write(row, i, j, formater.header)
        row += 1
        for result in self.result_wrapper:
            block = result.block
            snps = self.block_sites[block]
            Chr= self.chrinfo[snps[0]]
            for hap in result:
                line = self.parse_result(block, Chr, ','.join(snps), hap, result[hap])
                fmt = formater_type(line, [8], formater)
                for n, v in enumerate(line):
                    sheet.write(row, n, v, fmt[n])
                row += 1

        if self.cov_num and self.covar_result_wrapper:
            workbook = xlsxwriter.Workbook(os.path.join(self.reportdir, 'haplotype_correction.xlsx'))
            formater = Formater(workbook)
            sheet = workbook.add_worksheet('单倍型分析')
            sheet.set_row(0, 30)
            sheet_hap = workbook.add_worksheet('haplotype')
            sheet_hap.set_row(0, 30)
            self.sampleshap_to_excel(sheet_hap, formater)
            sheet_readme = workbook.add_worksheet('ReadMe')
            readmefile = os.path.join(self.basepath, 'ReadMetxt/readme_hap.txt')
            print_readme(sheet_readme, readmefile, formater)

            row = 0
            for i, j in enumerate(header):
                sheet.write(row, i, j, formater.header)
            row += 1
            for result in self.covar_result_wrapper:
                block = result.block
                snps = self.block_sites[block]
                Chr= self.chrinfo[snps[0]]
                for hap in result:
                    line = self.parse_result(block, Chr, ','.join(snps), hap, result[hap])
                    fmt = formater_type(line, [8], formater)
                    for n, v in enumerate(line):
                        sheet.write(row, n, v, fmt[n])
                    row += 1
        workbook.close()

    def sampleshap_to_excel(self, sheet, formater):
        header = list(self.sampleshaps.columns)
        header.insert(0, 'Sample')
        row = 0
        for i, v in enumerate(header):
            sheet.write(row, i, v, formater.header)
        row += 1
        for r in self.sampleshaps.index:
            sheet.write(row, 0, r, formater.normal)
            line = self.sampleshaps.loc[r]
            for n, v in enumerate(line):
                sheet.write(row, n+1, v, formater.normal)
            row += 1

    def hap_phase(self):
        output = os.path.join(self.resultdir, 'phase')
        subprocess.run([self.plink,
                        '--bfile', self.bedfile,
                        '--hap', self.hapfile,
                        '--hap-phase', '--allow-no-sex',
                        '--out', output,
                        '--noweb'])
        return output + '.phase-'

    def hap_freq(self):
        output = os.path.join(self.resultdir, 'freq')
        subprocess.run([self.plink,
                        '--bfile', self.bedfile,
                        '--hap', self.hapfile,
                        '--hap-freq', '--allow-no-sex',
                        '--out', output,
                        '--noweb'])
        return output + '.frq.hap'

    def load_chrinfo(self):
        dic = {}
        count = 0
        with open(self.snp_file, 'rt') as fh:
            for line in fh:
                count += 1
                if count == 1:
                    continue
                arr = line.split()
                dic[arr[0]] = arr[1]
        return dic

    def filter_low_freq_hap(self, freqfile):
        """Keep blocks that need to be handled, drop those with
        very low frequency.
        """
        dic = {}
        count = 0
        with open(freqfile, 'rt') as fh:
            for line in fh:
                count += 1
                if count == 1:
                    continue
                arr = line.split()
                dic.setdefault(arr[0], []).append(arr[1])
        return dic

    def record_hapfile(self):
        dic = {}
        with open(self.hapfile,'rt') as fh:
            for line in fh:
                arr = line.strip().split()
                dic[arr[1]] = arr[2:]
        return dic

    def load_sampleinfo(self):
        info_table = pd.read_table(self.info_file, header=0, index_col=0, sep='\t')
        info_table.iloc[:, 0].replace('case', 1, inplace=True)
        info_table.iloc[:, 0].replace('control', 0, inplace=True)
        self.sampleshaps.reindex(info_table.index)

        # handling correction
        if self.cov_num:
            header = info_table.columns
            cols = deepcopy(self.cov_num)
            cols.insert(0,0)
            tmp_tb = info_table.iloc[:, cols]
            head = list(tmp_tb.columns)
            head[0] = 'grp'
            tmp_tb.columns = head
            return tmp_tb

        smaple_pheno = info_table.iloc[:, 0]
        sample_pheno.columns = ['grp']
        return sample_pheno

    def load_phase(self, block, filename):
        """read phase file generated by plink --hap command and filter out
        meaningless records by PH and BEST values.
        """
        table = pd.read_table(filename, header=0, index_col=0, sep='\s+')
        nonas = table[table['PH'].notnull()]
        nonas = nonas[nonas['BEST'] != 0]
        self.sampleshaps[block] = nonas.apply(self.sample_haplotype, axis=1)
        return nonas

    def table_for_LR(self, phase_table, pheno_table, hap):
        """prepare dataset for LR analysis for each haplotype."""
        merged = pd.concat([pheno_table, phase_table], axis=1, join_axes=[pheno_table.index])
        # `ind_var` is actually current haplotype
        merged['ind_var'] = merged.apply(self.hap_numeralization, args=(hap,), axis=1)
        instance_count = merged.groupby(['grp', 'ind_var']).size()
        frq_case = self.cal_freq(instance_count[1])
        frq_ctl = self.cal_freq(instance_count[0])
        freq = dict(fcase=frq_case, fctl=frq_ctl)
        return merged, freq

    @staticmethod
    def cal_freq(series):
        fmt = '{0}({1})'
        indexes = series.index
        total = series.sum()
        if 2 in indexes:
            count = series[1] + series[2] * 2
            freq = count/(total * 2)
        else:
            try:
                count = series[1]
                freq = count/(total * 2)
            except:
                count = 0
                freq = 0
        return fmt.format(count, freq)

    def perform_LR(self, block, hap, merged_table, covar=False):
        lr_res_dir = os.path.join(self.resultdir, 'LR_result')
        dir_check(lr_res_dir)
        if self.cov_num and covar:
            output = os.path.join(lr_res_dir, '.'.join([block, hap, 'cov.result']))
            header = merged_table.columns
            col_names = [header[i] for i in self.cov_num]
            col_names.append('ind_var')
            formula = 'grp ~ ' + '+'.join(col_names)
            y, X = patsy.dmatrices(formula, merged_table)
            logit = LogitRegression(y, X)
            result = logit.gofit()
            with open(output, 'wt') as fh:
                fh.write(str(result.summary2()))
        else:
            output = os.path.join(lr_res_dir, '.'.join([block, hap, 'result']))
            y, X = patsy.dmatrices('grp ~ ind_var', merged_table)
            logit = LogitRegression(y, X)
            result = logit.gofit()
            with open(output, 'wt') as fh:
                fh.write(str(result.summary2()))
        # OR and 95%CI calculation
        conf = pd.DataFrame(result.conf_int())
        conf['OR'] = result.params
        conf.columns = ['L95', 'U95', 'OR']
        res = np.exp(conf)
        res['pvalue'] = result.pvalues
        # last row is 'ind_var', i.e hap variable
        return dict(res.iloc[-1])

    @staticmethod
    def hap_numeralization(series, hap):
        return list(series).count(hap)

    @staticmethod
    def sample_haplotype(series):
        return '/'.join(series[2:4])

    @staticmethod
    def parse_result(block, Chr, SNPs, hap, dic=None):
        line = [block, Chr, SNPs, hap]
        line.append(dic.get('fcase', ''))
        line.append(dic.get('fctl', ''))
        line.append(dic.get('OR', ''))
        line.append('-'.join([str(dic.get('L95', '')), str(dic.get('U95', ''))]))
        line.append(dic.get('pvalue', ''))
        return line


class HapResultWrapper(UserDict):
    """A class to handle analysis result for each block."""
    def __init__(self, block=None, hap=None, result=None):
        self.block = block
        self.data = {}

        if hap is not None and result is not None:
            self.data[hap] = result

    def __setitem__(self, key, value):
        self.data[str(key)] = value

    def __getitem__(self, key):
        try:
            return self.data[str(key)]
        except KeyError:
            return None













