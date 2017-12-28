"""
    _formation module
    ~~~~~~~~~~~~~~~~~

    This module works to normalize data input for plink.
"""

import os
import re
import errno
import functools

import pandas as pd

from .utils import dir_check, parse_column


class Formater:
    """Take in tab separated genotype file and phenotype file, and turns
    them into plink format.

    :param config: a config instance.
    """
    from .assoc_operator import plink_operator

    def __init__(self, asinst):
        self.config = asinst.config
        self.root_path = self.config.get('ROUTINE', None)
        self.geno_file = self.config.get('GENOFILE', None)
        self.info_file= self.config.get('INFOFILE', None)
        self.snp_file  = self.config.get('SNPFILE', None)
        self.gender_col= self.config.get('GENDER', None)
        self.cov_num = self.config.get('CORRECTION', None)
        self.pheno_num = self.config.get('PHENO', None)
        self.hap_cutofff = self.config.get('HAP_CUTOFF', None) or 100
        self.mdr_analysis = self.config.get('MDR') or False
        self.hap_treat = None

        self.tmpdir = os.path.join(self.root_path, 'tmp')
        dir_check(self.tmpdir)
        self.config['TMPDIR'] = self.tmpdir

        self.load_table()
        self.to_ped()
        self.to_map()

    def load_table(self):
        """Load genotype file and phenotype file into pandas DataFrame. And check if these files
        exists or matched."""
        try:
            geno_tab = pd.read_table(self.geno_file, header=0, index_col=0, sep='\t', low_memory=False)
            geno_tab.replace('\s+', '0 0', regex=True, inplace=True)
            geno_tab.replace('\/', ' ', regex=True, inplace=True)
            self.geno_tab = geno_tab.fillna('0 0')
        except IOError as e:
            if e.errno in (errno.ENOENT, errno.EISDIR):
                e.strerror = 'Unable to load geno_file <%s>' % e.strerror
                raise Exception()
            raise Exception('Melformed geno_file <%s>' % e.strerror)

        try:
            info_tab = pd.read_table(self.info_file, header=0, index_col=0, sep='\t')
            self.info_tab = info_tab.fillna('-9')
        except IOError as e:
            if e.errno in (errno.ENOENT, errno.EISDIR):
                e.strerror = 'Unable to load info_file <%s>' % e.strerror
                raise Exception()
            raise Exception('Melformed info_file <%s>' % e.strerror)

        shape_geno = geno_tab.shape
        shape_info = info_tab.shape
        if not shape_geno[0] == shape_info[0]:
            raise Exception('Samples not match: <geno: %s / pheno: %s>' %(shape_geno[0], shape_info[0]))

    def to_ped(self):
        """Turn genotype file into plink ped format."""
        self.snv_sites = self.geno_tab.columns
        self.info_tab.replace('case', 2, inplace=True)
        self.info_tab.replace('control', 1, inplace=True)
        # self.info_tab.iloc[:,0][self.info_tab.iloc[:,0] == 'case'] = 2
        # self.info_tab.iloc[:,0][self.info_tab.iloc[:,0] == 'control'] = 1

        shape_info = self.info_tab.shape
        merged = pd.concat([self.info_tab, self.geno_tab], axis=1)
        #merged.to_csv('merged.csv', header=True, index=True, sep='\t')
        shape_merged = merged.shape
        require_cols = [0] + list(range(shape_info[1], shape_merged[1]))
        require_cols_mdr = list(range(shape_info[1], shape_merged[1])) + [0]

        required_tab = merged.iloc[:, require_cols]

        if self.gender_col is not None and re.match(r'\d', str(self.gender_col)):
            required_tab.insert(0, 'gender', merged.iloc[:, self.gender_col - 1])
        else:
            required_tab.insert(0, 'gender', 0)

        required_tab.insert(0, 'MID',0)
        required_tab.insert(0, 'PID',0)
        required_tab.insert(0, 'IID',required_tab.index)

        ped_file = os.path.join(self.tmpdir, 'sample.ped')
        required_tab.apply(self.convert_dtype).to_csv(ped_file, header=False, index=True, sep='\t')
        self.config['PEDFILE'] = ped_file

        self.to_covar(merged)
        self.to_pheno(merged)
        if self.mdr_analysis:
            self.to_mdr(merged, require_cols_mdr)

    def to_map(self):
        """Turn genotype file into plink map format."""
        try:
            snp_info = {}
            with open(self.snp_file, 'rt') as fh:
                count = 0
                for line in fh:
                    count += 1
                    if re.match(r'^\s+$', line):
                        continue
                    arr = line.strip().split()
                    if count == 1:
                        if not set(arr) & set(['Gene', 'Chr', 'Position', 'ref', 'alt']):
                            raise Exception('Loss file header <%s>' % self.snp_file)
                        if len(arr) > 5 or re.match(r'GENE', arr[-1], re.I):
                            self.hap_treat = True
                            continue
                    snp_info[arr[0]] = arr
            if self.hap_treat and count <= self.hap_cutofff:
                self.hap_treat = True
                self.config['TREATHAP'] = True

            map_file = os.path.join(self.tmpdir, 'sample.map')
            with open(map_file, 'wt') as fh:
                for snv in self.snv_sites:
                    record = snp_info.get(snv, None)
                    if record is None:
                        raise Exception('Loss snv info: <%s>' % snv)
                    rs, chrs, pos, *rest = record
                    pos = self.parse_pos(pos)
                    fh.write('{0}\t{1}\t0\t{2}\n'.format(chrs, rs, pos))
            if self.hap_treat:
                hapfile = os.path.join(os.path.dirname(self.snp_file), 'raw_hap.txt')
                self.to_hap(hapfile, snp_info)
                self.config['RAW_HAP'] = hapfile

        except IOError as e:
            if e.errno in (errno.ENOENT, errno.EISDIR):
                e.strerror = 'Unable to load geno_file <%s>' % e.strerror
                raise Exception()

    def to_mdr(self, table, cols):
        filename = os.path.join(self.tmpdir, 'mdrdata.txt')
        required_tab = table.iloc[:, cols]
        required_tab.to_csv(filename, header=True, index=False, sep='\t')
        self.config['MDR'] = filename

    def to_covar(self, table):
        filename = os.path.join(self.tmpdir, 'covar.txt')
        if self.cov_num is not None and re.search(r'\d', str(self.cov_num)):
            cols = parse_column(self.cov_num)
            covar = table.iloc[:, cols]
            covar.insert(0, 'IID', table.index)
            covar.insert(0, 'FID', table.index)
            covar.apply(self.convert_dtype).to_csv(filename, header=True, index=False, sep='\t')
            self.config['COVARFILE'] = filename

    def to_pheno(self, table):
        filename = os.path.join(self.tmpdir, 'pheno.txt')
        if self.pheno_num is not None and re.search(r'\d', str(self.pheno_num)):
            cols = parse_column(self.pheno_num)
            covar = table.iloc[:, cols]
            covar.insert(0, 'IID', table.index)
            covar.insert(0, 'FID', table.index)
            covar.apply(self.convert_dtype).to_csv(filename, header=True, index=False, sep='\t')
            self.config['PHENOFILE'] = filename

    @staticmethod
    def convert_dtype(x):
        try:
            return x.astype(int)
        except:
            return x

    @staticmethod
    def to_hap(filename, snp_info):
        from collections import defaultdict
        tmp = defaultdict(list)
        for i in snp_info:
            gene = snp_info[i][-1]
            if len(snp_info[i]) < 6 or not gene.strip():
                continue
            tmp[gene].append(i)
        with open(filename, 'wt') as fh:
            for gene in tmp:
                fh.write('**\t{0}\t{1}\n'.format(gene, '\t'.join(tmp[gene])))

    @staticmethod
    def parse_pos(pos):
        """Handle indel regions cause plink map format only recognise base position.
        -  e.g. 180047739-180047739 will be converted to 180047739.

        :param pos: chromosome position of a snp site.
        """
        return pos.split('-')[0]

    @plink_operator('--file', '--make-bed')
    def make_bed(self):
        filename = os.path.join(self.tmpdir, 'sample')
        outname = os.path.join(self.tmpdir, 'bsample')
        self.config['BED'] =  outname
        return filename, outname




