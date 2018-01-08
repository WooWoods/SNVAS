"""
    stratify module
    ~~~~~~~~~~~~~~~

    Implements stratification.
"""

import os
from shutil import copyfile
from itertools import combinations

import pandas as pd
from .utils import dir_check, parse_column


class Stratify:
    """Class for stratification, it just prepares the files that needed by
    main program according to information provided by user. After this
    process, user may get a series of sub-projects which may meet expectatioin
    of user, and then the user need to check these sub-projects and run main
    program mannuly.

    :param assoc_inst: an instance of AssocStudy
    """
    def __init__(self, assoc_inst):
        self.config = assoc_inst.config
        self.root_path = self.config.get('ROUTINE', '')
        self.geno_file = self.config.get('GENOFILE', '')
        self.info_file= self.config.get('INFOFILE', '')
        self.snp_file  = self.config.get('SNPFILE', '')
        self.gender_col= self.config.get('GENDER', '')
        self.cov_num = self.config.get('CORRECTION', '')
        self.pheno_num = self.config.get('PHENO', '')
        self.stratify = self.config.get('STRATIFY', None)
        if self.stratify is None:
            raise Exception('Strafication information not provided.')
        self.hapfile = self.config.get('HAPFILE', '')
        self.chi_test = self.config.get('CHI_TEST', '')
        self.ttest = self.config.get('TTEST', '')

    def load_table(self):
        try:
            self.geno_tab = pd.read_table(self.geno_file, header=0, index_col=0, sep='\t')
            # geno_tab.replace('\/', ' ', regex=True, inplace=True)
            # self.geno_tab = geno_tab.fillna('0 0')
        except IOError as e:
            if e.errno in (errno.ENOENT, errno.EISDIR):
                e.strerror = 'Unable to load geno_file <%s>' % e.strerror
                raise Exception()
            raise Exception('Melformed geno_file <%s>' % e.strerror)

        try:
            info_tab = pd.read_table(self.info_file, header=0, index_col=0, sep='\t')

            if self.gender_col:
                info_tab.iloc[:, self.gender_col - 1].fillna('-9')
            self.info_tab = info_tab
            self.info_names = self.info_tab.columns
        except IOError as e:
            if e.errno in (errno.ENOENT, errno.EISDIR):
                e.strerror = 'Unable to load geno_file <%s>' % e.strerror
                raise Exception()
            raise Exception('Melformed geno_file <%s>' % e.strerror)

    def go(self):
        self.load_table()
        for group in self.stratify:
            combinates, col, n = self.strati_groups(group)
            for combinate in combinates:
                self.sample_by_combinate(combinate, col, n)

    def strati_groups(self, group):
        """decide stratification groups by cols and n provided by user.

        :param group: a tuple like (2, 1), in which 2 means cols in sample.info used
                      for stratification, 1 means n for combinations.
        """
        col, n = group
        col -= 1
        col_name = self.info_names[col]
        col_values = set(self.info_tab.iloc[:, col][self.info_tab.iloc[:, col].notnull()].values)
        combinates = list(combinations(col_values, n))
        if n > 1:
            if len(col_values) < 2:
                raise Exception('%s variable contains values %s less than 2, not enough for grouping' %(col_name, col_values))
        return combinates, col, n

    def sample_by_combinate(self, combinate, strati_col, n):
        """Extracting samples from sample.info according to combinate.

        :param combinate: combination of values of strati_col.
        :param strati_col: viaing which column to stratifying.
                e.g.
                  stratify by strati_col=2, which comes to be gender,
                  and n=1, that means stratifying by male or female,
                  combinate could be
                  (1, ) or (2, )
        """
        col_name = self.info_names[strati_col]
        if n == 1:
            dirname = col_name + '_' + str(combinate[0])
        else:
            dirname = col_name + '-VS-'.join(map(lambda x: str(x),combinate))
        strati_root = os.path.join(self.root_path, dirname)
        dir_check(strati_root)
        data_path = os.path.join(strati_root, 'data')
        dir_check(data_path)
        snpfile = os.path.join(data_path, 'anno.txt')
        copyfile(self.snp_file, snpfile)
        strati_infofile = os.path.join(data_path, 'sample.info')
        strati_genofile = os.path.join(data_path, 'sample.geno')
        if n == 1:
            tmp_info = self.info_tab[self.info_tab.iloc[:, strati_col].isin(combinate)]
        else:
            partition_info = self.info_tab[self.info_tab.iloc[:, strati_col].isin(combinate)]
            tmp_info = pd.concat([partition_info.iloc[:, strati_col], partition_info.iloc[:, 1: strati_col], partition_info.iloc[:, (strati_col + 1):]], axis=1)
            combinate = list(combinate)
            combinate.sort(reverse=True)

            tmp_info.iloc[:, 0].replace(combinate[0], 'case', inplace=True)
            tmp_info.iloc[:, 0].replace(combinate[1], 'control', inplace=True)
            # tmp_info.iloc[:, 0][tmp_info.iloc[:, strati_col] == combinate[0]] = 'case'
            # tmp_info.iloc[:, 0][tmp_info.iloc[:, strati_col] == combinate[1]] = 'control'
        samples = tmp_info.index
        tmp_geno = self.geno_tab.loc[samples]
        tmp_info.apply(self.convert_dtype).to_csv(strati_infofile, header=True, index=True, sep='\t')
        tmp_geno.to_csv(strati_genofile, header=True, index=True, sep='\t')
        self.print_configfile(strati_genofile, strati_infofile, snpfile, strati_root)

    def print_configfile(self, fgeno, finfo, fsnp, root_path):
        config_model ="""ROUTINE = '{rootpath}'
GENOFILE = '{genofile}'
INFOFILE = '{infofile}'
SNPFILE = '{snpfile}'
HAPFILE = '{hapfile}'
GENDER = '{gender}'
CORRECTION = {covar}
PHENO = {pheno}
CHI_TEST = {chi}
TTEST = {ttest}
        """
        tmpdict = dict(rootpath=root_path,
                genofile=fgeno,
                infofile=finfo,
                snpfile=fsnp,
                hapfile='',
                gender=self.gender_col,
                covar=self.cov_num,
                pheno=self.pheno_num,
                chi=self.chi_test,
                ttest=self.ttest)
        filename = os.path.join(root_path, 'config.ini')
        with open(filename, 'wt') as fh:
            fh.write(config_model.format_map(tmpdict))

    @staticmethod
    def convert_dtype(x):
        try:
            return x.astype(int)
        except:
            return x




