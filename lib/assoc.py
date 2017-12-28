"""
    assoc module
    ~~~~~~~~~~~~

    Implements the central functions of association study.
"""

import os
import re
import subprocess
from collections import namedtuple

from .config import Config
from .utils import dir_check, file_check


default_config = {
        'PLINK':  '/home/wuj/bin/software/plink_1.90_beta/plink',
        'annovar': '/home/pub/software/annovar/annotate_variation.pl',
        'basepath': os.path.abspath(os.path.dirname(__file__)),
        }


class AssocStudy:
    """An object that implements association studies."""
    from .assoc_operator import plink_operator, genetic_models

    config_class = Config

    def __init__(self, cfgfile=None):
        self.config = self.make_config()

        if os.path.isfile(cfgfile):
            self.config.from_pyfile(cfgfile)

    def make_config(self):
        """Used to create the config attribute."""
        return self.config_class(default_config)

    def batch_run(self):
        self.hwe()
        self.freq()
        self.freqcc()
        self.annotation()
        self.chitest()
        self.modelchi()
        self.logistic()

    def annotation(self):
        annovar = self.config.get('annovar')
        snpfile = self.config.get('SNPFILE', None)
        outdir = os.path.join(self.config.get('ROUTINE'), 'result/hwe')
        dir_check(outdir)
        library = self.library_prepare(snpfile, outdir)
        subprocess.run([annovar, '--hgvs', '-filter',
                        '-dbtype', '1000g2014oct_chbs',
                        '--buildver', 'hg19',
                        library, '/home/pub/database/Human/hg19/Annotation'])
        subprocess.run([annovar, '--hgvs', '--splicing_threshold',
                        '8', '--buildver', 'hg19', library,
                        '/home/pub/database/Human/hg19/Annotation'])

    @plink_operator('--bfile', '--hardy')
    def hwe(self):
        filename = self.config.get('BED', None)
        outdir = os.path.join(self.config.get('ROUTINE'), 'result/hwe')
        dir_check(outdir)
        outname = os.path.join(outdir, 'hwe')
        Options = namedtuple('Opts', 'filename outname covar pheno')

        options = Options(filename=filename, outname=outname, covar=None, pheno=None)

        return options

    @plink_operator('--bfile', '--freq')
    def freq(self):
        filename = self.config.get('BED', None)
        outdir = os.path.join(self.config.get('ROUTINE'), 'result/hwe')
        dir_check(outdir)
        outname = os.path.join(outdir, 'freq')
        Options = namedtuple('Opts', 'filename outname covar pheno')

        options = Options(filename=filename, outname=outname, covar=None, pheno=None)

        return options

    @plink_operator('--bfile', '--freq', 'case-control')
    def freqcc(self):
        filename = self.config.get('BED', None)
        outdir = os.path.join(self.config.get('ROUTINE'), 'result/hwe')
        dir_check(outdir)
        outname = os.path.join(outdir, 'freq')
        Options = namedtuple('Opts', 'filename outname covar pheno')

        options = Options(filename=filename, outname=outname, covar=None, pheno=None)

        return options

    @plink_operator('--bfile', '--assoc', '--adjust', '--ci', '0.95')
    def chitest(self):
        filename = self.config.get('BED', None)
        outdir = os.path.join(self.config.get('ROUTINE'), 'result/chi-test')
        dir_check(outdir)
        outname = os.path.join(outdir, 'chi')

        Options = namedtuple('Opts', 'filename outname covar pheno')
        options = Options(filename=filename, outname=outname, covar=None, pheno=None)
        return options

    @plink_operator('--bfile', '--assoc', 'fisher', '--adjust', '--ci', '0.95')
    def fishertest(self):
        filename = self.config.get('BED', None)
        outdir = os.path.join(self.config.get('ROUTINE'), 'result/fisher-test')
        dir_check(outdir)
        outname = os.path.join(outdir, 'fisher')

        Options = namedtuple('Opts', 'filename outname covar pheno')
        options = Options(filename=filename, outname=outname, covar=None, pheno=None)
        return options

    @plink_operator('--bfile', '--model', '--cell', '0')
    def modelchi(self):
        filename = self.config.get('BED', None)
        outdir = os.path.join(self.config.get('ROUTINE'), 'result/chi-test')
        dir_check(outdir)
        outname = os.path.join(outdir, 'model_chi')

        Options = namedtuple('Opts', 'filename outname covar pheno')
        options = Options(filename=filename, outname=outname, covar=None, pheno=None)
        return options

    @plink_operator('--bfile', '--model', 'fisher')
    def modelfisher(self):
        filename = self.config.get('BED', None)
        outdir = os.path.join(self.config.get('ROUTINE'), 'result/fisher-test')
        dir_check(outdir)
        outname = os.path.join(outdir, 'model_fisher')

        Options = namedtuple('Opts', 'filename outname covar pheno')
        options = Options(filename=filename, outname=outname, covar=None, pheno=None)
        return options

    @genetic_models('--logistic', '', 'dominant', 'recessive', 'hethom')
    def logistic(self):
        filename = self.config.get('BED', None)
        outdir = os.path.join(self.config.get('ROUTINE'), 'result/logistic-test')
        dir_check(outdir)
        outname = os.path.join(outdir, 'logistic')

        covar_file = self.config.get('COVARFILE', None)
        pheno_file = self.config.get('PHENOFILE', None)

        covar = None
        pheno = None

        if covar_file is not None and file_check(covar_file):
            covar = covar_file
        if pheno_file is not None and file_check(pheno_file):
            pheno = pheno_file

        Options = namedtuple('Opts', 'filename outname covar pheno')
        options = Options(filename=filename, outname=outname, covar=covar, pheno=pheno)
        return options

    @staticmethod
    def library_prepare(filename, outdir):
        output = os.path.join(outdir, 'library')
        count = 0
        with open(filename, 'rt') as fh,\
                open(output, 'wt') as foh:
            for line in fh:
                count += 1
                if re.match(r'^\s$', line):
                    continue
                if count == 1:
                    continue
                arr = line.strip().split()
                snp, Chr, pos, *rest = arr
                if re.search('-', pos):
                    start, end = pos.split('-')
                else:
                    start = pos
                    end = pos
                newline = '\t'.join([Chr, start, end] + rest + [snp])
                foh.write(newline + '\n')
        return output









