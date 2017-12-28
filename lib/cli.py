"""
    cli module
    ~~~~~~~~~~

    Command line options.
"""

import os
import argparse

from . import AssocStudy, Formater, MdrOperate, hap_analysis, reporter, PhenoIndepTest,\
        Stratify, LRanalysis, Chi_test


AP = argparse.ArgumentParser(
        description="ASkit, one-stop solution for case-control association analysis.",
        formatter_class=argparse.RawTextHelpFormatter,
        )

AP_subparsers = AP.add_subparsers(
        help="Sub-commands (use with -h for more info)"
        )

##########################################################################
### Batch
##########################################################################

def _batch_command(args):
    """Perform complete analysis including:
        - hwe analysis
        - chi-test of genotype and case/control
        - logit-test of genotype and case/control
        - chi-test and t-test of phenotype and case/control
        - association of phenotype and genotype
        - mdr anasysis for gene X gene interaction
        - haplotype analysis.
    Usage:
        ASkit.py batch -cfg config.ini
    """
    curr_case = AssocStudy(args.cfg)
    formater = Formater(curr_case)
    formater.make_bed()
    curr_case.batch_run()
    reporter(curr_case)

    mdr = MdrOperate(curr_case)
    mdr.go()
    hap_analysis(curr_case)


P_batch = AP_subparsers.add_parser('batch', help=_batch_command.__doc__)
P_batch.add_argument('-cfg', metavar='config file', required=True)
P_batch.set_defaults(func=_batch_command)


##########################################################################
### Plink_analysis
##########################################################################

def _plink_stage(args):
    """Perform plink association analysis including:
        - hwe analysis
        - plink association analysis
    Usage:
        ASkit.py plink -cfg config.ini
    """
    curr_case = AssocStudy(args.cfg)
    formater = Formater(curr_case)
    formater.make_bed()
    curr_case.batch_run()
    reporter(curr_case)

P_plink = AP_subparsers.add_parser('plink', help=_plink_stage.__doc__)
P_plink.add_argument('-cfg', metavar='config file',required=True)
P_plink.set_defaults(func=_plink_stage)

##########################################################################
### Plink_analysis
##########################################################################

def _plink_report(args):
    """Perform plink association analysis result report.
    Usage:
        ASkit.py report -cfg config.ini
    """
    curr_case = AssocStudy(args.cfg)
    reporter(curr_case)

P_report = AP_subparsers.add_parser('report', help=_plink_report.__doc__)
P_report.add_argument('-cfg', metavar='config file',required=True)
P_report.set_defaults(func=_plink_report)

##########################################################################
### MDR
##########################################################################

def _mdr_stage(args):
    """Perform mdr analysis.
    Usage:
        ASkit.py mdr -cfg config.ini
    """
    curr_case = AssocStudy(args.cfg)
    formater = Formater(curr_case)
    formater.make_bed()

    MdrOperate(curr_case)
    MdrOperate.go()

P_mdr = AP_subparsers.add_parser('mdr', help=_mdr_stage.__doc__)
P_mdr.add_argument('-cfg', metavar='config file',required=True)
P_mdr.set_defaults(func=_mdr_stage)


##########################################################################
### Haplotype
##########################################################################

def _hap_stage(args):
    """Perform haplotype analysis.
    Usage:
        ASkit.py hap -cfg config.ini
    """
    curr_case = AssocStudy(args.cfg)
    formater = Formater(curr_case)
    formater.make_bed()
    hap_analysis(curr_case)

P_hap = AP_subparsers.add_parser('hap', help=_hap_stage.__doc__)
P_hap.add_argument('-cfg', metavar='config file',required=True)
P_hap.set_defaults(func=_hap_stage)


##########################################################################
### PhenoTest
##########################################################################

def _pheno_stage(args):
    """Perform chi-test and t-test for phenotype-VS-case/control.
    Usage:
        ASkit.py pheno -cfg config.ini
    """
    curr_case = AssocStudy(args.cfg)
    formater = Formater(curr_case)
    pheno_test = PhenoIndepTest(curr_case)
    pheno_test.go()

P_pheno = AP_subparsers.add_parser('pheno', help=_pheno_stage.__doc__)
P_pheno.add_argument('-cfg', metavar='config file',required=True)
P_pheno.set_defaults(func=_pheno_stage)

##########################################################################
### PhenoTest
##########################################################################

def _stratify(args):
    """Preparing stratification sub-project by `STRATIFY` information
    provided in the config.ini file.
    Usage:
        ASkit.py strati -cfg config.ini
    """
    curr_case = AssocStudy(args.cfg)
    stratification = Stratify(curr_case)
    stratification.go()

P_strati = AP_subparsers.add_parser('strati', help=_stratify.__doc__)
P_strati.add_argument('-cfg', metavar='config file',required=True)
P_strati.set_defaults(func=_stratify)

##########################################################################
### Logistic Regressioin
##########################################################################

def _LR(args):
    """Logistic regression use files and R-style formula.
    :param filename: a <TAB> delimited text file.
    :param formula: R-style formula. e.g. y ~ a + b,
                    in which 'y', 'a', 'b' are column names of the file
    Usage:
        ASkit.py LR -filename /home/example.txt -formula y~a+b
    """
    filename = args.filename
    formula = args.formula
    result = LRanalysis(filename, formula)
    print(result)

P_LR = AP_subparsers.add_parser('LR', help=_LR.__doc__)
P_LR.add_argument('-filename', metavar='a <TAB> delimited text file to be analysised',required=True)
P_LR.add_argument('-formula', metavar='a formula like `y ~ a + b`',required=True)
P_LR.set_defaults(func=_LR)

##########################################################################
### Chisquare test
##########################################################################

def _chi_test(args):
    """2X2 contingency chi-square analysis and OR, 95 percent CI calculations.
    :param filename: a <TAB> delimited text file.
    :param y: variable name for groupby.
    :param x: another variable.
    Usage:
        ASkit.py -filename /home/example.txt -y group -x gender
    """
    filename = args.filename
    y = args.y
    x = args.x
    result = Chi_test(filename, y, x)
    print(result)

P_chi = AP_subparsers.add_parser('chi', help=_chi_test.__doc__)
P_chi.add_argument('-filename', metavar='a <TAB> delimited text file to be analysised',required=True)
P_chi.add_argument('-y', metavar='variable',required=True)
P_chi.add_argument('-x', metavar='another variable',required=True)
P_chi.set_defaults(func=_chi_test)

##########################################################################
### Shim for command-line execution
##########################################################################

def parse_args(args=None):
    """Parse the command line."""
    return AP.parse_args(args=args)

