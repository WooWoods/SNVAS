"""
    cli module
    ~~~~~~~~~~

    Command line options.
"""

import os
import argparse

from . import AssocStudy, Formater, MdrOperate, hap_analysis, reporter, PhenoIndepTest, Stratify


AP = argparse.ArgumentParser(
        description="SSRPA, a toolkit for population analysis with SSR markers.",
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
### MDR
##########################################################################

def _mdr_stage(args):
    """Perform mdr analysis."""
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
    """Perform haplotype analysis."""
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
    """Perform chi-test and t-test for phenotype-VS-case/control"""
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
    """Preparing stratification project by `STRATIFY` information
    provided in the config.ini file.
    """
    curr_case = AssocStudy(args.cfg)
    stratification = Stratify(curr_case)
    stratification.go()

P_strati = AP_subparsers.add_parser('strati', help=_stratify.__doc__)
P_strati.add_argument('-cfg', metavar='config file',required=True)
P_strati.set_defaults(func=_stratify)


##########################################################################
### Shim for command-line execution
##########################################################################

def parse_args(args=None):
    """Parse the command line."""
    return AP.parse_args(args=args)

