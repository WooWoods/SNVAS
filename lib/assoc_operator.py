"""
    plink_operator module
    ~~~~~~~~~~~~~~~~~~~~~

    Implements all kinds of plink operations.
"""

import os
import re
import subprocess
from functools import wraps

from .utils import dir_check
from .assoc import default_config
plink = default_config.get('PLINK')


def plink_operator(filetype, *args):
    def decorator(func):
        @wraps(func)
        def wrapper(*opts):
            filename, outname, *rest = func(*opts)
            commands = [plink, filetype, filename, '--out', outname, '--allow-no-sex'] + list(args)
            subprocess.run(commands)

        return wrapper
    return decorator

def genetic_models(*args):
    def decorator(func):
        @wraps(func)
        def wrapper(*opts):
            analysis, *models = args
            options = func(*opts)
            filename = options.filename
            outname = options.outname
            basedir = os.path.dirname(outname)
            covar = options.covar
            pheno = options.pheno
            for model in models:
                tmpname = outname + model
                if model.strip():
                    commands = [plink, '--bfile', filename, analysis, model,
                            '--adjust', '--ci', '0.95',
                            '--out', tmpname, '--allow-no-sex']
                else:
                    commands = [plink, '--bfile', filename, analysis,
                            '--adjust', '--ci', '0.95',
                            '--out', tmpname, '--allow-no-sex']
                subprocess.run(commands)

                if covar is not None:
                    outdir = os.path.join(basedir, 'logit_covar')
                    dir_check(outdir)
                    tmpname = os.path.join(outdir, 'logistic%s' %model)
                    if model.strip():
                        commands = [plink, '--bfile', filename, analysis, model,
                                '--adjust', '--ci', '0.95', '--covar', covar,
                                '--out', tmpname, '--allow-no-sex']
                    else:
                        commands = [plink, '--bfile', filename, analysis,
                                '--adjust', '--ci', '0.95', '--covar', covar,
                                '--out', tmpname, '--allow-no-sex']
                    subprocess.run(commands)
                if pheno is not None:
                    outdir = os.path.join(basedir, 'phenoassoc')
                    dir_check(outdir)
                    tmpname = os.path.join(outdir, 'logit_%s' %model)
                    if model.strip():
                        commands = [plink, '--bfile', filename, analysis, model,
                                '--adjust', '--ci', '0.95', '--pheno', pheno,
                                '--all-pheno', '--out', tmpname, '--allow-no-sex']
                    else:
                        commands = [plink, '--bfile', filename, analysis,
                                '--adjust', '--ci', '0.95', '--pheno', pheno,
                                '--all-pheno', '--out', tmpname, '--allow-no-sex']
                    subprocess.run(commands)
                if pheno is not None and covar is not None:
                    outdir = os.path.join(basedir, 'phenoassoc_covar')
                    dir_check(outdir)
                    tmpname = os.path.join(outdir, 'logit_%s' %model)
                    if model.strip():
                        commands = [plink, '--bfile', filename, analysis, model,
                                '--adjust', '--ci', '0.95', '--pheno', pheno,
                                '--all-pheno', '--covar', covar,
                                '--out', tmpname, '--allow-no-sex']
                    else:
                        commands = [plink, '--bfile', filename, analysis,
                                '--adjust', '--ci', '0.95', '--pheno', pheno,
                                '--all-pheno', '--covar', covar,
                                '--out', tmpname, '--allow-no-sex']
                    subprocess.run(commands)
        return wrapper
    return decorator

