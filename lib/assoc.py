"""
    assoc module
    ~~~~~~~~~~~~

    Implements the central functions.
"""

import os

from .config import Config


class AssocStudy:
    """An object that implements association studies."""

    config_class = Config

    default_config = {
            'PLINK':  '/home/wuj/bin/software/plink_1.90_beta/plink',

            }

    def __init__(self,):

        self.config = make_config()

    def make_config(self):
        """Used to create the config attribute."""

        return self.config_class(self.default_config)


