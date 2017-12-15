#!/usr/bin/env python3

import os
import sys
from lib import cli


if __name__ == '__main__':
    args = cli.parse_args()
    args.func(args)


