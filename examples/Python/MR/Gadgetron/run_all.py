#!/usr/bin/env python
import glob
import os
import sys

for i in glob.glob('*.py'):
    if os.path.abspath(__file__) == os.path.abspath(i):
        continue
    print(i)
    os.system(sys.executable + ' ' + i)
