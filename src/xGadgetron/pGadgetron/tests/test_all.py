# -*- coding: utf-8 -*-
"""test sets
v{version}

Usage:
  test_all [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

{author}

{licence}
"""
from sirf.Utilities import __licence__, RE_PYEXT
from glob import glob
from os import path
__version__ = "0.2.0"
__author__ = "Casper da Costa-Luis"


if __name__ == "__main__":
    from docopt import docopt
    args = docopt(__doc__.format(
        version=__version__, author=__author__, licence=__licence__),
                  version=__version__)

    record = args['--record']
    verbose = args['--verbose']

    failed = 0
    ntests = 0

    for script in glob('*.py'):
        if path.abspath(__file__) == path.abspath(script):
            continue
        print('\n\n--- running %s' % script)
        test = RE_PYEXT.sub("", script)
        main = RE_PYEXT.sub(".test_main(record, verbose, throw=False)", script)
        exec('import ' + test)
        f, n = eval(main)
        if f:
            print('    %d of %d tests failed' % (f, n))
        failed += f
        ntests += n

    if failed:
        import sys
        print('%d of %d tests failed' % (failed, ntests))
        sys.exit(failed)
    print('all %d tests passed' % ntests)
