"""infile reader for SOG buildbot.

Read the specified SOG infile and return the value for the specified key.

Commandline usage: python read_infile.py infile key

test_infile_reader.py is the test suite for this script.

:Author: Doug Latornell <dlatorne@eos.ubc.ca>
:Created: 2009-03-07
"""
import logging
import optparse
import os
import re
import sys


# Minimalist logging to sys.stderr.  Use logging.debug() for debug
# printing.
logging.basicConfig(format='%(message)s', level=logging.DEBUG)


def read_infile(infile, key):
    """Read the specified SOG infile and return the value for the
    specified key.

    :arg infile: path and filename of infile to read
    :type infile: string

    :arg key: key to return value for; whitespace in key is okay
    :type key: string

    :returns: value for the specified key
    :rtype: string
    """
    try:
        infile = open(infile).readlines()
    except IOError, msg:
        return_code = 2
        print msg
        raise SystemExit(return_code)
    # Key is separated from value by "+whitespace
    p = re.compile(r'"\s')
    for i, line in enumerate(infile):
        # Skip empty lines and comments
        if line.startswith('\n') or line.startswith('!'):
            continue
        # Split line into [key, value, comment]
        split_line = p.split(line)
        if split_line[0].strip('"')  == key:
            value = split_line[1].strip('"').strip().rstrip('\n')
            if value:
                # Value on the same line as key
                print value
            else:
                # Value on line after key
                print p.split(infile[i+1])[0].strip().strip('"')
        return_code = 0
    raise SystemExit(return_code)


if __name__== '__main__':
    # Parse the commandline and call read_infile
    usage = 'Usage: %prog infile key'
    parser = optparse.OptionParser(usage, prog=os.path.basename(__file__))
    options, args = parser.parse_args(args=sys.argv[1:])
    if len(args) < 2:
        return_code = 2
        parser.print_usage()
        raise SystemExit(return_code)
    infile = args[0]
    key = ' '.join(args[1:])
    read_infile(infile, key)
