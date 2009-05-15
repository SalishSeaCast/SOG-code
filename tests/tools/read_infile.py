"""infile reader for SOG buildbot.

test_infile_reader.py is the test suite for this script.

:Author: Doug Latornell <dlatorne@eos.ubc.ca>
:Created: 2009-03-07
"""
import re
import sys


def read_infile(argv):
    """infile reader for SOG buildbot.

    Read the specified SOG infile and return the output file names
    separated by newlines.
    """
    if len(argv) != 3:
        return_code = 1
        print 'usage: %(prog)s infile key' % {'prog':argv[0]}
        raise SystemExit(return_code)
    key = ' '.join(argv[2:])
    p = re.compile(r'"\s')
    try:
        infile = open(argv[1]).readlines()
    except IOError, msg:
        return_code = 1
        print msg
        raise SystemExit(return_code)
    for i, line in enumerate(infile):
        # Skip empty lines and comments
        if line.startswith('\n') or line.startswith('!'):
            continue
        # Split line into [key, value, comment]
        split_line = p.split(line)
        if split_line[0].strip('"')  == key:
            if split_line[1].strip('"').strip().rstrip('\n'):
                # Value on the same line as key
                print split_line[1].strip('"')
            else:
                # Value on line after key
                print p.split(infile[i+1])[0].strip().strip('"')
        return_code = 0
    raise SystemExit(return_code)


if __name__== '__main__':
    read_infile(sys.argv)
