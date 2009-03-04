"""diff checker for SOG buildbot.

diff the specified files and check the result.  If the files differ
only in the content of the RunDateTime or Date: header line, exit with
status 0 and no output.  Otherwise, exit with status 1 and the diff as
output.

test_check_diff.py is the test suite for this script.

:Author: Doug Latornell <dlatorne@eos.ubc.ca>
:Created: 2009-03-01
"""

import sys
from subprocess import Popen, PIPE, STDOUT


def check_diff(argv):
    """diff checker for SOG buildbot.

    diff the specified files and check the result.  If the files
    differ only in the content of the RunDateTime or Date: header
    line, exit with status 0 and no output.  Otherwise, exit with
    status 1 and the diff as output.
    """
    if len(argv) != 3:
        return_code = 1
        print 'usage: %(prog)s file1 file2' % {'prog':argv[0]}
        raise SystemExit(return_code)
    # diff the files
    proc = Popen(['diff', argv[1], argv[2]],
                 stderr=PIPE, stdout=PIPE)
    proc.wait()
    diff, stderr = proc.communicate()
    split_diff = diff.split('\n')
    # Check the result
    if stderr.strip():
        return_code = 1
        diff = stderr
    elif len(split_diff) == 1:
        return_code = 0
    elif len(split_diff) > 5:
        return_code = 1
    else:
        line1 = split_diff[1][1:].strip()
        line3 = split_diff[3][1:].strip()
        if ((line1.startswith('*RunDateTime') 
             and line3.startswith('*RunDateTime'))
            or (line1.startswith('Date =') and line3.startswith('Date ='))):
            return_code = 0
        else:
            return_code = 1
    # Return the status, and the diff if it was dirty
    if return_code == 1:
        print diff
    raise SystemExit(return_code)


if __name__== '__main__':
    check_diff(sys.argv)
