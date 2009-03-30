"""diff checker for SOG buildbot.

diff the specified files and check the result.  If the files differ
only in the content of the *RunDateTime header line, or Date, or 1 or
more output filename lines, exit with status 0 and no output.
Otherwise, exit with status 1 and the diff as output.

test_check_diff.py is the test suite for this script.

:Author: Doug Latornell <dlatorne@eos.ubc.ca>
:Created: 2009-03-01
"""
import sys
from subprocess import Popen, PIPE, STDOUT


def check_diff(argv):
    """diff checker for SOG buildbot.

    diff the specified files and check the result.  If the files
    differ only in the content of the *RunDateTime header line, or
    Date, or 1 or more output filename lines, exit with status 0 and
    no output.  Otherwise, exit with status 1 and the diff as output.
    """
    if len(argv) != 3:
        return_code = 1
        print 'usage: %(prog)s file1 file2' % {'prog':argv[0]}
        raise SystemExit(return_code)
    # diff the files
    proc = Popen(['diff', argv[1], argv[2]], stderr=PIPE, stdout=PIPE)
    diff, stderr = proc.communicate()
    split_diff = diff.split('\n')
    if not split_diff[-1]: split_diff.pop()
    diff_lines = len(split_diff)
    # Check the result
    if stderr.strip():
        # diff error
        return_code = 1
        diff = stderr
    elif diff_lines == 0:
        # empty diff
        return_code = 0
    else:
        # Might be timeseries, profile, halocline, or Hoffmueller
        # files that differ only by their *RunDateTime header lines,
        # or stdout files that differ only by Date, or 1 or more
        # output file name lines
        ignore = ('*RunDateTime',
                  'Date',
                  'standard physics timeseries', 
                  'user-defined physics timeseries',
                  'standard biology timeseries',
                  'user-defined biology timeseries',
                  'file for halocline results',
                  'profile file base (datetime will be added)',
                  'file for Hoffmueller results')
        ignored_lines = 0
        for line in split_diff:
            if (line[0] in '123456789'
                or line.startswith('---')
                or line[1:].lstrip().split(' =', 1)[0] in ignore
                or line[1:].lstrip().split(':', 1)[0] in ignore):
                ignored_lines += 1
        if ignored_lines == diff_lines:
            return_code = 0
        else:
            return_code = 1
    # Return the status, and the diff if it was dirty
    if return_code == 1:
        print diff
    raise SystemExit(return_code)


if __name__== '__main__':
    check_diff(sys.argv)
