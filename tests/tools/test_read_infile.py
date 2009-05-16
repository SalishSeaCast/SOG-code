"""Test suite for SOG buildbot infile reader.

Use `python test_read_infile.py` to run the test suite.

:Author: Doug Latornell <dlatorne@eos.ubc.ca>
:Created: 2009-03-07
"""
import os
import StringIO
import subprocess
import sys
import tempfile
import unittest
from read_infile import read_infile


# Name of module under test
module = 'read_infile.py'

class Test_read_infile_cli(unittest.TestCase):
    """Tests for the commandline interface to read_infile.
    """
    def test_too_few_args(self):
        """print usage and exit with status=2 if too few arguments.
        """
        cmd = ['python', os.path.join('.', module)]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 2)
        self.assertEqual(stdout.strip(), 
                         'Usage: %(module)s infile key' % globals())


class Test_read_infile_function(unittest.TestCase):
    """Tests for the read_infile function.
    """
    def test_bad_file(self):
        """print error and exit with status=2 if file can't be read.
        """
        cmd = ['python', os.path.join('.', module), 'badfile', 'key']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 2)
        self.assertEqual(stdout.strip(), 
                        "[Errno 2] No such file or directory: 'badfile'")


    def test_1_line_params(self):
        """print outfile and exit with status=0 for files with params on 1 line.
        """
        test_file = tempfile.NamedTemporaryFile()
        test_file.write("""\
"dt"		900	 "time step [s]"
! Implicit solver iteration limit
"max_iter"		30	"implicit solver max iterations"


! Time series output files
"std_phys_ts_out" "timeseries/std_phys_nov04-test.out" "standard physics timeseries"
"user_phys_ts_out" "timeseries/user_phys_nov04-test.out" "user-defined physics timeseries"

! Profiles output
"noprof"	1		"no. of profiles to print"
! *** It would be nice to replace the yr-day and day-sec lists with a list
! *** of dates/times for profiles output
"profday"	303		"yr-day for profile"
"proftime"	43200.		"day-sec for profile"
"haloclinefile"	"profiles/halo-nov04-test.out"	"file for halocline results"
"profile_base"	"profiles/nov04-test"	"profile file base (datetime will be added)"
! Hoffmueller diagram output (a collection of profiles at time intervals
! for contour or colourmap plotting)
"Hoffmueller file"	"profiles/hoff-nov04-test.dat"	"file for Hoffmueller results"\
""")
        test_file.flush()
        cmd = ['python', os.path.join('.', module),
               test_file.name, 'haloclinefile']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 0)
        self.assertEqual(stdout.strip(), "profiles/halo-nov04-test.out")


    def test_key_with_whitespace(self):
        """print outfile and exit with status=0 for key containing whitespace.
        """
        test_file = tempfile.NamedTemporaryFile()
        test_file.write("""\
"dt"		900	 "time step [s]"
! Implicit solver iteration limit
"max_iter"		30	"implicit solver max iterations"


! Time series output files
"std_phys_ts_out" "timeseries/std_phys_nov04-test.out" "standard physics timeseries"
"user_phys_ts_out" "timeseries/user_phys_nov04-test.out" "user-defined physics timeseries"

! Profiles output
"noprof"	1		"no. of profiles to print"
! *** It would be nice to replace the yr-day and day-sec lists with a list
! *** of dates/times for profiles output
"profday"	303		"yr-day for profile"
"proftime"	43200.		"day-sec for profile"
"haloclinefile"	"profiles/halo-nov04-test.out"	"file for halocline results"
"profile_base"	"profiles/nov04-test"	"profile file base (datetime will be added)"
! Hoffmueller diagram output (a collection of profiles at time intervals
! for contour or colourmap plotting)
"Hoffmueller file"	"profiles/hoff-nov04-test.dat"	"file for Hoffmueller results"\
""")
        test_file.flush()
        cmd = ['python', os.path.join('.', module),
               test_file.name, 'Hoffmueller file']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 0)
        self.assertEqual(stdout.strip(), "profiles/hoff-nov04-test.dat")


    def test_multiline_params(self):
        """print outfile and exit with status=0 for multiline param files.
        """
        test_file = tempfile.NamedTemporaryFile()
        test_file.write("""\
"dt"		900	 "time step [s]"
! Implicit solver iteration limit
"max_iter"		30	"implicit solver max iterations"


! Time series output files
"std_phys_ts_out"
        "timeseries/std_phys_nov04-test.out" 
        "standard physics timeseries"
"user_phys_ts_out" "timeseries/user_phys_nov04-test.out" "user-defined physics timeseries"

! Profiles output
"noprof"	1		"no. of profiles to print"
! *** It would be nice to replace the yr-day and day-sec lists with a list
! *** of dates/times for profiles output
"profday"	303		"yr-day for profile"
"proftime"	43200.		"day-sec for profile"
"haloclinefile"	"profiles/halo-nov04-test.out"	"file for halocline results"
"profile_base"	
        "profiles/nov04-test"	
	"profile file base (datetime will be added)"
! Hoffmueller diagram output (a collection of profiles at time intervals
! for contour or colourmap plotting)
"Hoffmueller file"	"profiles/hoff-nov04-test.dat"	"file for Hoffmueller results"\
""")
        test_file.flush()
        cmd = ['python', os.path.join('.', module),
               test_file.name, 'std_phys_ts_out']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 0)
        self.assertEqual(stdout.strip(), "timeseries/std_phys_nov04-test.out")


if __name__ == '__main__':
    unittest.main()
