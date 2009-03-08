"""Test suite for SOG buildbot infile reader.

Use `python test_read_infile.py` to run the test suite.

:Author: Doug Latornell <dlatorne@eos.ubc.ca>
:Created: 2009-03-07
"""
import sys
import StringIO
import unittest
import tempfile
from read_infile import read_infile


class Test_check_diff(unittest.TestCase):
    def assertRaises(self, excClass, callableObj, *args, **kwargs):
        """Monkey-patch unittest.assertRaises() to return exception
        instance.
        """
        try:
            callableObj(*args, **kwargs)
        except excClass, exc_instance:
            return exc_instance
        else:
            if hasattr(excClass,'__name__'): excName = excClass.__name__
            else: excName = str(excClass)
            raise self.failureException, "%s not raised" % excName


    def setUp(self):
        """Pre-test setup.
        """
        # Capture stdout
        sys.stdout = StringIO.StringIO()


    def tearDown(self):
        """Post-test tear-down.
        """
        # Restore stdout to default
        sys.stdout = sys.__stdout__


    def test_too_few_args(self):
        """print usage and exit with status=1 if too few arguments.
        """
        argv = 'read_infile.py'.split()
        exception = self.assertRaises(SystemExit, read_infile, argv)
        self.assertEqual(exception.code, 1)
        self.assertEqual(sys.stdout.getvalue(), 
                         'usage: %(prog)s infile key\n' % {'prog':argv[0]})


    def test_too_many_args(self):
        """print usage and exit with status=1 if too many arguments.
        """
        argv = 'read_infile.py infile key extra'.split()
        exception = self.assertRaises(SystemExit, read_infile, argv)
        self.assertEqual(exception.code, 1)
        self.assertEqual(sys.stdout.getvalue(), 
                         'usage: %(prog)s infile key\n' % {'prog':argv[0]})


    def test_bad_file(self):
        """print error and exit with status=1 if file can't be read.
        """
        argv = 'read_infile.py badfile key'.split()
        exception = self.assertRaises(SystemExit, read_infile, argv)
        self.assertEqual(exception.code, 1)
        self.assertEqual(sys.stdout.getvalue(), 
                        "[Errno 2] No such file or directory: 'badfile'\n")


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
        argv = ['read_infile.py', test_file.name, 'haloclinefile']
        exception = self.assertRaises(SystemExit, read_infile, argv)
        self.assertEqual(exception.code, 0)
        self.assertEqual(sys.stdout.getvalue(), """\
profiles/halo-nov04-test.out
""")


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
        argv = ['read_infile.py', test_file.name, 'haloclinefile']
        exception = self.assertRaises(SystemExit, read_infile, argv)
        self.assertEqual(exception.code, 0)
        self.assertEqual(sys.stdout.getvalue(), """\
profiles/halo-nov04-test.out
""")


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
        argv = ['read_infile.py', test_file.name, 'std_phys_ts_out']
        exception = self.assertRaises(SystemExit, read_infile, argv)
        self.assertEqual(exception.code, 0)
        self.assertEqual(sys.stdout.getvalue(), """\
timeseries/std_phys_nov04-test.out
""")


if __name__ == '__main__':
    unittest.main()
