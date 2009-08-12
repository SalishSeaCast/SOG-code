"""Test suite for SOG buildbot timeseries comparison graph generator.

Use `python test_compare_graphs.py` to run the test suite.

:Author: Doug Latornell <dlatorne@eos.ubc.ca>
:Created: 2009-07-28
"""
import numpy
import StringIO
import sys
import tempfile
import unittest
from datetime import datetime
from compare_graphs import Relation, SOG_Timeseries


class TestRelation(unittest.TestCase):
    """Unit tests for Relation class.
    """
    def setUp(self):
        # Create a test data file
        self.indep_data = [0.0, 0.5]
        self.dep_data = [12.1, 12.05]
        self.test_file = tempfile.NamedTemporaryFile()
        for pair in zip(self.indep_data, self.dep_data):
            self.test_file.write("%f %f\n" % pair)
        self.test_file.flush()
        # Mock the _read_header method
        def _read_header(fobj):
            field_names = ['depth', 'surface temperature']
            field_units = ['m', 'deg C']
            return field_names, field_units
        self._read_header = _read_header


    def test_init(self):
        """creation of Relation instance sets expected attributes
        """
        rel = Relation('test_datafile')
        self.assertEqual(rel.datafile, 'test_datafile')


    def test_read_header(self):
        """_read_header raises NotImplementedError
        """
        # Delete the _read_header mock set by self.setUp
        del self._read_header
        rel = Relation('test_datafile')
        self.assertRaises(NotImplementedError, rel._read_header, file)


    def test_read_data_bad_file(self):
        """read_data method raises SystemExit on read error
        """
        rel = Relation('badfile')
        sys.stdout = StringIO.StringIO()
        self.assertRaises(SystemExit, rel.read_data,
                          'depth', 'surface temperature')
        self.assertEqual(sys.stdout.getvalue().strip(),
                         "[Errno 2] No such file or directory: 'badfile'")
        sys.stdout = sys.__stdout__


    def test_read_data_bad_field_name(self):
        """read_data method raises SystemExit on bad field name
        """
        rel = Relation(self.test_file.name)
        rel._read_header = self._read_header
        sys.stdout = StringIO.StringIO()
        self.assertRaises(SystemExit, rel.read_data,
                          'deepness', 'surface temperature')
        self.assertEqual(sys.stdout.getvalue().strip(),
                         "Invalid field name: deepness, choices are:\n "
                         "['depth', 'surface temperature']")
        sys.stdout = sys.__stdout__


    def test_read_data(self):
        """read_data method sets independent & dependent data arrays
        """
        rel = Relation(self.test_file.name)
        rel._read_header = self._read_header
        rel.read_data('depth', 'surface temperature')
        self.assertEqual(rel.indep_units, 'm')
        self.assertEqual(rel.dep_units, 'deg C')
        numpy.testing.assert_equal(rel.indep_data, numpy.array(self.indep_data))
        numpy.testing.assert_equal(rel.dep_data, numpy.array(self.dep_data))


class TestSOG_Timeseries(unittest.TestCase):
    """Unit tests for SOG_Timeseries class.
    """
    def setUp(self):
        # Create a test data file
        self.test_file = tempfile.NamedTemporaryFile()
        self.test_file.write('''/
! SOG code standard time series output from physics model
! Time series of iteration count for each time step; mixing layer depth;
! velocity components, temperature, & salinity; at surface, and averaged over
! top 3 m of water column and surface par
*FromCode: $Source$
*RunDateTime: 2009-07-27 07:14:24
*InitialCTDDateTime: 2004-10-19 12:22:00
*FieldNames: time, iteration count, mixing layer depth, surface u velocity, 3 m avg u velocity, surface v velocity, 3 m avg v velocity, surface temperature, 3 m avg temperature, surface salinity, 3 m avg salinity, surface PAR
*FieldUnits: hr since 2004-10-19 00:00:00 LST, None, m, m/s, m/s, m/s, m/s, deg C, deg C, None, None, W/m2
*EndOfHeader
   12.3667    4    1.2622    0.0046    0.0011   -0.0491   -0.0151   12.0920   12.1608   22.8330   23.2946  104.9434
   12.6167    2    1.2915    0.0057    0.0017   -0.0645   -0.0266   12.1199   12.1743   22.8852   23.2923  103.4063/
''')
        self.test_file.flush()


    def test_parent_class(self):
        """SOG_Timeseries subclasses Relation
        """
        self.assertTrue(issubclass(SOG_Timeseries, Relation))


    def test_read_header(self):
        """_read_header returns expected lists & sets expected attributes
        """
        ts = SOG_Timeseries(self.test_file.name)
        f = open(ts.datafile, 'r')
        field_names, field_units = ts._read_header(f)
        self.assertEqual(field_names,
                         'time, iteration count, mixing layer depth, '
                         'surface u velocity, 3 m avg u velocity, '
                         'surface v velocity, 3 m avg v velocity, '
                         'surface temperature, 3 m avg temperature, '
                         'surface salinity, 3 m avg salinity, '
                         'surface PAR'.split(', '))
        self.assertEqual(field_units,
                         'hr since 2004-10-19 00:00:00 LST, None, m, m/s, '
                         'm/s, m/s, m/s, deg C, deg C, None, None, '
                         'W/m2'.split(', '))
        self.assertEqual(ts.run_datetime, datetime(2009, 7, 27, 7, 14, 24))
        self.assertEqual(ts.initial_CTD_datetime,
                         datetime(2004, 10, 19, 12, 22, 0))


    def test_read_data(self):
        """read_data method sets independent & dependent data arrays
        """
        ts = SOG_Timeseries(self.test_file.name)
        ts.read_data('time', 'mixing layer depth')
        self.assertEqual(ts.indep_units, 'hr since 2004-10-19 00:00:00 LST')
        self.assertEqual(ts.dep_units, 'm')
        numpy.testing.assert_equal(ts.indep_data,
                                   numpy.array([12.3667, 12.6167]))
        numpy.testing.assert_equal(ts.dep_data,
                                   numpy.array([1.2622, 1.2915]))


if __name__ == '__main__':
    unittest.main()

