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
from compare_graphs import Relation


class TestRelation(unittest.TestCase):
    """Unit tests for Relation class.
    """
    def setUp(self):
        # Create a test data file
        indep_data = [0.0, 0.5]
        dep_data = [12.1, 12.05]
        self.test_file = tempfile.NamedTemporaryFile()
        for pair in zip(indep_data, dep_data):
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
        self.assertEqual(rel.indep_data, numpy.array(rel.indep_data))
        self.assertEqual(rel.dep_data, numpy.array(rel.dep_data))


if __name__ == '__main__':
    unittest.main()

