import numpy


class Relation(object):
    """
    """
    def __init__(self, datafile):
        """
        """
        self.datafile = datafile

    
    def read_data(self, indep_field, dep_field):
        """
        """
        try:
            f = open(self.datafile, 'r')
        except IOError, msg:
            return_code = 2
            print msg
            raise SystemExit(return_code)
        # Read the file header
        (field_names, field_units) = self._read_header(f)
        # Translate the field names into column numbers
        for field in (indep_field, dep_field):
            if field not in field_names:
                return_code = 2
                msg = ("Invalid field name: %(field)s, "
                       "choices are:\n %(field_names)s" % vars())
                print msg
                raise SystemExit(return_code)
        indep_col = field_names.index(indep_field)
        dep_col = field_names.index(dep_field)
        # Data units
        self.indep_units = field_units[indep_col]
        self.dep_units = field_units[dep_col]
        # Read the data
        self.indep_data, self.dep_data = [], []
        try:
            for line in f:
                self.indep_data.append(float(line.split()[indep_col]))
                self.dep_data.append(float(line.split()[dep_col]))
        finally:
            f.close()
        # Transform data lists into NumPy arrays
        self.indep_data = numpy.array(self.indep_data)
        self.dep_data = numpy.array(self.dep_data)
