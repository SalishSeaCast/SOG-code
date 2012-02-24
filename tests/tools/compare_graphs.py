"""SOG buildbot timeseries comparison graph generator.

Usage:
  python compare_graphs.py run_results_dir ref_results_dir std_phys_file std_bio_file std_chem_file graph_file
"""
"""
test_compare_graphs.py is the test suite for this module.

:Author: Doug Latornell <dlatorne@eos.ubc.ca>
"""
# Standard library:
import getopt
import os
import sys
# Dateutil extensions to standard library datetime:
import dateutil.parser
# NumPy
import numpy
# matplotlib
import matplotlib.figure
import matplotlib.backends.backend_pdf
from matplotlib.font_manager import FontProperties


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


class Relation(object):
    """A Relation object has a pair of NumPy arrays containing the
    independent and dependent data values of a data set. It also has
    attributes that contain the filespec from which the data is read,
    and the units of the data arrays.

    This is a base class for implementing specific types of relations
    such as timeseries, or profiles. This class provides a
    :meth:`read_data` method, but expects classes that inherit from it
    to implement a :meth:`_read_header` method that returns lists of
    field names and field units read from the SOG data file header.
    """
    def __init__(self, datafile):
        """Create a Relation instance with its datafile attribute
        initialied.
        """
        self.datafile = datafile

    def _read_header(self, fobj):
        """This method is expected to be implemented by classes based
        on Relation.

        It must return lists of field names and field units from the
        SOG data file header.
        """
        raise NotImplementedError

    def read_data(self, indep_field, dep_field):
        """Read the data for the specified independent and dependent
        fields from the data file.

        Sets the indep_data and dep_data attributes to NumPy arrays,
        and the indep_units and dep_units attributes to units strings
        for the data fields.
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
                msg = ("Invalid field name: {0}, choices are:\n {1}"
                       .format(field, field_names))
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


class SOG_Timeseries(Relation):
    def _read_header(self, fobj):
        """Read a SOG timeseries file header, return the field_names
        and field_units lists, and set attributes with the
        run_datetime and initial_CTD_datetime values.
        """
        for line in fobj:
            line = line.strip()
            if line.startswith('*FieldNames:'):
                # Drop the *FieldNames: label and keep the
                # comma-delimited list
                field_names = line.split(': ', 1)[1].split(', ')
            if line.startswith('*FieldUnits:'):
                # Drop the *FieldUnits: label and keep the
                # comma-delimited list
                field_units = line.split(': ', 1)[1].split(', ')
            if line.startswith('*RunDateTime:'):
                datetime_str = ' '.join(line.split()[1:])
                self.run_datetime = dateutil.parser.parse(datetime_str)
            if line.startswith('*InitialCTDDateTime:'):
                datetime_str = ' '.join(line.split()[1:])
                self.initial_CTD_datetime = dateutil.parser.parse(datetime_str)
            if line.startswith('*EndOfHeader'):
                break
        return field_names, field_units


class GraphPage(object):
    def __init__(self, figsize=(8, 14), subplot_rows=4, subplot_cols=1):
        """Create a GraphPage instance with various attributes initialized.

        :arg figsize:
        :type figsize: 2-tuple of ints

        :arg subplot_rows:
        :type subplot_rows: int

        :arg subplot_cols:
        :type subplot_cols: int
        """
        self.fig = matplotlib.figure.Figure(figsize=figsize)
        self.subplot_rows = subplot_rows
        self.subplot_cols = subplot_cols
        self.legend_data = []

    def save_pdf(self, filename):
        """Save the graph page as a PDF.
        """
        filename = (filename if filename.lower().endswith('.pdf')
                    else filename + '.pdf')
        canvas = matplotlib.backends.backend_pdf.FigureCanvasPdf(self.fig)
        canvas.print_figure(filename)

    def _add_subplot(self, num):
        """Add a subplot to the figure.
        """
        return self.fig.add_subplot(self.subplot_rows, self.subplot_cols, num)

    def one_axis_subplot(self, num, field, run_ts, ref_ts, colour='k'):
        """Create a subplot with a single y-axis on the left.
        """
        ax_left = self._add_subplot(num)
        # Plot ref results as a black line with markers
        ref_ts.read_data('time', field)
        ref_line = ax_left.plot(
            ref_ts.indep_data, ref_ts.dep_data, 'k+-',
            label='Ref: {0.run_datetime}'.format(ref_ts))
        # Plot run results as a coloured line. Plot run results 2nd so
        # that the coloured line is visible when the run and ref
        # results are identical.
        run_ts.read_data('time', field)
        run_line = ax_left.plot(
            run_ts.indep_data, run_ts.dep_data, color=colour,
            label='Run: {0.run_datetime}'.format(run_ts))
        # Anchor x-axis at 0
        ax_left.set_xlim(xmin=0)
        # y-axis label
        ax_left.set_ylabel(
            '{0} [{1.dep_units}]'.format(field.title(), run_ts),
            color=colour, size='x-small')
        # Legend
        if not self.legend_data:
            self.legend_data = {
                'lines': (run_line, ref_line),
                'labels': ('Run: {0.run_datetime}'.format(run_ts),
                           'Ref: {0.run_datetime}'.format(ref_ts))}

    def two_axis_subplot(self, num, fields, run_ts, ref_ts,
                         colours=('b', 'g')):
        """Create a subplot with both left and right y-axis and one
        field plotted on each.
        """
        ax_left = self._add_subplot(num)
        ax_right = ax_left.twinx()
        matplotlib.axes.Axes(self.fig, ax_left.get_position(), sharex=ax_right)
        # Plot ref results as black lines with markers
        ref_ts.read_data('time', fields[0])
        ref_line = ax_left.plot(
            ref_ts.indep_data, ref_ts.dep_data, 'k+-',
            label='Ref: {0.run_datetime}'.format(ref_ts))
        ref_ts.read_data('time', fields[1])
        ax_right.plot(
            ref_ts.indep_data, ref_ts.dep_data, 'k+-',
            label='{0.run_datetime}'.format(ref_ts))
        # Plot run results as coloured lines. Plot run results 2nd so
        # that the coloured line is visible when the run and ref
        # results are identical.
        run_ts.read_data('time', fields[0])
        run_line = ax_left.plot(
            run_ts.indep_data, run_ts.dep_data, color=colours[0],
            label='Run: {0.run_datetime}'.format(run_ts))
        # Left y-axis label
        ax_left.set_ylabel(
            '{0} [{1.dep_units}]'.format(fields[0].title(), run_ts),
            color=colours[0], size='x-small')
        run_ts.read_data('time', fields[1])
        ax_right.plot(
            run_ts.indep_data, run_ts.dep_data, '-', color=colours[1],
            label='Run: {0.run_datetime}'.format(run_ts))
        # Right y-axis label
        ax_right.set_ylabel(
            '{0} [{1.dep_units}]'.format(fields[1].title(), run_ts),
            color=colours[1], size='x-small')
        # Anchor x-axis at 0
        ax_left.set_xlim(xmin=0)
        # Legend
        if not self.legend_data:
            self.legend_data = {
                'lines': (run_line, ref_line),
                'labels': ('Run: {0.run_datetime}'.format(run_ts),
                           'Ref: {0.run_datetime}'.format(ref_ts))}

    def legend(self):
        self.fig.legend(self.legend_data['lines'], self.legend_data['labels'],
                        'upper left', prop=FontProperties(size=8.0))


def make_graphs(run_results_dir, ref_results_dir, phys_results_file,
                bio_results_file, chem_results_file, graph_file):
    """Create a PDF with 4 graphs showing run and ref results for:

       * surface temperature
       * surface salinity
       * mixing layer depth
       * suface nitrate concentration
       * surface micro phytoplankton biomass
       * surface dissolved inorganic carbon concentration
       * surface dissolved oxygen concentration
    """
    pg = GraphPage()
    run_ts = SOG_Timeseries(os.path.join(run_results_dir, phys_results_file))
    ref_ts = SOG_Timeseries(os.path.join(ref_results_dir, phys_results_file))
    pg.two_axis_subplot(
        1, ('surface temperature', 'surface salinity'),
        run_ts, ref_ts, colours=('r', 'b'))
    pg.one_axis_subplot(
        2, 'mixing layer depth', run_ts, ref_ts, colour='m')
    run_ts = SOG_Timeseries(os.path.join(run_results_dir, bio_results_file))
    ref_ts = SOG_Timeseries(os.path.join(ref_results_dir, bio_results_file))
    pg.two_axis_subplot(
        3, ('surface nitrate concentration',
            'surface micro phytoplankton biomass'),
        run_ts, ref_ts, colours=('c', 'g'))
    run_ts = SOG_Timeseries(os.path.join(run_results_dir, chem_results_file))
    ref_ts = SOG_Timeseries(os.path.join(ref_results_dir, chem_results_file))
    pg.two_axis_subplot(
        4, ('surface DIC concentration', 'surface oxygen concentration'),
        run_ts, ref_ts, colours=('brown', 'orange'))
    pg.legend()
    pg.save_pdf(graph_file)


def main(argv=[__name__]):
    try:
        try:
            opts, args = getopt.getopt(argv[1:], 'h', ['help'])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ('-h', '--help'):
                print __doc__
                sys.exit(0)
        try:
            run_results_dir, ref_results_dir = args[:2]
            phys_results_file, bio_results_file, chem_results_file = args[2:5]
            graph_file = args[5]
        except IndexError:
            raise Usage('missing 1 or more arguments')
        make_graphs(run_results_dir, ref_results_dir, phys_results_file,
                    bio_results_file, chem_results_file, graph_file)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, 'for help use --help'
        sys.exit(2)


if __name__ == '__main__':
    main(sys.argv)
