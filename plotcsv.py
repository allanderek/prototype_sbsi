"""A simple python script to use gnuplot to plot a csv file.
   It may also accept a separator for example tab rather than
   a comma.
"""
import os
import sys
from subprocess import Popen
import argparse


def get_separator(filename, arguments):
  """Obtain the separator character either from the the command-line
     arguments or by looking at the file/using the default"""
  separator = arguments.sep
  if not separator:
    datafile = open(filename, "r")
    separator = ""
    # We're not catching StopIteration here since I guess
    # if don't find any of the possible separators we're humped.
    while separator == "":
      line = datafile.next()
      if "," in line:
        separator = ","
      elif "\t" in line:
        separator = "\t"
    datafile.close()
  return separator


def obtain_headers(datafile, separator):
  """Obtain the column names from the csv file"""
  csvfile = open(datafile, "rb")
  headrow = csvfile.next()
  while separator not in headrow:
    headrow = csvfile.next()

  headers = headrow.split(separator)

  csvfile.close()
  return headers


def set_gnuplot_options(gnuplotfile, arguments):
  """Set up some of the available options in the gnuplot file"""
  title = arguments.title
  if title:
    gnuplotfile.write("set title \"" + title + "\"\n")

  x_range = arguments.x_range
  if x_range:
    gnuplotfile.write("set xrange " + x_range + "\n")
  y_range = arguments.y_range
  if y_range:
    gnuplotfile.write("set yrange " + y_range + "\n")
  key = arguments.key
  if key:
    gnuplotfile.write("set key " + key + "\n")


class NameAliaser:
  """A simple class which effectively ensures we are using unique names.
     Everytime we ask for a unique name, given a desired name, we check
     if the desired name is already in use. If not then we can simply
     return the desired name. If so we add a number to the end of it"""
  def __init__(self):
    self.dictionary = dict()

  def get_unique_name(self, name):
    if name in self.dictionary:
      # If the name is already in the dictionary increase the count
      # and make a new name based on the given name plus the current
      # count.
      num_entries = self.dictionary[name]
      self.dictionary[name] = num_entries + 1

      generated_name = name + "_" + str(num_entries)

      # To absolutely ensure that we have a unique name, we recursively
      # call 'get_unique_name' on the generated name.
      return self.get_unique_name(generated_name)
    else:
      self.dictionary[name] = 1
      return name

def create_gnuplot_file(basename, arguments, datafiles):
  """From the headers obtained from the csv file, write a
     gnuplot file which will read the csv file and plot a timeseries"""

  gnufilename = basename + ".gnuplot"
  gnuplotfile = open (gnufilename, "w")

  epsfilename = basename + ".eps"
  gnuplotfile.write("set term post eps color\n")
  gnuplotfile.write("set output \"" + epsfilename + "\"\n")
  # This assumes that all the files have the same separator
  separator = get_separator(datafiles[0], arguments)
  gnuplotfile.write("set datafile separator \"" + separator + "\"\n")
  
  set_gnuplot_options(gnuplotfile, arguments)

  line_style = arguments.linestyle
  if not line_style:
    line_style = "lp"

  line_prefix = "plot "
  header_aliaser = NameAliaser()
  for datafile in datafiles:
    headers = obtain_headers(datafile, separator)
    columns = arguments.column
    mcolumns = arguments.mcolumn

    # Just a bit of debugging print(column_names)
  
    # The first column must be printed out a little differently to
    # the others, since we need to separate the lines, so we set
    # the prefix to what it will be for the first line, and then whenever
    # we print one out we just set it to what the prefix should be for
    # any other line (which is to end the previous line and indent).
    quote_data_file = "\"" + datafile + "\"" 
    for i in range (1, len(headers)):
      header = headers[i].lstrip().rstrip()
      header_title = header_aliaser.get_unique_name(header)
      if header != "Time" and ((not columns or header in columns) and
                               (not mcolumns or header not in mcolumns)):
        line = (line_prefix + quote_data_file + " using 1:" +
                str(i + 1) + " w " + line_style + 
                " title '" + header_title + "'")
        # If it's not currently the first line then this will simply
        # set the prefix to what it already is, so no harm done.
        line_prefix = ", \\\n  "
        gnuplotfile.write(line) 

  gnuplotfile.write("\n")
  gnuplotfile.close()
  return (gnufilename, epsfilename)
 

def run ():
  """Simply do all the work"""
  description = "Plot csv files using gnuplot"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="a csv file to plot")
  parser.add_argument('--sep', action='store')
  parser.add_argument('--basename', action='store')
  parser.add_argument('--title', action='store')
  parser.add_argument('--x_range', action='store')
  parser.add_argument('--y_range', action='store')
  parser.add_argument('--key', action='store')
  parser.add_argument('--linestyle', action='store')
  parser.add_argument('--column', action='append')
  parser.add_argument('--mcolumn', action='append')

  arguments = parser.parse_args()

  filenames    = arguments.filenames
  datafile     = filenames[0] 

  # We base the output file name on the first file name
  basename = arguments.basename
  if not basename:
    basename = os.path.splitext(datafile)[0]

  (gnufilename, epsfilename) = create_gnuplot_file(basename,
                                                   arguments,
                                                   filenames)

  # We now also actually run gnuplot
  gnuplot_command = [ "gnuplot", gnufilename ]
  gnuplot_process = Popen(gnuplot_command)
  gnuplot_process.communicate()

  epstopdf_command = [ "epstopdf", epsfilename ]
  epstopdf_process = Popen(epstopdf_command)
  epstopdf_process.communicate()


if __name__ == "__main__":
  run()
