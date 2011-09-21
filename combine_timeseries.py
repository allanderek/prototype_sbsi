"""A simple script to combine several timecourses together"""
import sys
import argparse
import timeseries
import plotcsv
import utils


def main():
  """Create the command-line parser, parse in all the timecourse files
     adding them together into one overall timecourse and finally out
     put the result to standard out"""
  description = """Add together a bunch of timeseries"""
  parser = argparse.ArgumentParser(add_help=True,
                                   description=description)
   # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="A function definition file to translate")
  parser.add_argument('--column', action=utils.ListArgumentAction,
                      help="Specify a column to be plotted")
  parser.add_argument('--mcolumn', action=utils.ListArgumentAction,
                      help="Specify a column not to be plotted")
  parser.add_argument('--output-separator', action='store',
                      help="Specify the separator to use in the output")
 
  arguments = parser.parse_args()

  timecourse = timeseries.get_timecourse_from_file(arguments.filenames[0])
  for filename in arguments.filenames[1:]:
    additional_timecourse = timeseries.get_timecourse_from_file(filename)
    timecourse.add_timeseries(additional_timecourse)

  for column_name in timecourse.get_column_names():
    if not plotcsv.should_plot_column(arguments, column_name):
      timecourse.remove_column(column_name)

  separator = arguments.output_separator
  if separator == "tab":
    separator = "\t"
  timecourse.write_to_file(sys.stdout, separator=separator)


 
if __name__ == "__main__":
  main()
