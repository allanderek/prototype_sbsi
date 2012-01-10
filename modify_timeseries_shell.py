""" We have certain scripts to modify timeseries, for example to add
    noise in simulation of experimental observations and to add
    derivatives for other purposes. These all share quite a large
    portion of their code so we generalise that here.
"""
import argparse
import sys

import timeseries
import utils

def create_arguments_parser():
  """ 
     Create an arguments parser to be used as a basis for
     all argument parsers for modifying timeseries scripts.
  """
  parser = argparse.ArgumentParser(add_help=False)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="the input files, should be exactly two")
  parser.add_argument('--column', action=utils.ListArgumentAction,
                      help="Specify a column to be in the output")
  parser.add_argument('--mcolumn', action=utils.ListArgumentAction,
                      help="Specify a column not to be in the output")
  utils.add_output_file_arg(parser)
  return parser



def run(arguments_parser, get_new_timecourse, needs_arguments):
  """
     The shell of the main method, essentially does all the busy work
     and accepts a function which will given a read-in timecourse,
     return a new timecourse.
  """
  arguments = arguments_parser.parse_args()
  if len(arguments.filenames) < 1:
    print ("Must provide at least one timeseries to add noise to")
    sys.exit(1)

  # We should do all this, once for each timecourse file.
  timecourse_file = arguments.filenames[0]
  timecourse = timeseries.get_timecourse_from_file(timecourse_file)

  all_names = timecourse.get_column_names()
  used_names = utils.get_non_ignored(all_names,
                                     arguments.column,
                                     arguments.mcolumn)
  for name in all_names:
    if name not in used_names:
      timecourse.remove_column(name)

  if needs_arguments:
    new_timecourse = get_new_timecourse(arguments, timecourse)  
  else:
    new_timecourse = get_new_timecourse(timecourse)  

   
  output_filename = utils.get_output_filename(timecourse_file, arguments,
                                              "_modified.csv")
  if output_filename == "stdout":
    output_file = sys.stdout
  else:
    output_file = open(output_filename, "w")

  # Write the new time course to the file
  new_timecourse.write_to_file(output_file)

  if output_filename != "stdout":
    output_file.close()

