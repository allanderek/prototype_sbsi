"""A script to add a noise function to time series data, generally useful
   for having a fake/dry-run of an optimisation problem. One can start from
   known parameters, generate timecourse data, use this script to add
   noise to that data and then attempt to see if you can optimise the model
   to obtain the original known parameters. If not then the problem is
   likely under-constrained
"""
import argparse
import sys

import timeseries
import random
import utils

def run():
  """perform the banalities of command-line argument processing and
     then go ahead and compare the parameter results to the
     initial parameter settings
  """ 
  description = "Add noise to a timeseries"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="the input files, should be exactly two")
  parser.add_argument('--column', action=utils.ListArgumentAction,
                      help="Specify a column to be in the output")
  parser.add_argument('--mcolumn', action=utils.ListArgumentAction,
                      help="Specify a column not to be in the output")
  parser.add_argument('--output-file', action='store',
                      help="Specify an output file location")

  arguments = parser.parse_args()

  if len(arguments.filenames) < 1:
    print ("Must provide at least one timeseries to add noise to")
    sys.exit(1)

  timecourse_file = arguments.filenames[0]
  timecourse = timeseries.get_timecourse_from_file(timecourse_file)

  all_names = timecourse.get_column_names()
  used_names = utils.get_non_ignored(all_names,
                                     arguments.column,
                                     arguments.mcolumn)
  for name in all_names:
    if name not in used_names:
      timecourse.remove_column(name)

  timecourse.apply_noise_function(dream_noise_function)
 
  output_filename = "noise_timecourse.csv" 
  if arguments.output_file:
    output_filename = arguments.output_file
  newfile = open (output_filename, "w")
  timecourse.write_to_file(newfile)
  newfile.close()

def dream_noise_function(orig_value):
  """Transforms the given value using the noise function used in the
     DREAM estimation of parameters competition of 2011"""
  guass_1 = random.gauss(0, 1)
  guass_2 = random.gauss(0, 1) 

  noise_value = orig_value + (0.1 * guass_1) + (0.2 * guass_2 * orig_value)
  return max(0, noise_value)

if __name__ == "__main__":
  run()
