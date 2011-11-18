"""A silly module to remove any columns from a timeseries which are
   entirely zeros. Some matlab enthusiasts use this to mean: no data
"""
import argparse
import timeseries

def process_file(filename):
  """Do the work for a single file"""
  timecourse = timeseries.get_timecourse_from_file(filename)

  timecourse.remove_zero_columns()
   
  # note, now it's not actually a new file, but the current one
  # so that we are actually overwriting the file 
  newfile = open (filename, "w")
  timecourse.write_to_file(newfile)
  newfile.close()



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
  arguments = parser.parse_args()

  for filename in arguments.filenames:
    process_file(filename)


if __name__ == "__main__":
  run()
