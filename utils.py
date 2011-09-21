"""A simple module to put some simple utility functions in"""
import os.path
import argparse

def change_filename_ext(filename, new_ext):
  """Returns a new file name based on the first but with the extension
     replaced by the given extension. The new_ext should include the
     '.' separator if that is desired"""
  basename = os.path.splitext(filename)[0]
  return basename + new_ext
 
def get_new_directory(desired_name):
  """Returns a newly created directory with the desired name. If that
     directory already exists it appends a number on to the name"""
  dir_name = desired_name
  number = 0
  while os.path.exists(dir_name):
    number += 1
    dir_name = desired_name + "_" + str(number)  
    
  # create directory and return path
  os.makedirs(dir_name)
  return dir_name


class ListArgumentAction(argparse.Action):
  """ A simple argument action to allow a list of arguments to be
      separated by a comma. This allows slightly more concise command
      lines. For example --column p1 --column p2 --column p3 can be
      reduced to --column p1,p2,p3
      To be used as: 
  parser.add_argument('--column', action=utils.ListArgumentAction,
                      help="Specify a column to be plotted")
  """
  def __call__(self, parser, namespace, values, option_string=None):
    """The action called when the parser encounters the associated
       command-line argument.
    """
    # First of all get the list of values specified as this argument
    list_values = values.split(",")
    # Then get the current list already specified, this allows for
    # example: --column p1,p2,p3 --column p4,p5,p6 or more likely just
    # that some user does not know about the comma separator and does
    # --column p1 --column p2 --column p3 ... etc
    current = getattr(namespace, self.dest)
    if current:
      # If there are values there already extend the list
      current.extend(list_values)
    else: 
      # otherwise set the list to the current set of values
      setattr(namespace, self.dest, list_values)


