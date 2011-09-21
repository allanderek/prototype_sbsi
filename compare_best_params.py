"""
Implements a utility to compare the results of optimisation, in particular
it compares two best parameters files.
"""
import sys
import argparse
import parameters

class ParamComparison:
  """A simple class to hold two values for the same named parameter"""
  def __init__(self, name, first, second):
    self.name = name
    self.first = first
    self.second = second

  def compare_ratio(self):
    """Returns the ratio between the two parameter values"""
    largest = max(self.first, self.second)
    smallest = min(self.first, self.second)
    return largest / smallest
  
def compare_parameters(first, second):
  """Compare the parameters in the two given lists which have equal names
  """
  comparison_list = []
  for param_name in first:
    if param_name in second:
      first_value = first[param_name]
      second_value = second[param_name]
      comparison = ParamComparison(param_name, first_value, second_value)
      comparison_list.append(comparison)

  sorted_list = sorted(comparison_list, key=lambda x : x.compare_ratio())
  for param_comparison in sorted_list:
    print (param_comparison.param_name + "\t" +
           str(param_comparison.first) + "\t" + 
           str(param_comparison.second))

def run():
  """perform the banalities of command-line argument processing and
     then get on with the proper work
  """ 
  description = "Compare results of two separate optimisations"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="the input files, should be exactly two")
  arguments = parser.parse_args()

  if len(arguments.filenames) != 2:
    print ("Must provide at least two best params files")
    sys.exit(1)

  best_params_lists = [ parameters.parse_param_file(f) 
                          for f in arguments.filenames ]

  # We can do better than simply comparing every other file to the
  # first parameter file but for now I'm going to do that.
  first_best_params = best_params_lists[0]
  for other_best_params in best_params_lists[1:]:
    compare_parameters(first_best_params, other_best_params)
  

if __name__ == "__main__":
  run()
