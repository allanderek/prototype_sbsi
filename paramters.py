"""Implements a module to and command-line utility to check that
   the parameter values that result from an optimisation. In particular
   we check them against the initial parameter conditions to see if
   any of the final parameter values are very close to either limit
   of their specified range, suggesting that perhaps the range is too
   tight."""
import sys
import argparse
import random

import sbml_parametiser

class Parameter:
  """The parameter class describes the description of a parameter
     to be optimised by the framework"""
  def __init__(self, name, default_value, low, high):
    self.name = name
    self.low = low
    self.high = high
    self.mutation_probability = 0.5
    self.default_value = default_value

  def set_mutation_probability(self, mut_prob):
    """Set the mutation probability to a new value"""
    self.mutation_probability = mut_prob

  def should_mutate(self):
    """Returns true if we should mutate this value of this parameter.
       This involves choosing a random number and comparing that to
       the mutation probability. In other words this is not something
       intrinsic to the parameter itself."""
    coin = random.random() 
    return coin < self.mutation_probability

  def get_value_range_size(self):
    """Return the size of the range from the lowest to the highest
       value in the allowed for this parameter"""
    return self.high - self.low

  def get_random_value_full_range(self):
    """Return a random value within the full range of the parameter"""
    return random.uniform(self.low, self.high)
 
def get_init_param_parameters(filename):
  """Open and parse the initial parameters file into a list of
     optimisable parameter definitions"""
  paramfile = open(filename, "r")
  # This needs to be somewhat more forgiving
  parameters = []
  for line in paramfile:
    columns  = line.split("\t")
    name     = columns[0]
    low      = float(columns[1])
    high     = float(columns[2])
    begin    = float(columns[3])
    param    = Parameter(name, begin, low, high)
    parameters.append(param)

  mut_prob = 1.0 / float(len(parameters))
  for param in parameters:
    param.set_mutation_probability(mut_prob)

  paramfile.close()
  return parameters


def check_parameters(params, best_params, arguments):
  """Checks initial parameter set up with the final best_params
     (see check_param_files) for more
  """
  tolerance = arguments.tolerance
  for param in params:
    if param.name not in best_params:
      print (param.name + " not in the best params definition") 
    else:
      upper_warn = param.high - (param.high * tolerance)
      # We should be careful, what if param.lower == 0?
      lower_warn = param.low + (param.low * tolerance)
      value = best_params[param.name]
      if value < lower_warn:
        print (param.name + " is close to the lower limit")
        print ("   value = " + str(value))
        print ("   limit = " + str(param.low))
      if value > upper_warn:
        print (param.name + " is close to the upper limit")
        print ("   value = " + str(value))
        print ("   limit = " + str(param.high))
 
def check_param_files(initparam_filename, 
                      bestparams_filename, 
                      arguments):
  """Check the two parameter files against each other, note they
     are not the same kind of parameter file, the first is a
     parameter initialisation file and the second is a set of 
     best params, the idea is to check that the best params are not
     stuck right up against the limit set by the initial parameter
     setup
  """
  params = get_init_param_parameters(initparam_filename)
  best_dict = dict()
  sbml_parametiser.parse_param_file(bestparams_filename, best_dict)
  check_parameters(params, best_dict, arguments)
  
def run():
  """perform the banalities of command-line argument processing and
     then go ahead and compare the parameter results to the
     initial parameter settings
  """ 
  description = "Print out an outline of an SBML file"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="the input files, should be exactly two")
  parser.add_argument("--tolerance", action="store", 
                      type=float, default=0.01,
    help="A value within 'value*tolerance' is considered equal")
  arguments = parser.parse_args()

  if len(arguments.filenames) != 2:
    print("Must provide two input files:")
    print("  an initial parameter file")
    print("  and a best params file")
    sys.exit(1)

  initparams_filename = arguments.filenames[0]
  bestparams_filename = arguments.filenames[1]

  check_param_files(initparams_filename, 
                    bestparams_filename, 
                    arguments)

if __name__ == "__main__":
  run()
