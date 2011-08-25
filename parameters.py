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
    # The search algorithm may change the 'default_value', however
    # it should never change initial_value, in this way we can check
    # the best params eventually reached with those of the initial
    # value settings.
    self.initial_value = default_value
    # If we are doing a structured experiment then we will need a
    # step size, this will by default be one tenth of the range from
    # low to high value. But can be overridden by a 5th column in the
    # init params file.
    self.step_size = (high - low) / 10.0

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
    # You can have an empty line or one which has no column data
    if len(columns) < 2:
      continue
    name     = columns[0]
    low      = float(columns[1])
    high     = float(columns[2])
    begin    = float(columns[3])
    param    = Parameter(name, begin, low, high)

    # The fifth column is optional so we only interpret it if it is there.
    if len(columns) >= 5:
      column_4 = columns[4].lstrip()
      if column_4:
        param.step_size = float(column_4)

    # Add this parameter to the list of returned parameters
    parameters.append(param)
 
  mut_prob = 1.0 / float(len(parameters))
  for param in parameters:
    param.set_mutation_probability(mut_prob)

  paramfile.close()
  return parameters


def check_parameter(param, value, arguments):
  """Checks a single parameter against the value assigned to it.
     Response with a FailedCheckResult if the parameter is too close
     to either of the two range limits. Responds with None if the
     parameter is fine.
  """
  tolerance = arguments.tolerance
  # The set of acceptable values is everything not too close to
  # the lower limit and also not too close to the upper limit.
  # However it's difficult to define 'too close', especially when
  # the upper and lower limits might be orders of magnitude apart.
  # So we take the range of values that are low enough, to be everything
  # below the initial value. Then we take the difference between the
  # initial value and the upper range and multiply that by the tolerance.
  # So let's say 1% is the tolerance. So we warn about values chosen
  # within the top 1% of everything between the initial value and the
  # upper limit. We do the similar thing for the lower limit. 
  # Notice that this means we might have different allowances for the
  # upper and lower limits (so a value might be allowed to be closer to
  # one than the other). This makes sense, since 1 away from 1000 is quite
  # close, whereas 1 away from 1 is pretty far away.
  # As an example consider:
  # x 10.0 100.0 1.0
  # If our tolerance is 1% then we can be as close as
  # (10.0 - 1.0) * 0.01 = 0.09 to the lower limit 1.0
  # but we must be at least (100.0 - 10.0) * 0.01 = 0.9 away from the
  # upper limit 100.0, so our range of exceptable values would be:
  # 1.09 - 99.1 not inclusive.
  upper_range = param.high - param.initial_value
  upper_warn  = param.high - (upper_range * tolerance)
  lower_range = param.initial_value - param.low
  lower_warn  = param.low + (lower_range * tolerance)

  if value <= lower_warn:
    result = FailedCheckResult(param, value)
    result.too_low = True
    return result
  # Note because we return when the value is too low, this means that
  # the value cannot be both too high and too low, this seems like an
  # acceptable deficiency since in that case either the tolerance is
  # too large or the range is too small.
  if value >= upper_warn:
    result = FailedCheckResult(param, value)
    result.too_high = True
    return result

  # Finally we check if a value is unchanged from its initial value
  # this is not so much to be considered a failure, since it might simply
  # be that the default value is one found in the literature and it
  # happens to be quite on the money.
  if value == param.initial_value:
    result = FailedCheckResult(param, value)
    result.unchanged = True
    return result

  # If everything is okay we return the None.
  return None

class FailedCheckResult:
  """A class for representing the result of failing a parameter check.
     self.value should be none if the parameter is not mentioned in
     the best params file.
  """
  def __init__(self, param, value):
    self.param = param
    # value might be none to indicate that it is not in the
    # best parameter dictionary.
    self.value = value
    self.too_high = False
    self.too_low = False
    self.unchanged = False

def check_parameters(params, best_params, arguments):
  """Checks initial parameter set up with the final best_params
     (see check_param_files) for more
  """
  failed_results = []
  for param in params:
    if param.name not in best_params:
      fail_result = FailedCheckResult(param, None)
      failed_results.append(fail_result)
    else:
      value = best_params[param.name]
      fail_result = check_parameter(param, value, arguments)
      if fail_result:
        failed_results.append(fail_result)
  return failed_results
 
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
  return check_parameters(params, best_dict, arguments)
  
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

  failed_results = check_param_files(initparams_filename, 
                                     bestparams_filename, 
                                     arguments)
  if len(failed_results) == 1:
    print ("No parameters were too close to their range limits")
  else:
    for fail_result in failed_results:
      param = fail_result.param
      if not fail_result.value:
        print (param.name + " is not in the best params file(s)")
      elif fail_result.unchanged:
        print (param.name + " has not changed from the initial value")
      else:  
        print (param.name + " is too close to its range limits")
        print ("   high limit: " + str(param.high))
        print ("   low  limit: " + str(param.low))
        print ("   value     : " + str(fail_result.value))
        print ("   too high  : " + str(fail_result.too_high))
        print ("   too low   : " + str(fail_result.too_low))

if __name__ == "__main__":
  run()
