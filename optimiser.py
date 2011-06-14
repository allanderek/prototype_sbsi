"""
This is a prototype of a new version of the sbsi optimisation
framework. The idea being that we are less coupled within the code
such that we combining together many smaller applications.
"""

import os
from subprocess import Popen, PIPE
import sys
import random
import math
import logging
import configuration

def get_options(option_name, arguments):
  """Gets all options of the form <option_name>=<argument> and returns
     the argument parts."""
  option_string = option_name + "="
  option_offset = len (option_string)
  options = [ option[option_offset:] for option in arguments
                if option.startswith(option_string) ]
  return options

def get_single_option(option_name, arguments):
  """The same as 'get_options' but insists on their being only one such,
     if there are no such arguments return the empty string"""
  options = get_options(option_name, arguments)
  if not options:
    return ""
  elif len(options) > 1:
    print ("Conflicting " + option_name +
           " arguments, choose one, exiting")
    sys.exit(1)
  else:
    return options[0]

def get_int_option_with_default(option_name, arguments, default_value):
  """A further option helper function. This time we are sure that the
     argument give to the option should be an integer number and we
     allow a default to be returned in the case that no such option is
     found."""
  value_str = get_single_option(option_name, arguments)
  if not value_str:
    return default_value
  else:
    return int(value_str)

class Timeseries:
  """This class represents a time series, that is the result of a
     numerical evaluation of the given model. It could also be the
     gold standard to which we are comparing the results in order to
     direct the search"""
 
  # The data structure for time series is questionable, in particular
  # it would be nice to at least just use arrays.

  def __init__(self, columns, rows):
    self.columns = columns
    self.rows = rows

  def get_column_names(self):
    """Returns the set of columns, however does not include the
       first column because that is the 'Time' column"""
    return self.columns[1:]

  def get_column_data(self, column_name, start=None, end=None):
    """Retrieve the data for just a single column"""
    # We should catch the exception (I think ValueError) in case
    # the given column name is not here.
    column_index = self.columns.index(column_name)

    column_data = []
    for this_row in self.rows:
      row_time = this_row[0]
      if ( (start == None or row_time >= start) and 
           (end == None or row_time <= end) ):
        column_data.append(this_row[column_index])
    return column_data

  def write_to_file(self, results_file):
    """Format the time series and write it to the given file"""
    results_file.write("# ")
    prefix = ""
    for column in self.columns:
      results_file.write(prefix)
      results_file.write(column)
      prefix = ", "
  
    results_file.write("\n")

    for row in self.rows:
      prefix = ""
      for value in row:
        results_file.write(prefix)
        results_file.write(str(value))
        prefix = ", "
      results_file.write("\n")

  def get_best_matching_time_row(self, gold_time):
    """Return the row which has the closest time to the given target
       time. Used for comparing timeseries, we have data in one row
       of one time series and we wish to find the row of this time
       series to which to compare it."""
    best_row_distance = sys.maxint
    best_row = self.rows[0]
    # The search for the best row could be made a bit faster by
    # only starting from the previous best row
    # We could also break as soon as the value is going higher.
    for row in self.rows:
      row_distance = abs(row[0] - gold_time)
      if row_distance < best_row_distance:
        best_row = row  
        best_row_distance = row_distance
    return best_row


def parse_csv(csv, separator):
  """Parse a comma-separated value file into a timeseries"""
  header_line = csv.next()
  while not separator in header_line:
    header_line = csv.next()
  headers = [ x.lstrip().rstrip()
              for x in header_line.split(separator) ]

  rows = [] 

  try:
    while 1: 
      line = csv.next()
      if line:
        values = [ float(x.rstrip()) for x in line.split(separator) ]
        rows.append(values)
  except StopIteration:
    pass

  return Timeseries(headers, rows)
 

class MultipleCostFunctions:
  """A cost function class that does no actual costing itself but
     adds together the cost functions of a list of cost functions"""
  def __init__(self):
    self.cost_functions = []

  def add_cost_function(self, cost_function):
    """Simply add a cost function to the list of function to be
       summed together"""
    self.cost_functions.append(cost_function)
  def compare_timeseries(self, candidate_ts):
    """Compare two time series, but using each child cost function
       to compare the two time series and then summing the totals"""
    cost = 0
    for comparison in self.cost_functions:
      cost += comparison.compare_timeseries(candidate_ts)
    return cost



class FFTcost:
  """A cost function based on the fast fourier transform"""
  def __init__(self, gold_ts):
    self.gold_standard_timeseries = gold_ts
    # This clearly isn't correct, since not all of the states
    # may be periodic
    self.state_names = gold_ts.get_column_names()


  @staticmethod
  def dft_single(x, k):
    inv = -1 
    N = len(x)
    Xk = 0
    for n in xrange(N) :
      Xk += x[n] * math.e**(inv * 2j * math.pi * k * n / N)
    Xk = Xk / N
    return abs(Xk)
 
  @staticmethod
  def dft(x, kman):
    N = len(x)

    # Rather than find the energy at all frequencies, just
    # figure out the energy at a few selected points 
    k1_energy = FFTcost.dft_single(x, 1)
    kman_energy = FFTcost.dft_single(x, kman)
    khalf_energy = FFTcost.dft_single(x, N / 2)
    kquarter_energy = FFTcost.dft_single(x, N / 4)
    energies = [ k1_energy, kman_energy, 
                 khalf_energy, kquarter_energy ]
    # energies = [ FFTcost.dft_single(x, n) for n in range(1, N / 2) ]
    # The cost then is how much of the energy at all the frequencies
    # is attributable to the energy at the frequency we are interested
    # in. This is essentially saying how relatively well represented is
    # the given frequency. Note there may be a corner case here in which
    # N is less than the number of frequencies we are considering but
    # I think that is not interesting for our purposes.
    sum_energies = sum(energies)
    # Must be positive since all energies are positive and the sum of
    # all the energies must be at least as large as the single energy
    # we are considering.
    cost = (sum_energies / kman_energy) - kman_energy
    logging.debug("dft_energy       = " + str(kman_energy))
    logging.debug("total_dft_energy = " + str(sum_energies))
    return cost

  def compare_timeseries(self, candidate_ts):
    """Evaluate the given time series based on this DFT cost"""
    # So the target period is actually 
    # Should automatically check if these are in the output region
    # of the solver since otherwise we will miss up.
    start_sampling = 100.0
    stop_sampling  = 300.0
    target_period  = 20.0
    num_cycles     = (stop_sampling - start_sampling) / target_period
    # 300 / 20 = 15 
    # num_cycles = 10 # 15 
    cost = 0
    for column_name in ["mRNA"]: # self.state_names:
      column_data = candidate_ts.get_column_data(column_name)
      
      dft_cost = self.dft(column_data, num_cycles)
      logging.debug("dft_cost = " + str(dft_cost))
      cost += dft_cost 
    return cost

# Returns the indexes into candidate time series which correspond
# to the columns of the gold standard time series. The first index
# is always Time in both gold and candidate time series.
def get_candidate_indexes(gold_ts, candidate_ts):
  """ Return a list of indexes into the candidate rows
      which correspond to the headers (not including time) of the
      gold standard columns/headers. We assume that 'time' is first the
      column, so we do include that as the first index, but we don't
      actually test for it, since sometimes the column is called 'T'
      and sometimes 'time' or 'Time' or whatever"""
  candidate_indexes = [0]
  for i in range(1, len(gold_ts.columns)):
    gold_column = gold_ts.columns[i]
    for j in range(1, len(candidate_ts.columns)):
      candidate_column = candidate_ts.columns[j]
      if candidate_column == gold_column:
        candidate_indexes.append(j)
        break
  return candidate_indexes


class SpecialCost:
  """A special cost function only for the Cholesterol model.
     However in time we do wish to allow for user cost functions
     and I believe writing them in Python is not a bad way to go"""
  def __init__(self, gold_standard):
    # self.factor = len(gold_standard.rows)
    self.gold_standard = gold_standard

  def compare_timeseries(self, candidate_ts):
    """cost the time series by checking if all indexes are within
       100 and 1000, other than of course the infamous Acetyl_CoA"""
    # all_rows = candidate_ts.rows
    # last_row = all_rows[len(all_rows) - 1] 
    cost = 0

    for gold_row in self.gold_standard.rows:
      # Find the best candidate row, here we assume that time is the
      # first column
      gold_time = gold_row[0]
      # Just skip if the gold standard has a zero time point
      if gold_time <= 0.0:
        continue
      best_row = candidate_ts.get_best_matching_time_row(gold_time)

      # Now that we've found the best row, we basically ignore the
      # gold standard and just check that all values within that
      # row are within 100 and 500
      for i in range(1, len(candidate_ts.columns)):
        candidate_column = candidate_ts.columns[i]
        if candidate_column.lstrip().rstrip() != "Acetyl_CoA":
          candidate_value = best_row[i]
          if candidate_value < 100:
            diff = 500 - candidate_value
            cost += diff * diff
          if candidate_value > 1500:
            diff = candidate_value - 1500
            cost += diff * diff
    return cost 

# This should be generalised into a range cost.
# It should take the following parameters,
# 1. The columns to check, perhaps all is permitted, or as in
#    the Cholesterol model, all but a given few.
# 2. The time range, rather than look at the gold standard times
#    we should specify a range of times within which we check
#    (perhaps also we can specify an interval, to avoid checking all
#     the output).
# 3. The range values, not that this will not be per column, but a
#    global range, the intention is that this is generally used in
#    conjunction with another cost function.
class CircadCost:
  """A special cost function only for the Circadian Clock model.
     However in time we do wish to allow for user cost functions
     and I believe writing them in Python is not a bad way to go"""
  def __init__(self, gold_standard):
    # self.factor = len(gold_standard.rows)
    self.gold_standard = gold_standard

  def compare_timeseries(self, candidate_ts):
    """cost the time series by checking if all indexes are below
       100"""
    # all_rows = candidate_ts.rows
    # last_row = all_rows[len(all_rows) - 1] 
    cost = 0

    for gold_row in self.gold_standard.rows:
      # Find the best candidate row, here we assume that time is the
      # first column
      gold_time = gold_row[0]
      # Just skip if the gold standard has a zero time point
      if gold_time <= 0.0:
        continue
      best_row = candidate_ts.get_best_matching_time_row(gold_time)

      # Now that we've found the best row, we basically ignore the
      # gold standard and just check that all values within that
      # row are within 100 and 500
      for i in range(1, len(candidate_ts.columns)):
        candidate_column = candidate_ts.columns[i]
        candidate_value = best_row[i]
        if candidate_value > 100:
          diff = candidate_value - 50
          cost += diff * diff
    return cost 


class X2cost:
  """ A chi-squared cost function which works by squaring the
      differences for each state at each time point and summing
      the results."""
  def __init__(self, gold_ts):
    self.gold_standard_timeseries = gold_ts

  def compare_timeseries(self, candidate_ts):
    """Compare two time series. The left, or first arguments is taken
       to be the gold standard against which to compare. This has the
       implication that everything in the gold standard is expected to
       be in the candidate time series but not vice versa."""
    gold_ts = self.gold_standard_timeseries
    candidate_indexes = get_candidate_indexes(gold_ts, candidate_ts)

    # Now we can compare each of the rows in the gold standard to
    # the row in the candidate time series.
    cost = 0
    for gold_row in gold_ts.rows:
      # Find the best candidate row, here we assume that time is the
      # first column
      gold_time = gold_row[0]
      # Just skip if the gold standard has a zero time point
      if gold_time <= 0.0:
        continue
      best_row = candidate_ts.get_best_matching_time_row(gold_time)

      # Now that we've found the best row, we should compare the 
      # two rows. We go from '1' because we are skipping over 'time'
      for i in range(1, len(gold_row)):
        gold_value      = gold_row[i]
        c_index = candidate_indexes[i]
        candidate_value = best_row[c_index]
        diff = gold_value - candidate_value
        cost += diff * diff
    return cost


# A solver should have the following protocol illustrated with
# the sundials cvodes solver as an example:
# initialise: translate sbml model to C and compile the C program
# parameterise: output a param_overrides file which the compiled C
#               program interprets
# solve: call the getRHS script (or actually we should probably
#        call the compiled C program directly)
# I think that's all for now.
class SbmlCvodeSolver:
  """A solver which uses translation of the sbml into a
     C program which, first reads in the parameter overrides
     and then calls the Sundials-cvodes solver to obtain
     the results."""
  def __init__(self, model_file, cflags_prefix):
    self.model_file = model_file
    self.model_exec = "model.exe"
    self.param_filename = os.path.join("UserModel", "param_overrides")
    self.cflags_prefix = cflags_prefix

  # This should actually check if the model file is newer than
  # the model executable and if not then we needn't recompile it.
  def initialise_solver(self):
    """Initialise the sbml solver, for this we require the translation
       of the SBML model into a C program, and the compilation of that
       C program."""
    sbml2c_command = [ "SBML2C", self.model_file]
    sbml2c_process = Popen(sbml2c_command)
    sbml2c_process.communicate()
    if sbml2c_process.returncode != 0:
      logging.error ("SBML2C command has failed")
      sys.exit(1)
 
    base_include = os.path.join(self.cflags_prefix, "include") 
    include_flags = [ "-I" + base_include,
                      "-I" + os.path.join(base_include, "libxml2"),
                      "-I" + os.path.join(base_include, "sbsi_numeric"),
                    ]
    lib_flags = [ "-L" + os.path.join(self.cflags_prefix, "lib") ]
    extra_flags = [ "-DWL=32", "-DNO_UCF", ] 
    c_compiler = "mpic++"
    first_c_command = [ c_compiler, "-o", "main_RHS_Model.o",
                        "-c", "main_RHS_Model.C"
                      ] + include_flags + lib_flags + extra_flags
    first_c_process = Popen(first_c_command)
    first_c_process.communicate()
    if first_c_process.returncode != 0:
      logging.error ("Failed to compile main_RHS_Model.C: ")
      logging.error (" ".join(first_c_command))
      sys.exit(1)

    snd_c_command = [ c_compiler, "-o", "UserModel/UserModel.o",
                      "-c", "UserModel/UserModel.C"
                    ] + include_flags + lib_flags + extra_flags
    snd_c_process = Popen(snd_c_command)
    snd_c_process.communicate()
    if snd_c_process.returncode != 0:
      logging.error ("Failed to compile user model: ")
      logging.error (" ".join(snd_c_command))
      sys.exit(1)

    trd_c_command = [ c_compiler, "-o", self.model_exec,
                      "./main_RHS_Model.o", "UserModel/UserModel.o",
                      "-lsbsi_numeric", "-lsbml", "-lpgapack",
                      "-lfftw3", "-lsundials_kinsol",
                      "-lsundials_nvecserial", "-lsundials_cvode",
                      "-lxml2"
                    ] + lib_flags + extra_flags
    trd_c_process = Popen(trd_c_command)
    trd_c_process.communicate()
    if trd_c_process.returncode != 0:
      logging.error ("Failed to compile user model: ")
      logging.error (" ".join(trd_c_command))
      sys.exit(1)

  def parameterise_model(self, dictionary):
    """Creates a param_overrides file from the given dictionary.
       The param_overrides file can be read by the main program
       which evaluates the model (see main_RHS_Model.C) to override
       the parameter values used in the model. This is the way in
       which we can parameterise the model, rather than changing the
       SBML model file and thus incurring the cost of a C compilation"""
    param_file = open(self.param_filename, "w")
    for (name, value) in dictionary.items():
      param_file.write(name + "\t" + str(value) + "\n")
    param_file.close()

  def solve_model(self, configuration):
    """Solve the given model and return a timeseries, or None
       if solving the model has failed for some reason"""

    # We need a relatively portable way to check that this command
    # succeeds. Because we are constrained for the moment to using
    # os.system, I'm a little unsure about simply checking the return
    # code. So for now I'm removing the results file and subsequently
    # I will check if it exists.
    # We should be able to update this to use subprocess, but for now
    # this is better than the alternative of reading in the previous
    # time series.
    # We have updated this with the user of subprocess, but I'm leaving
    # it in for now since it seems quite robust, but it would certainly
    # be a slight optimistation to avoid doing this.
    results_file_prefix = "model"
    results_file = results_file_prefix + "_RHS.dat"
    if os.path.exists(results_file):
      os.remove(results_file)

    # Obviously all of these things should come from the
    # configuration
    max_times    = 100000000000
    interval     = 0.01
    out_interval = 0.1
    atol         = 1.0e-14
    reltol       = 1.0e-4

    # This could also be run without mpi.
    mpirun_command = [ "mpirun",
                       self.model_exec,
                       "model_name", # Could get model name from xml file
                       str(configuration.t_final),
                       str(configuration.t_init),
                       str(max_times),
                       str(interval),
                       str(out_interval),
                       str(atol),
                       str(reltol),
                       results_file_prefix,
                     ]
    mpi_process = Popen(mpirun_command, stdout=PIPE)
    mpi_output = mpi_process.communicate()[0]
    logging.debug(mpi_output)
                    
    # So we check if the results file actually exists and if
    # not we assume it failed. Also now I can actually check
    # the return code
    if not os.path.exists(results_file) or mpi_process.returncode != 0:
      logging.warning("Model solving failed")
      return None 

    csv_file = open(results_file,  "r")
    timeseries = parse_csv(csv_file, "\t")
    csv_file.close()
    return timeseries

class BioPEPASolver:
  """A solver which has as input Bio-PEPA files and simply uses
     the command-line version of the Bio-PEPA eclipse plugin to
     solve the model"""
  def __init__(self, model_file, biopepajar):
    self.model_file = model_file
    (file_directory, basename) = os.path.split(model_file)
    self.paramed_file = os.path.join(file_directory, "p_" + basename)
    self.biopepajar = biopepajar

  # This should actually check if the model file is newer than
  # the model executable and if not then we needn't recompile it.
  def initialise_solver(self):
    """Initialise the Bio-PEPA solver, nothing special required here"""
    pass

  def parameterise_model(self, dictionary):
    """A simple and obviously broken method for parameterising a
       Bio-PEPA file. We assume that each parameter definition begins
       a line, with 'param_name = ....' and we simply replace the dots
       with the new value of the parameter."""
    biopepa_file = open(self.model_file, "r")
    output_file = open(self.paramed_file, "w")

    parameterised_names = []

    for line in biopepa_file:
      if '=' in line:
        name = line.partition('=')[0]
        name = name.lstrip().rstrip()
        if dictionary.has_key (name):
          new_value = dictionary[name]
          output_file.write(name + " = " + str(new_value) + " ;\n")
          parameterised_names.append(name)
        else:
          output_file.write(line)
      else:
        output_file.write(line)

    unparameterised = [ x for x in dictionary.keys()
                            if x not in parameterised_names ]
    if unparameterised:
      logging.error ("Failed to fully parameterise Bio-PEPA file")
      sys.exit(1)

    biopepa_file.close()
    output_file.close()

  def solve_model(self, configuration):
    """Solve the parameterised version of the model. This assumes
       that parameterise_model has already been called and that hence
       a parameterised version of the model exists in a file named
       self.paramed_file"""

    biopepa_command = [ "java",
                        "-jar",
                        "biopepa.jar",
                        "timeseries",
                        self.paramed_file,
                        "--no-warnings", 
                        "--solver",
                        "dopr-adaptive",
                        "--timeStep",
                        "0.01",
                        "--startTime",
                        str(configuration.t_init),
                        "--stopTime",
                        str(configuration.t_final),
                        "--dataPoints",
                        "11"
                        # "--output-file",
                        # csvfile 
                      ]
    # print (biopepa_command)
    biopepa_process = Popen(biopepa_command, stdout=PIPE)
    output = biopepa_process.communicate()[0]

    if biopepa_process.returncode != 0:
      logging.error ("biopepa process failed to return")
      sys.exit(1)
    output_lines = output.split("\n")
    timeseries = parse_csv(output_lines.__iter__(), ", ")

    return timeseries
 
def evaluate_individual(individual, configuration):
  """Evaluate an individual by first parameterising the model
     and then numerically evaluating it to produce a time series.
     Finally we run the configuration's cost function over the
     produced results"""
  solver = configuration.solver
  solver.parameterise_model(individual.dictionary)
  timeseries = solver.solve_model(configuration)
  individual.results = timeseries

  # So if solving the model didn't actually work then return
  # the maximum integer as the cost, this will just mean that this
  # parameterisation will be ignore with respect to the search for
  # good values. There isn't an awful lot else that we can do.
  if timeseries == None:
    return sys.maxint

  # If we really did get a timeseries then use the configuration's
  # cost function to obtain the cost for this set of parameters
  cost_function = configuration.cost_function
  cost = cost_function.compare_timeseries(timeseries)
  logging.info("Individual " + str(individual.number) +
               " Cost: " + str(cost))
  return cost

 

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
    if self.should_mutate():
      return random.uniform(self.low, self.high)
    else:
      return self.default_value
    
class Individual:
  """An individual is essentially a dictionary mapping the
     optimisable parameters to their chosen values. We also
     store the results of numerical analysis"""
  def __init__(self, number, dictionary):
    self.number = number
    self.dictionary = dictionary
    # These should be filled in when the individual is evaluated.
    self.results = None

def create_default_individual(number, parameters):
  """Create a default individual with all parameters set to their
     default values"""
  individual_dictionary = dict()
  for param in parameters:
    individual_dictionary[param.name] = param.default_value
  return Individual(number, individual_dictionary)


class Optimisation:
  """A class holding a single optimisation problem, essentially
     this just stores the data we have gained from the files for
     the model, parameters and the gold standard"""
  def __init__(self, model_file, gold_standard, parameters):
    self.model_file = model_file
    self.gold_standard = gold_standard
    self.parameters = parameters
  def initialise(self):
    """Initialise the optimisation process"""
    pass
  


def update_default_parameters(parameters, best_citizen):
  """When we find a new individual that has the best cost,
     we generally wish to update the parameters' default values to
     reflect the new best fit. This affects the way in which new
     values will be chosen"""
  for param in parameters:
    param.default_value = best_citizen.dictionary[param.name]
 
class SimulatedAnnealing:
  """A class defining the simulated annealing algorithm for
     optimisation"""
  def __init__(self):
    self.algorithm_name = "sa"
    self.evaluated = 0

  @staticmethod 
  def move_state_function(current_cost, cost, temperature):
    """ The temperature name for the argument refers to the simulated
        annealing origin in that of metallurgy. To be precise,
        the idea is that for large temperatures the evolution of the
        state of the system is sensitive to coarseer energy variations,
        while it is sensitive to finer energy variations
        when T is small"""
    # So if the current cost is not as good as the one being
    # evaluated then definitely move to the new individual
    if current_cost > cost:
      return True
    # If not, we may still wish to move 'uphill' in order to avoid
    # being stuck in a local minimum.
    move_probability = math.exp ((current_cost - cost) / temperature)
    coin = random.random()
    logging.debug("Uphill moves: ")
    logging.debug("    current_cost = " + str(current_cost))
    logging.debug("    cost = " + str(cost))
    logging.debug("    move_probability = " + str(move_probability))
    logging.debug("    coin = " + str(coin))
    # We should move if the coin is below the probability of moving
    return coin < move_probability

  @staticmethod
  def get_number_to_change(num_params):
    """Decide how many parameters to change, at least one and
       may be more with decreasing probability"""
    # First of all we will modify at least one.
    # After that we continually have a
    # 50/50 chance of adding one more parameter to change
    num_to_change = 1
    coin = random.choice([True, False])
    while coin and num_to_change < num_params:
      # note the inequality above will not allow us to change ALL
      # the parameters (unless there is only one) but I think that
      # that seems pretty reasonable.
      num_to_change += 1
      coin = random.choice([True, False])
    return num_to_change
  
  @staticmethod
  def mutate_parameter(change_param, current_value, temperature):
    """Given a parameter to change, decide on its new value."""
    # The temperature is evals/maxevals, so how far through
    # the maximum search range we are. Hence we want 1 minus that
    # to give us the proportion of the range that we should choose
    # from. 
    range_size = (1 - temperature) * (change_param.get_value_range_size())
    half_range = range_size / 2
    # Instead of just using max and min we could try to shift the
    # range of possible values so that we don't go over the highest
    # or under the lowest.
    highest = min(current_value + half_range, change_param.high)
    lowest = max(current_value - half_range, change_param.low)
    new_value = random.uniform(lowest, highest)
    return new_value

  def get_neighbour(self, parameters, 
                    current_individual,
                    num_evals,
                    temperature):
    """Based on the current state that we're in (and state here is
       simply an individual, ie. a chosen set of values for the
       optimisable parameters), find a new state (ie. a new individual)
       sufficiently close to the current individual but also different
       such that we can evaluate that and decide both whether it is
       a better fit and whether we should move to that state"""
    current_dictionary = current_individual.dictionary
    individual_dictionary = dict()

    # essentially copy over the entire current dictionary
    for param in parameters:
      individual_dictionary[param.name] = current_dictionary[param.name]

    # Decide how many parameters we should modify,
    num_to_change = self.get_number_to_change(len(parameters))

    # Now we know how many parameters to change so we simply
    # select a random sample of that many parameters and change
    # them.
    change_parameters = random.sample(parameters, num_to_change)
    for change_param in change_parameters:
      # Note that this does not make reference to the parameter's
      # default value, because the default value is that of the
      # current best citizen, which is reflected in current_value
      # and NOT within the parameter's original definition
      current_value = individual_dictionary[change_param.name]
      new_value = self.mutate_parameter(change_param,
                                        current_value, 
                                        temperature)
   
      individual_dictionary[change_param.name] = new_value
    return Individual(num_evals + 1, individual_dictionary)

  def run_optimisation(self, optimisation, configuration):
    """Run the optimisation algorithm using simulated annealing"""
    parameters = optimisation.parameters
    # The initial state is the default individual
    current_individual = create_default_individual(1, parameters)
    current_cost = evaluate_individual(current_individual, 
                                       configuration)
    # The initial bests then are the default:
    best_individual = current_individual
    best_cost = current_cost
    # How many evaluations we have done
    num_evals = 1
    max_evals = configuration.num_generations
    target_cost = configuration.target_cost
  
    while num_evals < max_evals and target_cost < best_cost:
      temperature = float(num_evals)/float(max_evals)
      new_individual = self.get_neighbour(parameters,
                                          current_individual,
                                          num_evals,
                                          temperature)   
      new_cost = evaluate_individual(new_individual,
                                     configuration)
      num_evals += 1
      # Should we select this new individual?
      if self.move_state_function(current_cost, new_cost, temperature):
        current_individual = new_individual
        current_cost = new_cost
      # If this new individual is the better than the best so far
      # then we save it as the best so far.
      if new_cost < best_cost:
        best_individual = new_individual
        best_cost = new_cost

    logging.info("Best cost = " + str(best_cost))
    return best_individual
    
    
class SimplestSearch:
  """Implements the optimisation process as the simplest kind of
     genetic algorithm search. For each generation we generate
     a population size amount of new candidates. Each candidate is
     a random new individual, taking the default value for each
     parameter and then with that parameter's mutation probability
     either leaving it as it is or mutating it to be a completely
     random value between that parameter's highest and lowest values.
     We then evaulate each individual and select the best, if it is
     better than the current best individual we update and then run
     a new generation and keep repeating until we have run the
     configured number of generations"""

  def __init__(self):
    self.algorithm_name = "simple"

  @staticmethod 
  def create_individual(number, parameters):
    """Create an individual with each optimisable parameter
       either kept as the current default value or given a new
       random number within the full range of that parameter"""
    individual_dictionary = dict()
    changed = 0
    for param in parameters:
      if param.should_mutate():
        new_value = param.get_random_value_full_range()
        individual_dictionary[param.name] = new_value
        changed += 1
      else:
        individual_dictionary[param.name] = param.default_value
    # After the for loop we make sure we have changed at least one
    # parameter, otherwise we choose one at random to change.
    if changed == 0:
      param = random.choice(parameters)
      new_value = param.get_random_value_full_range()
      individual_dictionary[param.name] = new_value
    return Individual(number, individual_dictionary)


  def run_generation(self, gen_number, optimisation, configuration):     
    """Run a single generation. We simply generate a population
       size number of random individuals and evaluate them all
       independently of each other. We then evaluate all the
       individuals and select the best. We return the best citizen
       and its cost"""
    generation = []
    for i in range (0, configuration.population_size):
      individual = self.create_individual(i, optimisation.parameters)
      generation.append(individual)

    lowest_cost = sys.maxint
    best_citizen = generation[0]
    for citizen in generation:
      cost = evaluate_individual(citizen,
                                 configuration)
      if cost < lowest_cost:
        lowest_cost = cost
        best_citizen = citizen

    logging.info ("Generation (" + str(gen_number) + 
                  ") has lowest cost: " + str(lowest_cost))
    logging.info ("=============================")
    return (lowest_cost, best_citizen)

  def run_optimisation(self, optimisation, configuration):
    """Run the optimisation process. Successively run generations,
       keeping a track of the best individual (parameter settings)
       encountered so far. Whenever we encounter a better individual
       we update all the parameters such that the default values are
       equal to the best found so far"""
    (lowest_cost, best_citizen) = self.run_generation(1, optimisation,
                                                      configuration)

    for i in range (2, configuration.num_generations + 1):
      # Modify the default values of all the parameters to be the
      # current best.
      for param in optimisation.parameters:
        param.default_value = best_citizen.dictionary[param.name]
      logging.info ("Beginning generation: " + str(i))
      (this_cost, this_mayor) = self.run_generation(i, optimisation,
                                                    configuration)
      if this_cost < lowest_cost:
        lowest_cost = this_cost
        best_citizen = this_mayor
    return best_citizen



def get_gold_standard_timeseries(filename):
  """Open the gold standard file and parse in as a comma-separated value
     file. Obtaining a time series which is returned"""
  gold_file = open(filename,  "r")
  gold_standard = parse_csv(gold_file, ", ")
  gold_file.close()
  return gold_standard

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
  
  return parameters


class Configuration:
  """A class to store the configuration in"""
  def __init__(self, optimisation):
    self.optimisation = optimisation
    self.num_generations = 5
    self.population_size = 5
    self.t_final = 1.0
    self.t_init  = 0.0

    home_dir = "/afs/inf.ed.ac.uk/user/a/aclark6/" 
    cflags_prefix = os.path.join (home_dir, "Source/svn-git-sbsi/install/")
    model_file = optimisation.model_file
    self.solver = SbmlCvodeSolver(model_file, cflags_prefix)

    self.search_algorithm = SimulatedAnnealing()
    self.target_cost = 0
    self.cost_function = None


  def set_solver(self, solver_name):
    """Set the solver of the configuration"""
    if solver_name == "biopepa":
      model_file = self.optimisation.model_file
      biopepajar = "biopepa.jar"
      self.solver = BioPEPASolver(model_file, biopepajar)
    else:
      print("Unrecognised solver name: " + solver_name)
      print("Please choose from: biopepa")
      print("Or nothing for the default sbml cvode solver")
      sys.exit(1)

  def set_search_agorithm(self, algorithm):
    """Set the search algorithm of the configuration"""
    if (self.search_algorithm and
        self.search_algorithm.algorithm_name == algorithm):
      pass
    elif algorithm == "simple":
      self.search_algorithm = SimplestSearch()
    elif algorithm == "sa" :
      self.search_algorithm = SimulatedAnnealing()
    # elif pga, particleswarm
    else:
      print ("Unrecognised algorithm name: " + algorithm)
      print ("Choose from 'sa' or 'simple'")
      sys.exit(1)

  @staticmethod
  def cost_function_of_name(function_name, gold_standard):
    """Return a cost function corresponding to the given name"""
    if function_name == "x2":
      return X2cost(gold_standard)
    elif function_name == "fft":
      return FFTcost(gold_standard)
    elif function_name == "special":
      return SpecialCost(gold_standard)
    elif function_name == "circad":
      return CircadCost(gold_standard)
    else:
      print ("Unrecognised cost function name: " + function_name)
      print ("Choose from: x2 fft")
      sys.exit(1)

  def set_cost_function(self, function_names, gold_standard):
    """Set the cost function depending on the given list
       of function names"""
    if not function_names:
      # The default of a chi squared cost function
      self.cost_function = X2cost(gold_standard)
    elif len(function_names) == 1:
      self.cost_function = self.cost_function_of_name(function_names[0],
                                                      gold_standard)
    else:
      multiple_costs = MultipleCostFunctions()
      for fname in function_names:
        cost_function = self.cost_function_of_name(fname, gold_standard)
        multiple_costs.add_cost_function(cost_function)
      self.cost_function = multiple_costs

def get_configuration(arguments, optimisation):
  """Return a configuration based on the command line arguments"""
  configuration = Configuration(optimisation)
  generations = get_single_option("generations", arguments)
  if generations:
    configuration.num_generations = int(generations)
  population_size = get_single_option("population", arguments)
  if population_size:
    configuration.population_size = int(population_size)
  stop_time = get_single_option("stop-time", arguments)
  if stop_time:
    configuration.t_final = float(stop_time)
  algorithm = get_single_option("algorithm", arguments)
  if algorithm:
    configuration.set_search_agorithm(algorithm)
  target_cost = get_single_option("target-cost", arguments)
  if target_cost:
    configuration.target_cost = float(target_cost)
  cost_functions = get_options("cost-function", arguments)
  configuration.set_cost_function(cost_functions,
                                  optimisation.gold_standard)
  solver = get_single_option("solver", arguments)
  if solver:
    configuration.set_solver(solver)
 
  return configuration

def initialise_logger(arguments):
  """Initialise the logging system, depending on the arguments
     which may set the log level and a log file"""
  log_level = get_single_option("loglevel", arguments)
  if log_level:
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
      print ("The log level must be one of the following:")
      print ("    debug, info, warning, error, critical")
      print ("Exiting")
      sys.exit(1) 
  else:
    numeric_level = logging.INFO

  # We could also change the format of the logging messages to
  # something like: format='%(levelname)s:%(message)s'

  log_file = get_single_option("logfile", arguments)
  if log_file:
    logging.basicConfig(filename=log_file, level=numeric_level)
  else:
    logging.basicConfig(level=numeric_level)



def get_optimisation_definition(filenames):
  """Return the definition of the optimisation based on the
     filenames given in the command line arguments"""
  model_file = filenames[0]
  gold_standard_file = filenames[1]
  init_param_file = filenames[2]
 
  gold_standard = get_gold_standard_timeseries(gold_standard_file)
  parameters    = get_init_param_parameters(init_param_file) 
  
  optimisation = Optimisation(model_file, gold_standard, parameters)
  optimisation.initialise()
  return optimisation



def run():
  """Process all the command line arguments and get going with the
     optimisation""" 
  # The command line arguments not including this script itself
  arguments    = [ x for x in sys.argv if not x.endswith(".py") ]
  filenames    = [ x for x in arguments if not '=' in x ]

  # I really just want to take in one file name which is a configuration
  # file which contains the others, but for now I'll do this:
  if len(filenames) < 3:
    print ("You must provide three files:")
    print ("   the model file")
    print ("   the gold standard data file")
    print ("   the initial parameter file")
    sys.exit(1)

  optimisation = get_optimisation_definition(filenames)

  configuration = get_configuration(arguments, optimisation)
  initialise_logger(arguments)
  configuration.solver.initialise_solver()

  algorithm = configuration.search_algorithm
  best_citizen = algorithm.run_optimisation(optimisation, configuration)
  
  print("After: " + str(configuration.num_generations) + " generations:")
  
  best_parameters_fname = "best_params"
  best_parameters_file = open(best_parameters_fname, "w") 
  for param in optimisation.parameters:
    line = param.name + ": " + str(best_citizen.dictionary[param.name])
    print (line)
    best_parameters_file.write(line)
    best_parameters_file.write("\n")

  best_parameters_file.close()

  if best_citizen.results:
    best_timeseries_fname = "best_timeseries.csv"
    best_timeseries_file = open(best_timeseries_fname, "w")
    best_timeseries = best_citizen.results
    best_timeseries.write_to_file(best_timeseries_file)
    best_timeseries_file.close()
    print ("Best timeseries written to file: " + best_timeseries_fname)
  else:
    print ("No best time series")

if __name__ == "__main__":
  random.seed()
  run()



