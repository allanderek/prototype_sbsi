"""
This is a prototype of a new version of the sbsi optimisation
framework. The idea being that we are less coupled within the code
such that we combining together many smaller applications.
"""
import sys
import argparse
import os.path
import random
import math
import logging
import shutil

# import configuration
import solve_model
import timeseries
import parameters
import range_cost
import fft_cost
import utils

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
    else:
      print ("gold column: " + gold_column + "not found in candidate ts")
      raise StandardError
  return candidate_indexes

class SpecialCost(range_cost.RangeCost):
  """A special cost function only for the Cholesterol model.
     However in time we do wish to allow for user cost functions
     and I believe writing them in Python is not a bad way to go"""
  def __init__(self, gold_standard):
    super(SpecialCost, self).__init__()
    self.set_ignored_columns("Acetyl_CoA")
    self.set_lower_limit(100)
    self.set_upper_limit(1500)

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
      assert (len(gold_row) <= len(candidate_indexes))
      for i in range(1, len(gold_row)):
        gold_value      = gold_row[i]
        c_index = candidate_indexes[i]
        candidate_value = best_row[c_index]
        diff = gold_value - candidate_value
        cost += diff * diff
    number_of_data_points = gold_ts.number_of_data_points()
    normalised = cost / number_of_data_points
    return normalised


class Results:
  """A class to hold a single set of results associated with a given
     gold standard and model file"""
  def __init__(self, model_data, timecourse):
    self.timecourse = timecourse
    self.model_data = model_data

  def results_filename (self):
    """Return a suitable filename for the results contained here in"""
    model_file = self.model_data.model_file
    basename = os.path.basename(model_file)
    result_filename = utils.change_filename_ext(basename, ".csv")
    return result_filename

  def write_to_file(self, csvfile):
    """Writes the time course results to the given csv file,
       note that the file should already be open for writing"""
    if self.timecourse:
      self.timecourse.write_to_file(csvfile)

# It still seems strange that this method is orphanless in the 
# middle of the optimiser, perhaps it should go in ModelData?
def evaluate_individual(individual, configuration):
  """Evaluate an individual by first parameterising the model
     and then numerically evaluating it to produce a time series.
     Finally we run the configuration's cost function over the
     produced results"""
  optimisation = configuration.optimisation
  individual.results = []
  total_cost = 0
  for model_data in optimisation.model_datas:
    solver = model_data.solver
    solver.parameterise_model(individual.dictionary)
    timecourse = solver.solve_model(model_data.solver_config)
    result = Results(model_data, timecourse)
    individual.results.append(result)

    # So if solving the model didn't actually work then return
    # the maximum integer as the cost, this will just mean that this
    # parameterisation will be ignore with respect to the search for
    # good values. There isn't an awful lot else that we can do.
    if timecourse == None:
      configuration.monitor.increase_failed_solves()
      return sys.maxint
    else:
      configuration.monitor.increase_success_solves()

    # If we really did get a timeseries then use the configuration's
    # cost function to obtain the cost for this set of parameters
    cost_function = model_data.cost_function
    cost = cost_function.compare_timeseries(timecourse)
    total_cost += cost

  logging.info("Individual " + str(individual.number) +
               " Cost: " + str(total_cost))
  return cost

   
class Individual:
  """An individual is essentially a dictionary mapping the
     optimisable parameters to their chosen values. We also
     store the results of numerical analysis"""
  def __init__(self, number, dictionary):
    self.number = number
    self.dictionary = dictionary
    # These should be filled in when the individual is evaluated.
    self.results = []

  def defined_parameter_names(self):
    """Return the parameter names held within the
       individual's dictionary"""
    return self.dictionary.keys()

  def get_param_value(self, name):
    """Returns the value of a parameter within this individual"""
    return self.dictionary[name] 

def create_default_individual(number, params):
  """Create a default individual with all parameters set to their
     default values"""
  individual_dictionary = dict()
  for param in params:
    individual_dictionary[param.name] = param.default_value
  return Individual(number, individual_dictionary)


def update_default_parameters(params, best_citizen):
  """When we find a new individual that has the best cost,
     we generally wish to update the parameters' default values to
     reflect the new best fit. This affects the way in which new
     values will be chosen"""
  for param in params:
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

  def get_neighbour(self, params, 
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
    for param in params:
      individual_dictionary[param.name] = current_dictionary[param.name]

    # Decide how many parameters we should modify,
    num_to_change = self.get_number_to_change(len(params))

    # Now we know how many parameters to change so we simply
    # select a random sample of that many parameters and change
    # them.
    change_parameters = random.sample(params, num_to_change)
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
    params = optimisation.parameters
    # The initial state is the default individual
    current_individual = create_default_individual(1, params)
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
      new_individual = self.get_neighbour(params,
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
    return OptimiserResults(best_individual, best_cost)
    
    
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
  def create_individual(number, params):
    """Create an individual with each optimisable parameter
       either kept as the current default value or given a new
       random number within the full range of that parameter"""
    individual_dictionary = dict()
    changed = 0
    logging.debug("Creating Individual: " + str(number))
    for param in params:
      if param.should_mutate():
        new_value = param.get_random_value_full_range()
        individual_dictionary[param.name] = new_value
        logging.debug(param.name + " = " + str(new_value))
        changed += 1
      else:
        individual_dictionary[param.name] = param.default_value
    # After the for loop we make sure we have changed at least one
    # parameter, otherwise we choose one at random to change.
    if changed == 0:
      logging.debug("Zero would have been changed, forcing choice")
      param = random.choice(params)
      new_value = param.get_random_value_full_range()
      individual_dictionary[param.name] = new_value
      logging.debug(param.name + " = " + str(new_value))
    return Individual(number, individual_dictionary)


  def run_generation(self, gen_number, optimisation, configuration):     
    """Run a single generation. We simply generate a population
       size number of random individuals and evaluate them all
       independently of each other. We then evaluate all the
       individuals and select the best. We return the best citizen
       and its cost"""
    generation = []
    # Now if the user has specified the population size then
    # we respect that, otherwise we take a default equal to the
    # number of parameters we are trying to optimise for.
    if configuration.population_size != None:
      population_size = configuration.population_size
    else:
      population_size = len(optimisation.parameters)
    for i in range (0, population_size):
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

    record_frequency = configuration.record_frequency
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

      if record_frequency and (i % record_frequency) == 0:
        optim_results = OptimiserResults (best_citizen, lowest_cost)
        results_dir = "ckpoint_" + str(i)
        optim_results.report_results(configuration,  results_dir)
    return OptimiserResults (best_citizen, lowest_cost)


class StructuredExperiment:
  """Instead of running an optimisation in order to attempt to find
     the 'best' citizen we instead run a grid of experiments and
     report all the results. (We could also return the best_citizen).
  """
  def __init__(self):
    self.algorithm_name = "experiment"
    # We keep a count of how many individuals we have evaluated.
    self.num_individuals = 0
    # The results database will map from 'individuals' to results
    # where results is an evaluation, that is a score/cost.
    self.results_database = dict()

  # The standard library insists that 'range' arguments are integers.
  # This somewhat makes sense because of the rounding issues with floats.
  # Should range(0.0,2.2,1.1) return [0.0,1.1] or [0.0,1.1,2.199999999]?
  # There's no way to be certain. But since we really want this behaviour
  # we just write our own range function.
  @staticmethod
  def float_range(start, stop, step):
    """The same as the standard library 'range' function except that
       it works over floating point numbers rather than simply integers.
    """
    value = start
    while True:
      if value >= stop:
        return
      yield value
      value += step

  def range_parameter(self, ind_dict, params, configuration):
    """ A recursive function which must range over
        the curent parameter's range, plus recursively range over all
        other parameters. This function actually returns the evaluated
        model results via placing them in this object's results database
    """
    # Nothing to do if the list of parameters is empty
    if not params:
      return
    param = params[0]
    other_params = params[1:]
    values = self.float_range(param.low, param.high, param.step_size)
    if other_params:
      for value in values:
        # Note that I don't actually need to copy the dictionary
        # I only need to do that for the leaf nodes.
        this_ind_dict = ind_dict # ind_dict.copy()
        this_ind_dict[param.name] = value
        self.range_parameter(this_ind_dict, other_params, configuration)
    else:
      for value in values:
        this_ind_dict = ind_dict.copy()
        this_ind_dict[param.name] = value
        self.num_individuals += 1
        individual = Individual(self.num_individuals, this_ind_dict)
        cost = evaluate_individual(individual, configuration)
        # This has a bit of a bad code smell to it, I'm comparing
        # against a value used to indicate that no timeseries was
        # obtained (probably because the solver failed)
        if cost != sys.maxint:
          self.results_database[individual] = cost

      

  def run_optimisation(self, optimisation, configuration):
    """The main entry point of the StructuredExperiment algorithm. """
    results_filename = "results_db.csv"
    self.range_parameter(dict(), optimisation.parameters, configuration)
    results_file = open (results_filename,  "w")
    write_results_database_file(self.results_database, results_file)
    results_file.close()
    return ExperimentResults(results_filename)

def write_results_database_file(results_db, results_file):
  """Format the results database and write to file"""
  # This depends on all of the individuals' having the same parameter
  # names in their dictionaries.
  result_items = results_db.items()
  first_ind = result_items[0][0]

  first_columns = first_ind.defined_parameter_names()
  results_file.write("# cost") 
  for name in first_columns:
    results_file.write(", " + name)
  results_file.write("\n")
  for individual, cost in result_items:
    results_file.write(str(cost))
    for name in first_columns:
      results_file.write(", " + str(individual.get_param_value(name)))
    results_file.write("\n")
 

def get_gold_standard_timeseries(filename):
  """Open the gold standard file and parse in as a comma-separated value
     file. Obtaining a time series which is returned"""
  return timeseries.get_timecourse_from_file(filename)

class Monitor:
  """ A class built to monitor or more rather log the progress of
      the optimisation routine. This is more than 'logging' in that
      we are recording actual data, such as the success or otherwise
      of the simulation runs."""
  def __init__(self):
    self.successful_solves = 0
    self.failed_solves = 0

  def increase_success_solves(self, amount=1):
    """Increase by the given amount (default 1) the number of
       successfully solved simulation runs"""
    self.successful_solves += amount

  def increase_failed_solves(self, amount=1):
    """Decrease by the given amount (default 1) the number of
       failed simulation runs""" 
    self.failed_solves += amount

  def get_successful_solves(self):
    """Return the sofar recorded number of successful simulation runs"""
    return self.successful_solves
  def get_failed_solves(self):
    """Return the sofar recorded number of failed simulation runs"""
    return self.failed_solves


class ModelData:
  """A class holding the experimental data for the given model.
     Essentially this just olds a model file and a gold standard file"""
  def __init__(self, model_file, gold_standard_file, arguments):
    self.model_file = model_file
    self.gold_standard_file = gold_standard_file
    self.gold_standard = get_gold_standard_timeseries(gold_standard_file)
    self.solver_config = arguments
    self.solver = solve_model.get_solver(model_file, arguments)
    self.cost_function = None

  def initialise_solver(self):
    """Initialise the solver associated with the model data"""
    self.solver.initialise_solver()

  @staticmethod
  def cost_function_of_name(function_name, gold_standard):
    """Return a cost function corresponding to the given name"""
    if function_name == "x2":
      return X2cost(gold_standard)
    elif function_name == "fft":
      return fft_cost.FFTcost(gold_standard)
    elif function_name == "special":
      return SpecialCost(gold_standard)
    elif function_name == "range":
      return range_cost.RangeCost(gold_standard)
    else:
      print ("Unrecognised cost function name: " + function_name)
      print ("Choose from: x2 fft")
      sys.exit(1)

  def set_cost_function(self, function_names):
    """Set the cost function depending on the given list
       of function names"""
    if not function_names:
      # The default of a chi squared cost function
      self.cost_function = X2cost(self.gold_standard)
    elif len(function_names) == 1:
      self.cost_function = self.cost_function_of_name(function_names[0],
                                                      self.gold_standard)
    else:
      multiple_costs = MultipleCostFunctions()
      for fname in function_names:
        cost_function = self.cost_function_of_name(fname,
                                                   self.gold_standard)
        multiple_costs.add_cost_function(cost_function)
      self.cost_function = multiple_costs



class Optimisation:
  """A class holding a single optimisation problem, essentially
     this just stores the data we have gained from the files for
     the model, parameters and the gold standard"""
  def __init__(self, params_file, model_datas):
    self.params_file = params_file
    self.parameters = parameters.get_init_param_parameters(params_file)
    self.model_datas = model_datas

  def get_final_gold_standard_time(self):
    """Returns the final time of any gold standard data file within
       this optimisation. This is generally useful for deciding how
       long the model should be solved for (how long in simulation time)
    """
    final_time = 0.0
    for model_data in self.model_datas:
      model_final_time = model_data.gold_standard.get_final_time()
      final_time = max(model_final_time, final_time)
    return final_time


  def initialise(self):
    """Initialise the optimisation process"""
    pass

  def initialise_solvers(self):
    """Initialise all of the solvers associated with each of the
       model datas we have"""
    for model_data in self.model_datas:
      model_data.initialise_solver()
 
class Configuration:
  """A class to store the configuration in"""
  def __init__(self, arguments, optimisation):
    self.optimisation = optimisation
    self.num_generations = 5
    self.population_size = None

    if not arguments.stop_time:
      arguments.stop_time = optimisation.get_final_gold_standard_time()

    self.search_algorithm = SimplestSearch()
    self.target_cost = 0
    self.cost_function = None
    self.monitor = Monitor()
    self.record_frequency = None

  def report_on_solves(self):
    """Report to the log, information about the simulation runs,
       in particular the number which failed and succeeded."""
    successes = self.monitor.get_successful_solves()
    failures = self.monitor.get_failed_solves()
    ratio = float (successes + failures)
    success_percentage = int(float(100 * successes) / ratio)
    logging.info("Successful solves: " + str(successes))
    logging.info("Failed solves: " + str(failures))
    logging.info("Success percentage: " + str(success_percentage))

  def report_on_best_params(self, best_params):
    """Report to the log information about the best parameters found.
       Essentially reporting on whether or not any of the parameters
       seem too close to their range limits or if any of them have
       not been changed from their initial values
    """
    class Arguments:
      """A dummy class just to provide the required arguments to
         the parameters.check_parameters method.
      """
      def __init__(self):
        self.tolerance = 0.01
    arguments = Arguments()
    params = self.optimisation.parameters
    failed_results = parameters.check_parameters(params, 
                                                 best_params, 
                                                 arguments)
    # This is bad, I've basically copy and pasted this from parameters.py
    # what I should do is have a 'report_on_parameters' in parameters.py
    # which also using the logging module.
    if len(failed_results) == 1:
      logging.info ("No parameters were too close to their range limits")
    else:
      for fail_result in failed_results:
        param = fail_result.param
        if not fail_result.value:
          logging.warning (param.name + 
                           " is not in the best params file(s)")
        elif fail_result.unchanged:
          logging.warning (param.name + 
                           " has not changed from the initial value")
        else:  
          logging.warning (param.name + 
                           " is too close to its range limits")
          logging.warning ("   high limit: " + str(param.high))
          logging.warning ("   low  limit: " + str(param.low))
          logging.warning ("   value     : " + str(fail_result.value))
          logging.warning ("   too high  : " + str(fail_result.too_high))
          logging.warning ("   too low   : " + str(fail_result.too_low))

  def set_search_agorithm(self, algorithm):
    """Set the search algorithm of the configuration"""
    if (self.search_algorithm and
        self.search_algorithm.algorithm_name == algorithm):
      pass
    elif algorithm == "simple":
      self.search_algorithm = SimplestSearch()
    elif algorithm == "sa" :
      self.search_algorithm = SimulatedAnnealing()
    elif algorithm == "experiment":
      self.search_algorithm = StructuredExperiment()
    # elif pga, particleswarm
    else:
      print ("Unrecognised algorithm name: " + algorithm)
      print ("Choose from 'sa' or 'simple'")
      sys.exit(1)

def get_configuration(arguments, optimisation):
  """Return a configuration based on the command line arguments"""
  configuration = Configuration(arguments, optimisation)

  configuration.num_generations = arguments.generations
  configuration.population_size = arguments.population
  configuration.set_search_agorithm(arguments.algorithm)
  configuration.target_cost = float(arguments.target_cost)
  if arguments.record_freq and arguments.record_freq != 0:
    configuration.record_frequency = arguments.record_freq
 
  return configuration


def get_optimisation_definition(arguments):
  """Return the definition of the optimisation based on the
     filenames given in the command line arguments"""
  filenames = arguments.filenames
  params_file = filenames[0]
  model_datas = [] 
  # -1 for the params file, and step two because we will get
  # two at a time. Note that we should check if the number of
  # file names is odd, it should be since it should be the
  # params file, plus a list of pairs. 
  for index in range(1, len(filenames) - 1, 2):
    model_file = filenames[index]
    gold_file = filenames[index + 1]
    model_data = ModelData(model_file, gold_file, arguments)
    model_data.set_cost_function(arguments.cost_function)
    model_datas.append(model_data)
  
  optimisation = Optimisation(params_file, model_datas)
  optimisation.initialise()
  return optimisation


def create_arguments_parser(add_help):
  """Create the parser for the command-line arguments"""
  # Create the solve model parser which will give us some arguments
  # controlling the solver (start/stop_time etc)
  solve_model_parser = solve_model.create_arguments_parser(False)

  description = "Perform an optimisation for parameter values"
  parser = argparse.ArgumentParser(add_help=add_help,
                                   parents=[solve_model_parser],
                                   description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
    help="The input files: model gold_standard initparams")
  cost_function_choices = ["x2", "fft", "special", "circad" ]
  parser.add_argument('--cost-function', action='append',
                      choices=cost_function_choices, 
    help="Set the cost function(s) to be used to evaluate individuals")
  parser.add_argument('--target-cost', action='store',
                      type=int, default=0,
    help="Set the target cost to meet, default = 0")
  parser.add_argument('--generations', action='store',
                      type=int, default=10,
    help="Set the number of generations to evolve")
  # There is no default for the population size, this is such that
  # the algorithms may distinguish between a population size set by
  # the user and the default. Hence the each algorithm can make its own
  # mind up about what the default should be, for example in simple it
  # is equal to the number of parameters.
  parser.add_argument('--population', action='store',
                      type=int, # no default, None signals default
    help="Set the population size of each generation")
  algorithm_choices = ["sa", "simple", "experiment"]
  parser.add_argument('--algorithm', action='store',
                      choices=algorithm_choices, default="simple",
    help="Select the genetic algorithm to deploy")
  parser.add_argument('--results-dir', action='store',
    help="Select a directory in which to store the results")
  parser.add_argument('--record-freq', action='store',
                      type=int, # no default, None signals no freq
    help="Specify how often results should be recorded")
               
 
  return parser


class OptimiserResults:
  """A class for the reporting of results of a run of an optimisation
     as distinct from reporting of the results of a flat experiment"""
  def __init__(self, best_citizen, best_cost):
    self.best_citizen = best_citizen
    self.best_cost = best_cost

  def print_best_citizen_results(self, results_dir):
    """Print the results for the best citizen if there are any"""
    if self.best_citizen.results:
      for results in self.best_citizen.results:
        basename = results.results_filename()
        best_timeseries_fname = os.path.join(results_dir, basename)
        best_timeseries_file = open(best_timeseries_fname, "w")
        results.write_to_file(best_timeseries_file)
        best_timeseries_file.close()
        print ("Best timeseries written to file: " + 
                best_timeseries_fname)
    else: 
      print ("No best candidate found")



  def report_results(self, configuration, results_dir):
    """Report on the results of a run of the optimisation"""
    optimisation = configuration.optimisation
    configuration.report_on_best_params(self.best_citizen.dictionary)
    
    print("After: " + str(configuration.num_generations) + " generations:")

    if not results_dir:
      params_dir = os.path.dirname(optimisation.params_file)
      results_dir = os.path.join(params_dir, "results")
    results_dir = utils.get_new_directory(results_dir)
  
    best_parameters_fname = os.path.join(results_dir, "best_params")
    best_parameters_file = open(best_parameters_fname, "w") 
    for param in optimisation.parameters:
      value = self.best_citizen.dictionary[param.name]
      line = param.name + ": " + str(value)
      print (line)
      best_parameters_file.write(line)
      best_parameters_file.write("\n")

    best_parameters_file.close()

    best_cost_fname = os.path.join(results_dir, "best_cost")
    best_cost_file = open(best_cost_fname, "w")
    best_cost_file.write(str(self.best_cost))
    best_cost_file.close()

    self.print_best_citizen_results(results_dir)

class ExperimentResults:
  """A class representing what we should do with the results of an
     a flat experiments (as opposed to an optimisation report)"""
  def __init__(self, results_file):
    self.results_file = results_file

  def report_results(self, configuration, results_dir):
    """Report on the results of the experiment to the user"""
    print ("Results data base written to: " + self.results_file)
    self.bundle_result(configuration, results_dir)

  def bundle_result(self, configuration, results_directory):
    """A simple function to bundle the results into the results
       directory."""
    if not results_directory:
      results_directory = utils.get_new_directory("results")
    shutil.copy2(self.results_file, results_directory)
    optimisation = configuration.optimisation
    shutil.copy2(optimisation.params_file, results_directory)
    # These are untested since I updated optimisation to hold
    # multiple model_datas rather than a single model file and
    # gold standard file.
    for model_data in optimisation.model_datas:
      shutil.copy2(model_data.model_file, results_directory)
      shutil.copy2(model_data.gold_standard_file, results_directory)
    print ("Find the results in: " + results_directory)
   

def run():
  """Process all the command line arguments and get going with the
     optimisation""" 
  # Parse in the command-line arguments
  parser    = create_arguments_parser(True)
  arguments = parser.parse_args()

  # I really just want to take in one file name which is a configuration
  # file which contains the others, but for now I'll do this:
  if len(arguments.filenames) < 3:
    print ("You must provide at least three files:")
    print ("   the initial parameter file")
    print ("   the model file")
    print ("   the gold standard data file")
    sys.exit(1)

  optimisation = get_optimisation_definition(arguments)
  configuration = get_configuration(arguments, optimisation)
  solve_model.initialise_logger(arguments)
  optimisation.initialise_solvers()

  algorithm = configuration.search_algorithm
  results = algorithm.run_optimisation(optimisation, configuration)

  configuration.report_on_solves()
  results.report_results(configuration, arguments.results_dir) 

if __name__ == "__main__":
  random.seed()
  run()

