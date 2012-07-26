"""A module/utility for numerically solving a model. Contains at least
   a solver using the sundials library converting SBML files into C
   programs, and also the biopepa command-line client
   for Bio-PEPA models"""
import argparse
import sys
import os.path
from subprocess import Popen, PIPE
import logging
import random
import math

import timeseries
import utils
import plotcsv
import parameters
import sbml_parametiser

import outline_sbml
import numpy as np
from scipy.integrate import odeint
import xml.dom.minidom
 
from solver_errors import SolverError, SolverModelError
import external_solvers

def get_time_grid(configuration):
  """From a solver configuration return the time points which should
     be returned from the solver
  """
  start_time = configuration.start_time
  stop_time = configuration.stop_time
  out_interval = configuration.out_interval
  # putting stop beyond the actual stop time means that the output
  # will actually include the stop time. Note that in some cases this
  # may result in stop_time + out_interval actually appearing in the
  # output as well, see numpy.arange documentation.
  return np.arange(start=start_time, stop=stop_time + out_interval,
                   step=out_interval)


class SBMLSolver(object):
  """Essentially an abstract class for SBML solvers which actually do
     the work internally. In particular we are trying to put any of the
     common work in here such that we don't repeat ourselves.
  """
  def __init__(self, model_file):
    self.model_file = model_file
    self._model = None

    self._params = None
    self._species = None
    self._reactions = None
    self._init_assigns = None

  # We use a bunch of property values here, this is so that separate
  # methods can refer to the list of reactions or species or the model,
  # but we don't parse them in from the model xml each time, and also we
  # do not need to pass them around as arguments. In theory we could just
  # initialise all of these and we wouldn't need to check but this is
  # a rather elegant way of ensuring that we never forget to initialise
  # one of these, or that we initialise them twice, or attempt to use
  # before initialising. This does mean that if you have code in a loop
  # which say reference the set of reactions you should alias it first:
  # reactions = self.reactions
  # for 1 in range (1000000):
  #   fire_all_the_reactions(reactions)

  @property
  def model(self):
    """Return the associated model, if the model file has already been
       parsed then we return the stored parsed model, if not then we first
       parse the model file and store the result before returning it.
    """
    if self._model:
      return self._model
    else:
      dom = xml.dom.minidom.parse(self.model_file)
      model = dom.getElementsByTagName("model")[0]
      self._model = model
      return model

  @property
  def params(self):
    """Return the parameters of the associated model, may even involve the
       parsing of the model file if this has not yet occurred.
    """
    if self._params:
      return self._params
    else:
      params = outline_sbml.get_list_of_parameters(self.model)
      self._params = params
      return params

  @property
  def species(self):
    """ Returns the species of the model again possibly parsing the model,
        again the species are then stored such that subsequent calls will
        not involve any re-parsing.
    """
    if self._species:
      return self._species
    else:
      species = outline_sbml.get_list_of_species(self.model)
      self._species = species
      return species

  @property
  def reactions(self):
    """Returns the reactions of the associated model, may parse the model
       file if it hasn't occured already
    """
    if self._reactions:
      return self._reactions
    else:
      reactions = outline_sbml.get_list_of_reactions(self.model)
      self._reactions = reactions
      return reactions

  @property
  def init_assigns(self):
    """ Returns the initial assignments of the associated model"""
    if self._init_assigns:
      return self._init_assigns
    else:
      init_assigns = outline_sbml.get_list_of_init_assigns(self.model)
      self._init_assigns = init_assigns
      return init_assigns


  def initialise_solver(self):
    """Initialise the scipy ode solver, nothing special required here"""
    pass

  def set_parameter_file(self, param_file):
    """This method is here to conform with the solver api, but cannot
       actually be used with anything but 'None' and hence will never
       have an effect.
    """
    # Ideally we should be able to parse in the param file and then
    # simply call parameterise_model
    if param_file:
      logging.error("Cannot use a parameter file with the scpiy solver")
      logging.info ("Model file : " + self.model_file)
      logging.error("Try --solver cvodes")
      raise ValueError

  def parameterise_model(self, param_file):
    """Given a parameter file modify the model (possibly in file)
       to obtain a new model with the parameters as specified in the
       param_file
    """
    dom = xml.dom.minidom.parse(self.model_file)
    model = dom.getElementsByTagName("model")[0]
    parameter_dict = parameters.parse_param_file(param_file)
    sbml_parametiser.parameterise_model(model, parameter_dict)

  def initialise_populations(self, population_dictionary,
                             rounded_species_names=None):
    """Initialises the population dictionary based on the SBML file.
       In particular we include, parameter definitions,
       species definitions and initial assignments.
       The third argument is the list of species names which you wish to
       have an integer population. So if you are implementing an ODE solver
       you probably don't wish to have any names here as a non-integer
       population is okay for an ODE solver. However if implementing a
       stochastic algorithm most populations should be integers. If you put
       the name in this list then here we choose an integer either above or
       below the given value with probability relative to how close, so for
       example the value 1.25 will 3/4 of the time be set to 1 and 1/4 of
       the time be set to 2.
    """
    # To avoid always having to check for None:
    if rounded_species_names == None:
      rounded_species_names = []
    # This superficially solves a problem with species taking a negative
    # value. What happened was if we had a value such as A = 120.3343
    # if A is decreased (always by 1 say) then eventually it gets to the
    # point where it the population of A is 0.3343, which means a reaction
    # in which A is involved is still enabled since the rate is not zero.
    # So we make all initial populations of species to be integer values.
    # However this does not completely solve the problem since the
    # stoichiometry may be greater than one, in which case if there is
    # exactly 1 left, then the reaction fires and now the population is -1.
    # I think the way to do this, without significantly affecting the speed
    # of the simulation is to first determine all the reactions in which
    # each species is a reactant. Then we take the least common multiple
    # of all the stoichiometries and then make sure that the initial
    # population is a multiple of these values. But actually that doesn't
    # work either. As long as 1 is a stoichiometry then at any time the
    # population may be any value, including 1. In particular if the
    # species is a product of a reaction with stoichiometry 1 and its
    # initial population is 0, then if we fire that reaction it will have
    # popultion 1, if it is then a reactant in a reaction with a
    # stoichiometry above 1 we're in trouble. So this needs some thought
    # and it may be that we can only solve this during the simulation. This
    # is slightly annoying as you may well have a value which can go below
    # zero.
    def round_initial_value(name, value):
      """If a value is a species then we assume that it should be given
         an integer value, this rounds the value appropriately if the
         name is in the list of species names.
      """
      if not name in rounded_species_names:
        return value
      else:
        floored = math.floor(value)
        remainder = value - floored
        # A value such as 2.3 should sometimes be rounded to 2 and
        # sometimes be rounded to 3. We take a random number and if it is
        # higher than the remainder part we round down and otherwise we
        # round up.
        if random.random() < remainder:
          return floored + 1
        else:
          return floored

    for param in self.params:
      # If the parameter value is not set here, we assume that it
      # is set in an initial assignment (if it isn't then we will later
      # fail if the parameter is used).
      if param.value:
        value = round_initial_value(param.name, float(param.value))
        population_dictionary[param.name] = value

    species = self.species
    for spec in species:
      if spec.initial_amount:
        population_dictionary[spec.name] = float(spec.initial_amount)

    init_assigns = self.init_assigns
    for init_assign in init_assigns:
      name   = init_assign.variable
      expr   = init_assign.expression
      value  = expr.get_value(environment=population_dictionary)
      value  = round_initial_value(name, value)
      population_dictionary[name] = value


  def optimise_expressions(self):
    """Optimises the expressions within the model, so in particular the
       reaction rate expressions are first parsed and then reduced
       (according to known constant values). This means that during the
       numerical evaluation they are far faster to evaluate than
       re-examining the xml. We also do this for initial assignments which
       might not seem like a big improvement but if we are performing many
       stochastic simulations it can be quite a speed-up.
    """
    # As a speed up, let's parse the reaction kinetic laws, we should
    # potentially do this in SBMLSolver as it may help the ODE solver
    # as well.
    # We now increase our optimisation by reducing the expressions based
    # upon a mapping for constant values. This allows us to reduce say:
    # R * (factor ^ 2)
    # to something like:
    # R * 4
    # if 'factor' is a constant with value '2'.
    # So first we must calculate the constant dictionary.
    constant_dictionary = dict()
    for param in self.params:
      # If the parameter value is not set here, we assume that it
      # is set in an initial assignment (if it isn't then we will later
      # fail if the parameter is used).
      if param.value and param.constant:
        constant_dictionary[param.name] = float(param.value)

    # Note that here we could also put the results into the constant
    # dictionary if they reduce to a single value and it is not a
    # species (which may change value at initialise time due to the
    # probabilistic rounding of non-integer populations, see above).
    for init_assign in self.init_assigns:
      # expr   = outline_sbml.parse_expression(init_assign.expression)
      # hmm, they already are expressions, we can still reduce them.
      expr   = init_assign.expression
      expr   = expr.reduce(constant_dictionary)
      init_assign.expression = expr

    for reaction in self.reactions:
      expr = outline_sbml.parse_expression(reaction.kinetic_law)
      expr = expr.reduce(constant_dictionary)
      reaction.kinetic_law = expr


   
class ScipyOdeSbmlSolver(SBMLSolver):
  """A class which implements an ODE based solver for SBML models.
     Based on scipy
  """
  # TODO: add use of the progress indicator here.
  def __init__(self, model_file):
    super(ScipyOdeSbmlSolver, self).__init__(model_file)

  def solve_model(self, configuration):
    """Solve the model, this does not call 'parameterise_model' so
       that should be called first if it is needed
    """
    population_dictionary = dict()
    self.initialise_populations(population_dictionary)
    # And also optimise the reaction rate expressions
    self.optimise_expressions()

    # As we have noted because 'self.reactions' is really a property
    # which does a check etc, it's worth aliasing it here to speed up
    # 'get_rhs'.
    reactions = self.reactions
    def get_rhs(current_pops, time):
      """The main function passed to the solver, it calculates from the
         current populations of species, the rate of change of each
         species. Also given a 'time' which may be used in the equations
         Essentially then solves for each ODE the right hand side of
         the ode at the given populations and time.
      """
      results = [0] * len(current_pops)
      population_dictionary["time"] = time
      for index in range(len(species_names)):
        population_dictionary[species_names[index]] = current_pops[index]
      for reaction in reactions:
        expr = reaction.kinetic_law
        rate = expr.get_value(environment=population_dictionary)
        # A negative rate is always possible, but it is usually because
        # a given species has been allowed to reduce to a negative
        # population. It is not quite clear what we should do in this
        # circumstance.
        # if rate < 0:
        #   print ("This rate is negative: " + str(rate))
        #   print ("  " + expr.show_expr())
        #   for name, value in population_dictionary.items():
        #     print (name + " = " + str(value))
          # print ("Exiting!")
          # sys.exit(1)
          # rate = 0.0

        reactants = reaction.reactants
        products = reaction.products
        for reactant in reactants:
          reactant_index = species_names.index(reactant.name) 
          results[reactant_index] -= (rate * reactant.stoich)
        for product in products:
          product_index = species_names.index(product.name)
          results[product_index] += (rate * product.stoich)
        
      return results  

    # We must set up the initial 'array' (here actually just a list,
    # perhaps we can speed this up using an actual array) we must initial
    # each species with the initial population, but we have already
    # initialised the population dictionary so all we need to do is for
    # each species look it up in the population dictionary and we need
    # not mind whether this initialisation is due to the species spec
    # or an initial assignment.
    species_names = [ s.name for s in self.species ]
    initials = [0] * len(species_names)
    for name, value in population_dictionary.items():
      # If the name is not in the species name then it is likely a
      # parameter or other value we need not concern outselves with here.
      if name in species_names:
        index = species_names.index(name)
        initials[index] = value
   
    # The time grid, I'm going to be honest I don't think I understand this.
    time_grid  = get_time_grid(configuration)
    # Solve the ODEs
    soln = odeint(get_rhs, initials, time_grid)

    timecourse = timeseries.Timeseries(species_names, soln)
    # For debugging purposes we'll just quickly plot the results
    # timecourse.plot_timecourse()

    return timecourse
   

def exponential_delay(mean):
  """From the given average length samples randomly to give an
     exponentially distributed delay.
     Remember, the argument here is the average *delay* so if you have
     a *rate* then take the recipriocal to get the average delay, eg:
     delay = exponential_delay(1.0 / rate)
  """
  return -mean * math.log(random.random())

class SimpleProgressIndicator(object):
  """This is a very simple progress indicator which I will use
     for keeping the user updated about the amount of a simulation done
  """
  def __init__(self, stop_amount):
    self.stop_amount = stop_amount
    self.total_stop_amount = stop_amount
    self.last_reported = 0
    self.precision = 100.0
    self._runs = 1
    self.completed_runs = 0

  @property
  def runs(self):
    """Return the total number of simulation runs which will be
       performed
    """
    return self._runs 

  @runs.setter
  def runs(self, runs):
    """Set the total number of simulation runs to be completed"""
    self._runs = runs
    self.total_stop_amount = self.stop_amount * runs

  def done_work(self, amount):
    """To be called upon completing some amount of work, the amount
       reported is the total amount done, so for example it would be
       the time in a simulation.
    """
    total_done = amount + (self.stop_amount * self.completed_runs)
    proportion_done = (total_done / self.total_stop_amount) * self.precision
    proportion_done = math.floor(proportion_done)
    percentage_done = (100.0 / self.precision) * proportion_done

    if percentage_done > self.last_reported:
      self.last_reported = percentage_done
      print (str(percentage_done) + "% done")


class SimulationRecorder(object):
  """A class to record simulation events, this is a simple recorder which
     only builds up a time series and only records those rows which
     we ultimately wish to output. But one could understand more
     interesting recorders which build up a file for say traviando.
  """
  def __init__(self, configuration):
    self.next_time_point = configuration.start_time
    self.interval   = configuration.out_interval
    self.last_time_point = configuration.stop_time
    self.species_names   = []
    self.timecourse_rows = []

  def record_event(self, time, population_dictionary):
    """Decide if we wish to record a new row (or rows) in the timeseries
       and do so if required.
    """
    # Using <= here for the first condition, since if the current time
    # exactly hits a time point then we probably wish to record it as is.
    # The second condition has less than because if we attempt to use
    # self.next_time_point <= self.last_time_point we suffer if the
    # interval doesn't exactly divide the stop time. Because then we do
    # not record the time at exactly the stop time, this still doesn't
    # do this quite correctly but we at least get a data point beyond
    # the stop time. We could ensure that we get exactly the final data
    # point by doing:
    #  this_row = [ min(self.next_time_point, self.last_time_point) ]
    while (self.next_time_point <= time and
           self.next_time_point < self.last_time_point + self.interval) :

      this_row = [ self.next_time_point ]
      for name in self.species_names:
        this_row.append(population_dictionary[name])
      self.timecourse_rows.append(this_row)
      self.next_time_point += self.interval

  def get_results(self):
    """Get the results so far, usually called after the end of the
       simulation
    """
    column_names = [ "Time" ] + self.species_names
    timecourse = timeseries.Timeseries(column_names, self.timecourse_rows)
    return timecourse


def dice_roll_choose_element(elements, values, sum_of_values):
  """Chooses from the list of elements based on their values and the sum
     of values. This is intended to choose a reaction from a list of
     reactions which have different rates. The larger the rate the more
     likely the reaction is to be chosen.
  """
  dice_roll = random.uniform(0, sum_of_values)
  accumulated_probability = 0.0
  index = 0
  for value in values:
    accumulated_probability += value
    if dice_roll < accumulated_probability:
      return elements[index]
    index += 1
  else:
    message = "Very bad, got to the end of a list of value to choose from"
    raise ValueError (message)

def log_simulation_state(pop_dictionary, reactions):
  """Logs, using logging.debug, the state given to it"""
  for reaction in reactions:
    expr = reaction.kinetic_law
    rate = expr.get_value(environment=pop_dictionary)
    logging.debug(reaction.format_reaction())
    logging.debug(str(rate))

  for name, value in pop_dictionary.items():
    logging.debug(name + " = " + str(value))

class StochasticSimulationSolver(SBMLSolver):
  """A class which implements a solver as stochastic simulation algorithm
  """
  def __init__(self, model_file):
    super(StochasticSimulationSolver, self).__init__(model_file)
    self.progress_indicator = None

  def solve_model(self, configuration):
    """Solve the model, this does not call 'parameterise_model' so
       that should be called first if it is needed
    """
    species_names = [ s.name for s in self.species ]
    # Again it's worth aliasing self.reactions because that is ultimately
    # a property which does a little checking before actually returning
    # the set of reactions which we do not wish to do every iteration of
    # the main loop. In fact we wouldn't do that any way, we would do it
    # each time we simulated the model though, which might be many many
    # times.
    reactions = self.reactions
    # Speed up the simulation by optimising the reaction rate expressions
    self.optimise_expressions()

    # Note that for the progress indicator we don't have to worry about
    # subtracting the start time from the stop time, since the start time
    # only affects what we record, not what we actually have to simulate.
    stop_time = configuration.stop_time
    self.progress_indicator = SimpleProgressIndicator(stop_time)
    self.progress_indicator.precision = configuration.progress_points
    self.progress_indicator.runs = configuration.runs

    timecourses = []
    for _ in range(configuration.runs):
      timecourse = self.simulate_model(configuration,
                                       species_names,
                                       reactions)
      timecourses.append(timecourse)
      self.progress_indicator.completed_runs += 1

    average_timecourse = timeseries.average_timeseries(timecourses)
    # average_timecourse.write_to_file(sys.stdout)
    return average_timecourse
    
  def simulate_model(self, configuration, species_names, reactions):
    """Simulate the model once, solve_model_calls this to perform one
       simulation run, it may call it several times and then average the
       results.
    """
    # Initialising the populations here, means that unfortunately we
    # do initialise the populations for each simulation run if there are
    # multiple runs. However recall that for stochastic simulation where
    # the initial population is something like: 69.3, we have to
    # probabilistically choose between 69 and 70. So we could do this
    # initialisation only once but we would have to decouple the logic
    # to probabilistically round the initial populations, this may actually
    # be worthwhile since then we don't need to worry about doing that for
    # the ODE simulator.
    population_dictionary = dict()
    self.initialise_populations(population_dictionary,
                                rounded_species_names=species_names)
    sim_recorder = SimulationRecorder(configuration)
    sim_recorder.species_names = species_names
    # We have to start the simulation at 0.0 and not
    # configuration.start_time, no matter where we start recording
    time = 0.0
    while time < configuration.stop_time:
      # add up all the rates of all the reactions
      rates = []
      for reaction in reactions:
        expr = reaction.kinetic_law
        rate = expr.get_value(environment=population_dictionary)
        # A negative rate is always possible, but it is usually because
        # a given species has been allowed to reduce to a negative
        # population. It is not quite clear what we should do in this
        # circumstance.
        if rate < 0:
          print ("This rate is negative: " + str(rate))
          print ("  " + expr.show_expr())
          for name, value in population_dictionary.items():
            print (name + " = " + str(value))
          rate = 0.0
        rates.append(rate)

      overall_rate = sum(rates)
      if overall_rate <= 0:
        print ("No reactions are live")
        time = configuration.stop_time
        break
      overall_delay = exponential_delay(1.0 / overall_rate)
      time += overall_delay
      self.progress_indicator.done_work(time)

      # Call the simulation recorder to, possibly record the current
      # populations, note that we do this *prior* to actually updating
      # the populations but *after* we have udpated the time for the
      # current delay. This is correct since the current populations will
      # be correct right up until the end of this delay so any time points
      # which occur between the current time before this delay and the time
      # plus this delay should be recorded as the current populations not
      # the updated ones. eg, if the current time was 0.23342 and the delay
      # was 0.3 making the 'time' now 0.53342, and lets say we're
      # recording every 0.1 time units, then then current populations
      # should be recorded for 0.3, 0.4 and 0.5 (0.2 should have already
      # been recorded previously)
      sim_recorder.record_event(time, population_dictionary)


      if time < 0:
        print ("--------")
        print (overall_rate)
        print (overall_delay)
        print (time)
        sys.exit(1)

      chosen_reaction = dice_roll_choose_element(reactions,
                                                 rates,
                                                 overall_rate)
      # Update population dictionary based on chosen reaction
      for reactant in chosen_reaction.reactants:
        population_dictionary[reactant.name] -= reactant.stoich
      for product in chosen_reaction.products:
        population_dictionary[product.name] += product.stoich
   
    # This is a somewhat annoying quirk of the logging facility,
    # I'd really like to query what the current level is and only
    # call this function based on that, this would avoid needlessly
    # scanning over the state and calling logging.debug when the loglevel
    # is set lower.
    log_simulation_state(population_dictionary, reactions)
      
    # End of the simulation, get the time course.
    timecourse = sim_recorder.get_results()
    # timecourse.write_to_file(sys.stdout)
    # For debugging purposes we'll just quickly plot the results
    # timecourse.plot_timecourse()

    return timecourse

def initialise_logger(arguments):
  """Initialise the logging system, depending on the arguments
     which may set the log level and a log file"""
  log_level = arguments.loglevel
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

  log_file = arguments.logfile
  if log_file:
    logging.basicConfig(filename=log_file, level=numeric_level)
  else:
    logging.basicConfig(level=numeric_level)


def create_arguments_parser(add_help):
  """Create the command-line arguments parser, takes in a single argument
     which is true if the created parser should add the --help flag. This
     should essentially only be true for a command that does not wish
     to add this parser as a parent parser (ie. if you do not wish to add
     further arguments)"""
  description = "Solve a single model one single time"
  epilog_usage_info = """
Some arguments are only valid for some solvers, for example the --runs
is only valid for the ssa solver, it doesn't really make any sense for
the deterministic solvers

The progress indicator is only implemented for some of the solvers,
essentially those that we have implemented natively here rather than
calling out to cvodes or the Bio-PEPA Eclipse compiler.
"""
  parser = argparse.ArgumentParser(add_help=add_help, 
                                   description=description,
                                   epilog = epilog_usage_info)

  solver_choices = ["cvodes", "biopepa", "scipy-ode", "ssa"]
  parser.add_argument('--solver', action='store',
                      choices=solver_choices,
                      help="Set the solver for numerical analysis")

  parser.add_argument('--start-time', action='store',
                      type=float, default=0.0,
                      help="Set the initial time of the numerical analysis")
  # Stop time cannot have a default value since the optimiser would
  # assume that the user has explictly set the stop time. We don't want
  # this because the optimiser should use the last data point as the
  # stop_time unless this is overridden by the user.
  parser.add_argument('--stop-time', action='store',
                      type=float, # default=1.0,
                      help="Set the stop time of the numerical analysis")
  parser.add_argument('--runs', action='store',
                      type=int, default=1,
    help="Number of indepenent runs, only used by the ssa solver")
  parser.add_argument('--progress-points', action='store',
                      type=int, default=1000,
    help="Number of progress indicator updates, 1000 gives 0.1% 0.2% etc")
  parser.add_argument('--reltol', action='store',
                      type=float, default=1.0e-6,
                      help="Set the solver's relative tolerance")
  parser.add_argument('--atol', action='store',
                      type=float, default=1.0e-6,
                      help="Set the solver's absolute tolerance")
  parser.add_argument('--interval', action='store',
                      type=float, default=0.0001,
                      help="Set the solver's internal time interval")
  parser.add_argument('--out-interval', action='store',
                      type=float, default=0.1,
                      help="Set the interval of the result's timecourse")
  parser.add_argument('--max-times', action='store',
                      type=int, default=10000000000,
                      help="Set the maximum number of computed times")
  parser.add_argument('--column', action=utils.ListArgumentAction,
                      help="Specify a column to be plotted")
  parser.add_argument('--mcolumn', action=utils.ListArgumentAction,
                      help="Specify a column not to be plotted")

  parser.add_argument('--plot-results', action='store_true',
                      help="Plot the resulting timeseries to a pdf file")
  log_choices = [ "info", "warning", "error", "critical", "debug" ]
  parser.add_argument('--loglevel', action='store',
                      choices=log_choices, default='info',
                      help="Set the level of the logger")
  parser.add_argument('--logfile', action='store',
                      help="The file to output the log to")
  return parser


def get_solver(filename, arguments):
  """Return the solver to be used based on the --solver flag given
     or, if that is not around, then the extension of the model file"""
  solver_name = arguments.solver
  if not solver_name:
    extension = os.path.splitext(filename)[1]
    if extension == ".xml" or extension == ".sbml":
      solver_name = "cvodes"
    elif extension == ".biopepa":
      solver_name = "biopepa"
    else:
      message = "Unknown filetype, should be: .xml, .sbml or .biopepa"
      logging.error (message)
      raise SolverError(message)
 
  if solver_name == "cvodes":
    # I've at least taken out the specificity to me, but now we are
    # a bit specific to Linux
    which_command = [ "which", "SBML2C" ]
    which_process = Popen(which_command, stdout=PIPE)
    which_output = which_process.communicate()[0]

    if which_process.returncode != 0:
      message = "which(SBML2C) process failed to return"
      logging.error (message)
      raise SolverModelError(message)

    bin_dir = os.path.dirname(which_output)
    install_dir = os.path.dirname(bin_dir)
    cflags_prefix = install_dir 
    solver = external_solvers.SbmlCvodeSolver(filename, cflags_prefix)
    return solver
  elif solver_name == "biopepa":
    biopepajar = "biopepa.jar"
    solver = external_solvers.BioPEPASolver(filename, biopepajar)
    return solver
  elif solver_name == "scipy-ode":
    solver = ScipyOdeSbmlSolver(filename)
    return solver
  elif solver_name == "ssa":
    solver = StochasticSimulationSolver(filename)
    return solver
  else:
    logging.error("Unknown solver name: " + solver_name)
    sys.exit(1)

def run():
  """Perform the banalities of command line processing then get on
     with the actual work"""
  parser    = create_arguments_parser(True)
   # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="an sbml file to solve numerically")
  parser.add_argument('--param-file', action='store',
    help="Provide a parameter file to override parameters in the sbml file")
  arguments = parser.parse_args()
  # stop_time cannot have a default because the optimiser needs to know
  # if the user has explicitly set the stop time or not, see the
  # add_argument call for 'stop_time' in create_arguments_parser.
  # So instead of a default value we check if it has been
  # set and if not we set it:
  if not arguments.stop_time:
    arguments.stop_time = 1.0

  initialise_logger(arguments)
  configuration = arguments

  for filename in arguments.filenames:
    try:
      if not os.path.exists(filename):
        logging.error("Model file: " + filename + " does not exist")
        continue
      solver = get_solver(filename, arguments)
      solver.initialise_solver()
      solver.set_parameter_file(arguments.param_file)
      timecourse = solver.solve_model(configuration)
      all_columns = timecourse.get_column_names()
      used_names = utils.get_non_ignored(all_columns, 
                                         arguments.column,
                                         arguments.mcolumn)
   
      if timecourse:
        # First remove any columns which are going to be plotted.
        for tc_column in all_columns:
          if tc_column not in used_names:
            timecourse.remove_column(tc_column)

        results_filename = utils.change_filename_ext(filename, ".csv")
        results_file = open(results_filename, "w")
        timecourse.write_to_file(results_file)
        results_file.close()

        if arguments.plot_results:
          plotcsv.run(argument_strings=[results_filename])  

      else:
        logging.error ("solving the model failed, no timeseries to report")
        continue
    except SolverError:
      # We should have already logged the error, now we just wish to
      # allow ourselves to continue and solve the remaining models
      continue

if __name__ == "__main__":
  run()
