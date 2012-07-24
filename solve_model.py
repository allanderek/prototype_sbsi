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


import biopepa
import timeseries
import utils
import plotcsv
import parameters
import sbml_parametiser

# for scipy ode solver which may be moved to a file on its own
import outline_sbml
import numpy as np
from scipy.integrate import odeint
import xml.dom.minidom

class SolverError(Exception):
  """A simple exception to be raised when we recognise that the model
     cannot be solved successfully. This allows callers, such as the
     optimisation routine, to catch this kind of error rather than
     fail completely. Note that this should be used to indicate an
     error in the actual solving of a particular instance of a model,
     not that the model itself cannot be solved. See SolverModelError
     for that.
  """
  def __init__(self, message):
    Exception.__init__(self)
    self.message = message

  def get_message(self):
    """Return the stored message"""
    return self.message

class SolverModelError(SolverError):
  """A simple subclass of SolverError defined above. Essentially,
     for the purposes of implementing many solves of a particular
     model we wish to distinguish between errors which mean that the
     model cannot be solved with the particular parameters, meaning
     for example if you change the parameters it might then be solvable
     and one in which it will never be solvable, for example in the
     Cvodes solver if we fail to be able to convert from SBML to C then
     we will never be able to solve that particular model.
  """
  pass


def run_command(command):
  """ A utility function to run a given command, log any ouput
      to standard out as debugging information, and any output to
      stderr as warning and return the returncode of the process.
  """
  process = Popen(command, stdout=PIPE, stderr=PIPE)
  output, errorout = process.communicate()

  if output:
    logging.debug(output)
  if errorout:
    logging.warning("The following command produced output on stderr")
    logging.warning(" ".join(command))
    logging.warning(errorout)

  return process.returncode

def run_command_with_error(command, exception):
  """Runs the given command as in 'run_command' however, raises
     the given exception if the return code indicates an error
     in running the given command.
  """
  return_code = run_command(command)
  if return_code != 0:
    logging.error(exception.get_message())
    logging.error(" ".join(command))
    raise exception

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
     the results. Additionally this also accepts as input a biopepa
     file, in which case we will first run the command to export the
     biopepa file to an SBML file.
  """
  def __init__(self, model_file, cflags_prefix):
    self.model_file = model_file
    self.model_exec = utils.change_filename_ext(model_file, ".exe")
    self.param_filename = None
    self.cflags_prefix = cflags_prefix

  def convert_biopepa(self):
    """Perform the conversion from biopepa to sbml if necessary """
    basename, extension = os.path.splitext(self.model_file)
    if extension == ".biopepa":
      new_model_file = basename + ".xml"
      biopepa_command = [ "java",
                          "-jar",
                          "biopepa.jar",
                          "export-sbml",
                          self.model_file,
                          "--output-file",
                          new_model_file,
                          "--no-warnings", 
                        ]
      error_message = "biopepa conversion to xml failed, command was:"
      error_exception = SolverModelError(error_message)
      run_command_with_error(biopepa_command, error_exception)
      self.model_file = new_model_file

  def run_sbml2c(self):
    """As part of the initialisation of the model we must run the
       sbml2c program over the model file
    """
    sbml2c_command = [ "SBML2C", self.model_file]
    sbml2c_returncode = run_command(sbml2c_command)
    if sbml2c_returncode != 0:
      logging.error ("SBML2C command has failed")
      exception = SolverModelError("SBML2C command has failed")
      raise exception

  # This should actually check if the model file is newer than
  # the model executable and if not then we needn't recompile it.
  def initialise_solver(self):
    """Initialise the sbml solver, for this we require the translation
       of the SBML model into a C program, and the compilation of that
       C program. Additionally if the model is a biopepa model file,
       then we must convert the biopepa into xml first.
    """
    self.convert_biopepa()
    self.run_sbml2c()
     

    # Obtain the directory in which the model file is, this is
    # necessary since we assume other files such as the main_RHS_Model.C
    # file are in the same directory.
    model_dir = os.path.dirname(self.model_file)

    base_include = os.path.join(self.cflags_prefix, "include") 
    include_flags = [ "-I" + base_include,
                      "-I" + os.path.join(base_include, "libxml2"),
                      "-I" + os.path.join(base_include, "sbsi_numeric"),
                    ]
    lib_flags = [ "-L" + os.path.join(self.cflags_prefix, "lib") ]
    extra_flags = [ "-DWL=32", "-DNO_UCF", ] 
    c_compiler = "mpic++"

    # Create the simple default main model C file
    main_rhs_model_cpath = os.path.join(model_dir, "main_RHS_Model.C")
    main_rhs_model_file = open(main_rhs_model_cpath, "w")
    main_rhs_model_file.write("#include \"UserModel/UserModel.h\"\n")
    main_rhs_model_file.write("#include <MainRHSTemplate.h>\n")
    main_rhs_model_file.close()


    # Run the first C Command
    first_c_command = [ c_compiler, "-o", 
                        os.path.join(model_dir, "main_RHS_Model.o"),
                        "-c",
                        main_rhs_model_cpath, 
                      ] + include_flags + lib_flags + extra_flags

    error_exception = SolverModelError("Failed to compile main_RHS_Model.C")
    run_command_with_error(first_c_command, error_exception)

 
    # Run the Second C command
    snd_c_command = [ c_compiler, "-o",
                      os.path.join(model_dir, "UserModel/UserModel.o"),
                      "-c",
                      os.path.join(model_dir, "UserModel/UserModel.C"),
                    ] + include_flags + lib_flags + extra_flags

    error_exception = SolverModelError("Failed to compile user model: ")
    run_command_with_error(snd_c_command, error_exception)

   
    # Run the final C command to link together the executable
    trd_c_command = [ c_compiler, "-o", self.model_exec,
                      os.path.join(model_dir, "main_RHS_Model.o"),
                      os.path.join(model_dir, "UserModel/UserModel.o"),
                      "-lsbsi_numeric", "-lsbml", "-lpgapack",
                      "-lfftw3", "-lsundials_kinsol",
                      "-lsundials_nvecserial", "-lsundials_cvode",
                      "-lxml2"
                    ] + lib_flags + extra_flags
    error_exception = SolverModelError("Failed to compile user model: ") 
    run_command_with_error(trd_c_command, error_exception)

  def set_parameter_file(self, param_filename):
    """Set the parameter file name, this allows us to use a parameter
       filename without going through the 'parameterise_model' method.
    """
    self.param_filename = param_filename

  def parameterise_model(self, dictionary):
    """Creates a param_overrides file from the given dictionary.
       The param_overrides file can be read by the main program
       which evaluates the model (see main_RHS_Model.C) to override
       the parameter values used in the model. This is the way in
       which we can parameterise the model, rather than changing the
       SBML model file and thus incurring the cost of a C compilation"""
    model_dir  = os.path.dirname(self.model_file)
    self.param_filename = os.path.join(model_dir,
                          os.path.join("UserModel", "param_overrides"))

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
    model_dir = os.path.dirname(self.model_exec)
    results_file_prefix = os.path.join(model_dir, "model")
    results_file = results_file_prefix + "_RHS.dat"
    if os.path.exists(results_file):
      # I'm not sure about this try-except block, I'm trying to
      # overcome what appears to be a bug in sbsi-numerics which creates
      # these weird files with nothing in them and setuid's with ????
      try:
        os.remove(results_file)
      except OSError:
        pass

    # This could also be run without mpi.
    mpirun_command = [ "mpirun",
                       self.model_exec,
                       "model_name", # Could get model name from xml file
                       str(configuration.stop_time),
                       str(configuration.start_time),
                       str(configuration.max_times),
                       str(configuration.interval),
                       str(configuration.out_interval),
                       str(configuration.atol),
                       str(configuration.reltol),
                       results_file_prefix,
                     ]
    if self.param_filename:
      mpirun_command.append(self.param_filename)
    mpi_returncode = run_command(mpirun_command)
                    
    # So we check if the results file actually exists and if
    # not we assume it failed. Also now I can actually check
    # the return code
    if not os.path.exists(results_file) or mpi_returncode != 0:
      message = "Model solving failed"
      logging.warning(message)
      raise SolverError(message)

    logging.debug("stop-time == " + str(configuration.stop_time))
    logging.debug("start-time == " + str(configuration.start_time))
    logging.debug("max-times == " + str(configuration.max_times))
    logging.debug("interval == " + str(configuration.interval))
    logging.debug("out-interval == " + str(configuration.out_interval))
    logging.debug("atol == " + str(configuration.atol))
    logging.debug("reltol == " + str(configuration.reltol))

    # We should remove the parameter overrides file here in case
    # we wish to simply solve the model on its own later.
    # This is actually solved because we have updated the solver file
    # such that it doesn't simply look to see if UserModel/param_overrides
    # file is there, but instead demands that it is specified on the
    # command-line if one is to be used.
    # if os.path.exists(self.param_filename):
    #   os.remove(self.param_filename)

    csv_file = open(results_file,  "r")
    timecourse = timeseries.parse_csv(csv_file, "\t")
    csv_file.close()
    return timecourse


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
    self.model = None

  def get_model(self):
    """Return the associated model, if the model file has already been
       parsed then we return the stored parsed model, if not then we first
       parse the model file and store the result before returning it.
    """
    if self.model:
      model = self.model
    else:
      dom = xml.dom.minidom.parse(self.model_file)
      model = dom.getElementsByTagName("model")[0]
      self.model = model

    return model



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

    model = self.get_model()
    for param in outline_sbml.get_list_of_parameters(model):
      # If the parameter value is not set here, we assume that it
      # is set in an initial assignment (if it isn't then we will later
      # fail if the parameter is used).
      if param.value:
        value = round_initial_value(param.name, float(param.value))
        population_dictionary[param.name] = value

    species = outline_sbml.get_list_of_species(model)
    for spec in species:
      if spec.initial_amount:
        population_dictionary[spec.name] = float(spec.initial_amount)

    init_assigns = outline_sbml.get_list_of_init_assigns(model)
    for init_assign in init_assigns:
      name   = init_assign.variable
      expr   = init_assign.expression
      value  = expr.get_value(environment=population_dictionary)
      value  = round_initial_value(name, value)
      population_dictionary[name] = value


  def optimise_reaction_rates(self, reactions):
    model = self.get_model()
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
    for param in outline_sbml.get_list_of_parameters(model):
      # If the parameter value is not set here, we assume that it
      # is set in an initial assignment (if it isn't then we will later
      # fail if the parameter is used).
      if param.value and param.constant:
        constant_dictionary[param.name] = float(param.value)

    for reaction in reactions:
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
    model = self.get_model()

    reactions = outline_sbml.get_list_of_reactions(model)
    
    population_dictionary = dict()
    self.initialise_populations(population_dictionary)
    # And also optimise the reaction rate expressions
    self.optimise_reaction_rates(reactions)


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
    species_names = [ s.name 
                      for s in outline_sbml.get_list_of_species(model)]
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


def get_time_points(configuration):
  """From the configuration, start_time, stop_time and out_interval
     figure out the list of times which should be reported.
  """
  new_times = []
  current_time = configuration.start_time
  stop_time = configuration.stop_time
  out_interval = configuration.out_interval
  while current_time < stop_time:
    new_times.append(current_time)
    current_time += out_interval 

  # This ensures that the actual stop time is included.
  new_times.append(stop_time)

  return new_times

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
    return self._runs 

  @runs.setter
  def runs(self, runs):
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
  def __init__(self, time_points):
    self.time_points     = time_points
    self.next_time_index = 0
    self.next_time       = time_points[0]
    self.species_names   = []
    self.timecourse_rows = []

  def record_event(self, time, population_dictionary):
    """Decide if we wish to record a new row (or rows) in the timeseries
       and do so if required.
    """
    while self.next_time <= time:
      this_row = [ self.next_time ]
      for name in self.species_names:
        this_row.append(population_dictionary[name])
      self.timecourse_rows.append(this_row)
      self.next_time_index += 1
      # This is a little fragile if the stop_time is not in the
      # time_points then unfortunately we will 
      if self.next_time_index >= len(self.time_points):
        break
      else:
        self.next_time = self.time_points[self.next_time_index]

  def get_results(self):
    """Get the results so far, usually called after the end of the
       simulation
    """
    column_names = [ "Time" ] + self.species_names
    timecourse = timeseries.Timeseries(column_names, self.timecourse_rows)
    return timecourse



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
    model = self.get_model()
    species_names = [ s.name 
                        for s in outline_sbml.get_list_of_species(model)
                    ]
    reactions = outline_sbml.get_list_of_reactions(model)
    # Speed up the simulation by optimising the reaction rate expressions
    self.optimise_reaction_rates(reactions)

    # Note that for the progress indicator we don't have to worry about
    # subtracting the start time from the stop time, since the start time
    # only affects what we record, not what we actually have to simulate.
    stop_time = configuration.stop_time
    self.progress_indicator = SimpleProgressIndicator(stop_time)
    self.progress_indicator.precision = configuration.progress_points
    self.progress_indicator.runs = configuration.runs

    timecourses = []
    for i in range(configuration.runs):
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

    # now the actual ssa algorithm
    time = configuration.start_time
    new_times = get_time_points(configuration)
    sim_recorder = SimulationRecorder(new_times)
    sim_recorder.species_names = species_names

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
          # print ("Exiting!")
          # sys.exit(1)
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

      if time < 0:
        print ("--------")
        print (overall_rate)
        print (overall_delay)
        print (time)
        sys.exit(1)

      dice_roll = random.uniform(0, overall_rate)
      accumulated_probability = 0.0
      for index in range(len(rates)):
        accumulated_probability += rates[index]
        if dice_roll < accumulated_probability:
          chosen_reaction = reactions[index]
          break
      else:
        message = "Very bad, got to the end of all the reactions????"
        raise Exception(message)

      # Update population dictionary based on chosen reaction
      for reactant in chosen_reaction.reactants:
        population_dictionary[reactant.name] -= reactant.stoich
      for product in chosen_reaction.products:
        population_dictionary[product.name] += product.stoich
      
      # Call the simulation record to record the time row, though it
      # may only call selected time rows, eg the ones we wish to output.
      sim_recorder.record_event(time, population_dictionary)

    # End of the simulation, get the time course.
    timecourse = sim_recorder.get_results()
    # timecourse.write_to_file(sys.stdout)
    # For debugging purposes we'll just quickly plot the results
    # timecourse.plot_timecourse()

    return timecourse


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

  def set_parameter_file(self, param_file):
    """This method is here to conform with the solver api, but cannot
       actually be used with anything but 'None' and hence will never
       have an effect.
    """
    # Ideally we should be able to parse in the param file and then
    # simply call parameterise_model
    if param_file:
      logging.error("Cannot use a parameter file with the Bio-PEPA solver")
      logging.info ("Model file : " + self.model_file)
      logging.error("Try --solver cvodes")
      raise ValueError

  def parameterise_model(self, dictionary):
    """A simple and obviously broken method for parameterising a
       Bio-PEPA file. We assume that each parameter definition begins
       a line, with 'param_name = ....' and we simply replace the dots
       with the new value of the parameter."""
    biopepa.dumb_parameteriser.parameterise_model_file(
              dictionary,
              self.model_file,
              self.paramed_file)

  def solve_model(self, configuration):
    """Solve the parameterised version of the model. This assumes
       that parameterise_model has already been called and that hence
       a parameterised version of the model exists in a file named
       self.paramed_file"""

    start_time = configuration.start_time
    stop_time = configuration.stop_time
    data_points = int((stop_time - start_time) / 
                       configuration.out_interval)

    biopepa_command = [ "java",
                        "-jar",
                        "biopepa.jar",
                        "timeseries",
                        self.paramed_file,
                        "--no-warnings", 
                        "--solver",
                        "dopr-adaptive",
                        "--timeStep",
                        str(configuration.interval),
                        "--startTime",
                        str(start_time),
                        "--stopTime",
                        str(stop_time),
                        "--dataPoints",
                        str(data_points),
                        # "--output-file",
                        # csvfile 
                      ]


    biopepa_returncode = run_command(biopepa_command)
    biopepa_process = Popen(biopepa_command, stdout=PIPE, stderr=PIPE)
    output, errorout = biopepa_process.communicate()
 
    if errorout:
      logging.error ("The biopepa process produced output on stderr")
      logging.error (errorout)
      logging.error ("The biopepa command was: ")
      logging.error (" ".join(biopepa_command))

    if biopepa_returncode != 0:
      logging.error ("biopepa process failed to return")
      sys.exit(1)
    # A bug in pylint causes it to complain about this, it incorrectly
    # thinks that 'output' is of type list, but it is of type string.
    output_lines = output.split("\n")
    timecourse = timeseries.parse_csv(output_lines.__iter__(), ", ")

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
    solver = SbmlCvodeSolver(filename, cflags_prefix)
    return solver
  elif solver_name == "biopepa":
    biopepajar = "biopepa.jar"
    solver = BioPEPASolver(filename, biopepajar)
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
