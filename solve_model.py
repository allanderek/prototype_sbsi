"""A module/utility for numerically solving a model. Contains at least
   a solver using the sundials library converting SBML files into C
   programs, and also the biopepa command-line client
   for Bio-PEPA models"""
import argparse
import sys
import os.path
from subprocess import Popen, PIPE
import logging


import biopepa
import timeseries
import utils
import plotcsv
import parameters
import sbml_parametiser

# for scipy ode solver which may be moved to a file on its own
import outline_sbml
import numpy as np
import matplotlib.pyplot as plt
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
  # putting stop beyone the actual stop time means that the output
  # will actually include the stop time. Note that in some cases this
  # may result in stop_time + out_interval actually appearing in the
  # output as well, see numpy.arange documentation.
  return np.arange(start=start_time, stop=stop_time + out_interval,
                   step=out_interval)
 
class ScipyOdeSbmlSolver(object):
  """A class which implements an ODE based solver for SBML models.
     Based on scipy
  """
  def __init__(self, model_file):
    self.model_file = model_file
    self.model = None

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


  def solve_model(self, configuration):
    """Solve the model, this does not call 'parameterise_model' so
       that should be called first if it is needed
    """
    if self.model:
      model = self.model
    else:
      dom = xml.dom.minidom.parse(self.model_file)
      model = dom.getElementsByTagName("model")[0]

    species = outline_sbml.get_list_of_species(model)
    reactions = outline_sbml.get_list_of_reactions(model)
    species_names = [ s.name for s in species ]

    population_dictionary = dict()
    for param in outline_sbml.get_list_of_parameters(model):
      population_dictionary[param.name] = float(param.value)
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
        reactants = reaction.reactants
        products = reaction.products
        expr_evaluator = outline_sbml.ExprEvaluator()
        expr_evaluator.name_mapping = population_dictionary
        expr_evaluator.visit_maths(reaction.kinetic_law) 
        rate = expr_evaluator.get_results()
        reactant_indices = [ species_names.index(reactant.name) 
                               for reactant in reactants ]
        # for reactant_index in reactant_indices:
        #   rate *= current_pops[reactant_index]
        for reactant_index in reactant_indices:
          results[reactant_index] -= rate
        for product in products:
          product_index = species_names.index(product.name)
          results[product_index] += rate
        
      return results  

    # The initial conditions
    initials = [0] * len(species)
    for spec in species:
      if spec.initial_amount:
        index = species_names.index(spec.name)
        initials[index] = float(spec.initial_amount)


    init_assigns = outline_sbml.get_list_of_init_assigns(model)
    for init_assign in init_assigns:
      name = init_assign.variable
      index = species_names.index(name)
      initials[index] = init_assign.expression.get_value()
    
    # The time grid, I'm going to be honest I don't think I understand this.
    time_grid  = get_time_grid(configuration)
    # Solve the ODEs
    soln = odeint(get_rhs, initials, time_grid)

    # For debugging purposes we'll just quickly plot the results
    plt.figure()
    for index in range(len(species_names)):
      name = species_names[index]
      timecourse = soln[:, index]
      plt.plot(time_grid, timecourse, label=name)

    plt.xlabel('Time')
    plt.ylabel('Population')
    plt.title('Straight-forward MM')
    plt.legend(loc=0)
    plt.show()




    timecourse = timeseries.Timeseries(species_names, soln)
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
  parser = argparse.ArgumentParser(add_help=add_help, 
                                   description=description)

  parser.add_argument('--solver', action='store',
                      choices=["cvodes", "biopepa", "scipy-ode" ],
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
