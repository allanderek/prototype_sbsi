"""A module/utility for numerically solving a model. Contains at least
   a solver using the sundials library converting SBML files into C
   programs, and also the biopepa command-line client
   for Bio-PEPA models"""
import argparse
import sys
import os.path
from subprocess import Popen, PIPE
import logging
import timeseries
import utils

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
    self.param_filename = os.path.join("UserModel", "param_overrides")
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
      biopepa_process = Popen(biopepa_command)
      biopepa_process.communicate()
      if biopepa_process.returncode != 0:
        logging.error("biopepa conversion to xml failed, command was:")
        logging.error(biopepa_command)
      self.model_file = new_model_file
 
  # This should actually check if the model file is newer than
  # the model executable and if not then we needn't recompile it.
  def initialise_solver(self):
    """Initialise the sbml solver, for this we require the translation
       of the SBML model into a C program, and the compilation of that
       C program. Additionally if the model is a biopepa model file,
       then we must convert the biopepa into xml first.
    """
    self.convert_biopepa()
     
    sbml2c_command = [ "SBML2C", self.model_file]
    sbml2c_process = Popen(sbml2c_command)
    sbml2c_process.communicate()

    if sbml2c_process.returncode != 0:
      logging.error ("SBML2C command has failed")
      sys.exit(1)
 
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

    first_c_command = [ c_compiler, "-o", 
                        os.path.join(model_dir, "main_RHS_Model.o"),
                        "-c",
                        os.path.join(model_dir, "main_RHS_Model.C"),
                      ] + include_flags + lib_flags + extra_flags
    first_c_process = Popen(first_c_command)
    first_c_process.communicate()

    if first_c_process.returncode != 0:
      logging.error ("Failed to compile main_RHS_Model.C: ")
      logging.error (" ".join(first_c_command))
      sys.exit(1)

    snd_c_command = [ c_compiler, "-o",
                      os.path.join(model_dir, "UserModel/UserModel.o"),
                      "-c",
                      os.path.join(model_dir, "UserModel/UserModel.C"),
                    ] + include_flags + lib_flags + extra_flags
    snd_c_process = Popen(snd_c_command)
    snd_c_process.communicate()
    if snd_c_process.returncode != 0:
      logging.error ("Failed to compile user model: ")
      logging.error (" ".join(snd_c_command))
      sys.exit(1)

    self.model_exec = os.path.join(model_dir, self.model_exec)
    trd_c_command = [ c_compiler, "-o", self.model_exec,
                      os.path.join(model_dir, "main_RHS_Model.o"),
                      os.path.join(model_dir, "UserModel/UserModel.o"),
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
    model_dir = os.path.dirname(self.model_exec)
    results_file_prefix = os.path.join(model_dir, "model")
    results_file = results_file_prefix + "_RHS.dat"
    if os.path.exists(results_file):
      os.remove(results_file)

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
    mpi_process = Popen(mpirun_command, stdout=PIPE)
    mpi_output = mpi_process.communicate()[0]
    logging.debug(mpi_output)
                    
    # So we check if the results file actually exists and if
    # not we assume it failed. Also now I can actually check
    # the return code
    if not os.path.exists(results_file) or mpi_process.returncode != 0:
      logging.warning("Model solving failed")
      return None 

    # We should remove the parameter overrides file here in case
    # we wish to simply solve the model on its own later.
    if os.path.exists(self.param_filename):
      os.remove(self.param_filename)

    csv_file = open(results_file,  "r")
    timecourse = timeseries.parse_csv(csv_file, "\t")
    csv_file.close()
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
    # print (biopepa_command)
    biopepa_process = Popen(biopepa_command, stdout=PIPE)
    output = biopepa_process.communicate()[0]

    if biopepa_process.returncode != 0:
      logging.error ("biopepa process failed to return")
      sys.exit(1)
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
                      choices=["cvodes", "biopepa"],
                      help="Set the solver for numerical analysis")

  parser.add_argument('--start_time', action='store',
                      type=float, default=0.0,
                      help="Set the initial time of the numerical analysis")
  # Stop time cannot have a default value since the optimiser would
  # assume that the user has explictly set the stop time. We don't want
  # this because the optimiser should use the last data point as the
  # stop_time unless this is overridden by the user.
  parser.add_argument('--stop_time', action='store',
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
  parser.add_argument('--out_interval', action='store',
                      type=float, default=0.1,
                      help="Set the interval of the result's timecourse")
  parser.add_argument('--max_times', action='store',
                      type=int, default=10000000000,
                      help="Set the maximum number of computed times")
  parser.add_argument('--column', action='append',
                      help="Specify a column to be plotted")
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
      logging.error ("Unknown filetype, should be: .xml, .sbml or .biopepa")
      sys.exit(1)
 
  if solver_name == "cvodes":
    # Clearly this bit should not be so specific to me.
    home_dir = "/afs/inf.ed.ac.uk/user/a/aclark6/" 
    cflags_prefix = os.path.join (home_dir,
                                  "Source/svn-git-sbsi/install/")
    solver = SbmlCvodeSolver(filename, cflags_prefix)
    return solver
  elif solver_name == "biopepa":
    biopepajar = "biopepa.jar"
    solver = BioPEPASolver(filename, biopepajar)
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
    solver = get_solver(filename, arguments)
    solver.initialise_solver()
    timecourse = solver.solve_model(configuration)
    if timecourse:
      if arguments.column:
        tc_columns = timecourse.get_column_names()
        for tc_column in tc_columns:
          if tc_column not in arguments.column:
            timecourse.remove_column(tc_column)

      results_filename = utils.change_filename_ext(filename, ".csv")
      results_file = open(results_filename, "w")
      timecourse.write_to_file(results_file)
      results_file.close()
    else:
      print ("Error: solving the model failed, no timeseries to report")
      sys.exit(1) 

if __name__ == "__main__":
  run()
