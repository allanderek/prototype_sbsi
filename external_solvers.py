""" A helper module for solve_model which implements external solvers,
    that is solvers in which we do not actually perform any of the work
    ourselves but rather call out to third-party programs. Note this is
    not simply calling a library function but literally a separate process
"""
import sys
import os.path
from subprocess import Popen, PIPE
import logging

import biopepa
import timeseries
import utils

from solver_errors import SolverError, SolverModelError

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


