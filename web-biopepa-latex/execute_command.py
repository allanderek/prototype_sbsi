""" Executes the given command with the given timeout """
import sys
import parcon

import biopepa.biopepa_to_latex as biopepa_to_latex
import biopepa.biopepa_to_sbml as biopepa_to_sbml

import timeouts
import datastore

def convert_to_latex(timeout, database, arguments):
  """ The method for converting a model to LaTeX """
  convert_id = int(arguments[0])
  model_source = datastore.get_model_source(database, convert_id)

  # Basically makes the method the same but with a timeout of
  # the given timeout.
  convert_method = timeouts.timeout(timeout)
  convert = convert_method(biopepa_to_latex.convert_source)

  try:
    new_text = convert(model_source)
    datastore.update_model_field(database, convert_id, "latex", new_text)
  except parcon.ParseException as parse_exception:
    message = "Operation to convert to LaTeX failed due to a parse error: "
    message += str(parse_exception)
    datastore.update_model_field(database, convert_id, "latex", "")
    datastore.add_error(database, convert_id, message)
  except timeouts.TimeoutError:
    message = "Operation to convert to LaTeX timed-out"
    datastore.update_model_field(database, convert_id, "latex", "")
    datastore.add_error(database, convert_id, message)


def convert_to_sbml(timeout, database, arguments):
  """ The method for converting a model to SBML """
  convert_id = int(arguments[0])
  model_source = datastore.get_model_source(database, convert_id)

  # Basically makes the method the same but with a timeout of
  # the given timeout.
  convert_method = timeouts.timeout(timeout)
  convert = convert_method(biopepa_to_sbml.convert_source)

  try:
    new_text = convert(model_source)
    datastore.update_model_field(database, convert_id,
                                 "modelsbml", new_text)
  except parcon.ParseException as parse_exception:
    message = "Operation to convert to SBML failed due to a parse error: "
    message += str(parse_exception)
    datastore.update_model_field(database, convert_id, "modelsbml", "")
    datastore.add_error(database, convert_id, message)
  except timeouts.TimeoutError:
    message = "Operation to convert to SBML timed-out"
    datastore.update_model_field(database, convert_id, "modelsbml", "")
    datastore.add_error(database, convert_id, message)

def run():
  """ The simple main method"""
  # script_name = sys.argv[0]
  command_name = sys.argv[1]
  timeout = int(sys.argv[2])
  arguments = sys.argv[3:]

  # This kind of assumes that all of these commands will execute
  # something on the database, which is probably true, for those that
  # don't we can probably have a completely separate module for this.
  database = datastore.connect_db()

  # This could be done like this but I prefer explicitly writing it out.
  # globals()[command_name](timeout, database, arguments)
  if command_name == "convert_to_sbml":
    convert_to_sbml(timeout, database, arguments)
  elif command_name == "convert_to_latex":
    convert_to_latex(timeout, database, arguments)
  else:
    print ("Command unknown: " + command_name)
  database.close()
 

if __name__ == "__main__":
  run ()

