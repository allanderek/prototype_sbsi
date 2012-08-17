""" Executes the given command with the given timeout """
import sys
import sqlite3
import parcon

import biopepa.biopepa_to_latex as biopepa_to_latex
import biopepa.biopepa_to_sbml as biopepa_to_sbml

import timeouts

def connect_db():
  """A simple method to connect to the database"""
  return 

def get_model_source(database, convert_id):
  """Obtain the model source from the database using the model's id"""
  cur = database.execute('select modelsource from entries where ident=?',
                         [convert_id])
  model_source = cur.fetchone()
  cur.close()
  return model_source

def add_error(database, convert_id, message):
  """Add an error to the given model"""
  database.execute('update entries set errors=? where ident=?',
                   [message, convert_id])
  database.commit()

def update_model_field(database, model_id, field, value):
  """Update the given field of model entry with the given model identifier
     with the new value provided.
  """
  database.execute("update entries set " + field + "=? where ident=?",
                   [value, model_id])
  database.commit()

def convert_to_latex(timeout, database, arguments):
  """ The method for converting a model to LaTeX """
  convert_id = int(arguments[0])
  model_source = get_model_source(database, convert_id)

  # Basically makes the method the same but with a timeout of
  # the given timeout.
  convert_method = timeouts.timeout(timeout)
  convert = convert_method(biopepa_to_latex.convert_source)

  try:
    new_text = convert(model_source[0])
    update_model_field(database, convert_id, "latex", new_text)
  except parcon.ParseException as parse_exception:
    message = "Operation to convert to LaTeX failed due to a parse error: "
    message += parse_exception.message
    add_error(database, convert_id, message)
  except timeouts.TimeoutError:
    message = "Operation to convert to LaTeX timed-out"
    add_error(database, convert_id, message)


def convert_to_sbml(timeout, database, arguments):
  """ The method for converting a model to SBML """
  convert_id = int(arguments[0])
  model_source = get_model_source(database, convert_id)

  # Basically makes the method the same but with a timeout of
  # the given timeout.
  convert_method = timeouts.timeout(timeout)
  convert = convert_method(biopepa_to_sbml.convert_source)

  try:
    new_text = convert(model_source[0])
    update_model_field(database, convert_id, "modelsbml", new_text)
  except parcon.ParseException as parse_exception:
    message = "Operation to convert to SBML failed due to a parse error: "
    message += parse_exception.message
    add_error(database, convert_id, message)
  except timeouts.TimeoutError:
    message = "Operation to convert to SBML timed-out"
    add_error(database, convert_id, message)

def run():
  """ The simple main method"""
  # script_name = sys.argv[0]
  command_name = sys.argv[1]
  timeout = int(sys.argv[2])
  arguments = sys.argv[3:]

  # This kind of assumes that all of these commands will execute
  # something on the database, which is probably true, for those that
  # don't we can probably have a completely separate module for this.
  database = sqlite3.connect('/home/aclark6/tmp/web-biopepa-latex.db')

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

