"""A module to retain all of the code necessary for interfacing with
   the database. I've called it datastore so as not to use up the name
   database.
"""
import sqlite3
from contextlib import closing

def connect_db(databasefile=None):
  """A simple method to connect to the database"""
  if databasefile == None:
    databasefile = "/tmp/web-biopepa-latex.db"
  return sqlite3.connect(databasefile)

def init_db(app):
  """A simple method to initialise the database"""
  with closing(connect_db()) as database:
    with app.open_resource('schema.sql') as schema_file:
      database.cursor().executescript(schema_file.read())
    database.commit()


def get_model_field(database, model_id, field_name):
  """Return the given field value of a model's entry"""
  if database == None:
    database = connect_db()
  select_command = "select " + field_name + " from entries where ident=?"
  cur = database.execute(select_command, [model_id])
  field_value = cur.fetchone()
  cur.close()
  return field_value[0]


def get_model_source(database, model_id):
  """Obtain the model source from the database using the model's id"""
  return get_model_field(database, model_id, "modelsource")

def add_error(database, convert_id, message):
  """Add an error to the given model"""
  database.execute('update entries set errors=? where ident=?',
                   [message, convert_id])
  database.commit()

def update_model_field(database, model_id, field, value):
  """Update the given field of model entry with the given model identifier
     with the new value provided.
  """
  if database == None:
    database = connect_db()
  database.execute("update entries set " + field + "=? where ident=?",
                   [value, model_id])
  database.commit()


