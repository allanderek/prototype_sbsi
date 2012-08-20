""" A web application for translating Bio-PEPA files into LaTeX"""
import argparse
from subprocess import Popen
from flask import Flask, request, g, redirect, url_for, \
     abort, render_template, flash
import flask.ext.login as flasklogin

import biopepa.biopepa_to_latex as biopepa_to_latex
import datastore

# create our little application :)
# app = Flask("web-biopepa-latex")
app = Flask(__name__)

class Configuration(object):
  """A simple object to hold the configuration"""
  def __init__(self):
    # pylint complains about all capitals but that is exactly what the
    # configuration from_object does so we disable the warning here.
    # pylint: disable=C0103
    self.DEBUG = True
    self.SECRET_KEY = "development key of biopepa latex"
configuration_object = Configuration()
app.config.from_object(configuration_object)


# Create the login code for our little application
login_manager = flasklogin.LoginManager()
login_manager.init_app(app)
login_manager.login_view = "login"


class User(flasklogin.UserMixin):
  """The class defining a user object."""
  def __init__(self, ident):
    self.ident = ident
    self.name = None

  def get_id(self):
    """Return the id of this user"""
    return self.ident

  def __repr__(self):
    return "%d/%s" % (self.ident, self.name)

# callback to reload the user object        
@login_manager.user_loader
def load_user(userid):
  """The call back the login manager uses to actually load the user object.
     This essentially must take the unicode userid and return a User object
     which represents that user.
  """
  g.db = datastore.connect_db()
  cur = g.db.execute("select name from users where ident=?", [userid])
  users = cur.fetchall()
  if not users:
    return flasklogin.AnonymousUser()

  user_tuple = users[0]
  user = User(userid)
  user.name = user_tuple[0]
  
  return user


@app.before_request
def before_request():
  """Occurs before every request, we make a connection to the database"""
  g.db = datastore.connect_db()

@app.teardown_request
def teardown_request(_exception):
  """At the end of dealing with every request we close the connection to
     the database."""
  g.db.close()


@app.route('/')
@flasklogin.login_required
def show_entries():
  """ Handles requests for the show entries, that is shows all the public
      models to the user.
  """
  # First get the current user:
  current_user = flasklogin.current_user 
  owner_id = current_user.get_id()
  fields = "ident, visibility, title, errors, modelsource, latex, modelsbml"
  from_part = "from entries"
  where_part = "where owner_id=?"
  order_part = "order by ident desc"
  command = " ".join(["select", fields, from_part, where_part, order_part])
  cur = g.db.execute(command, [owner_id])
  entries = [ dict(ident=row[0], 
                   visibility=row[1],
                   title=row[2],
                   errors=row[3],
                   modelsource=row[4],
                   latex=row[5],
                   modelsbml=row[6])
              for row in cur.fetchall()]
  return render_template('show_entries.html',
                         entries=entries,
                         user=current_user)


def create_insert_command(table_name, names):
  """Creates an insert command for the data base given a set of names
     for the element to create and their associated values"""
  # So lets say names is something like ["title", "text", "owner"]
  # data_form will be "(title, text, owner)"
  data_form  = "(" + ", ".join(names) + ")"
  # and question_marks would simply be a list of 3 '?'s
  question_marks = ["?"] * len(names)
  # and data_holes would be "(?, ?, ?)"
  data_holes = "(" + ", ".join(question_marks) + ")"
  # So we should return in our example case:
  # "insert into table (title, text, owner) values (?,?,?)
  return " ".join(["insert into", table_name,
                   data_form,
                   "values", data_holes])

@app.route('/add', methods=['POST'])
@flasklogin.login_required
def add_entry():
  """Handles request to add a new entry"""
  current_user = flasklogin.current_user 
  
  owner_id = current_user.get_id()
  visibility = "private"
  title = request.form['title']
  modelsource = request.form['modelsource']
  
  insert_command = create_insert_command("entries", [ "owner_id",
                                                      "visibility",
                                                      "title",
                                                      "modelsource" ]) 
  values = [ owner_id, visibility, title, modelsource ]
  g.db.execute(insert_command, values)
  g.db.commit()
  flash('New entry was successfully posted')
  return redirect(url_for('show_entries'))

@app.route('/convert-to-latex', methods=['POST'])
def convert_to_latex():
  """Handles conversion requests, that is a request to convert a model
     from Bio-PEPA to LaTeX
  """
  convert_id = request.form["convert_id"]
  new_text = "Calculating ..."
  datastore.update_model_field(None, convert_id, "latex", new_text)
  command = [ "python", "web-biopepa-latex/execute_command.py", 
              "convert_to_latex", "5", convert_id ]
  _process = Popen(command)
  return redirect(url_for('show_entries'))
 
@app.route('/convert-to-sbml', methods=['POST'])
def convert_to_sbml():
  """Handles conversion to sbml requests, that is a request to convert a
     model from Bio-PEPA to SBML
  """
  convert_id = request.form["convert_id"]
  new_text = "Calculating ..."
  datastore.update_model_field(None, convert_id, "modelsbml", new_text)
  command = [ "python", "web-biopepa-latex/execute_command.py", 
              "convert_to_sbml", "5", convert_id ]
  _process = Popen(command)
  return redirect(url_for('show_entries'))
 
  


@app.route('/delete', methods=['POST'])
def delete_entry():
  """Handles requests to delete a given entry/model"""
  convert_id = request.form["convert_id"]
  g.db.execute('delete from entries where ident=?', [convert_id])
  g.db.commit()
  # flash ("You just deleted the identity: " + convert_id)
  return redirect(url_for('show_entries'))

@app.route('/no-store-biopepa-latex', methods=['GET', 'POST'])
def no_store_biopepa_latex():
  """Handles requests, both POST and GET to the no store biopepa latex,
     allowing users to convert their models without necessarily uploading
     their model to be stored on our server
  """
  error = None
  model_source = None
  latex = None
  if request.method == 'POST':
    model_source = request.form['modelsource']
    latex = biopepa_to_latex.convert_source(model_source)
  current_user = flasklogin.current_user 
  return render_template('raw-to-latex.html', 
                         source=model_source,
                         latex=latex, error=error,
                         user=current_user)

def goto_anonymous_login(error=None):
  """A utility method to send the user to the login page as an anonymous
     user. Essentially call this whenever the user requires to be logged
     in for some action but isn't currently.
  """
  return render_template('login.html',
                         error=error,
                         user=flasklogin.AnonymousUser())


@app.route('/login', methods=['GET', 'POST'])
def login():
  """The login handler, handles both the GET where it display the login
     in form, and the POST where it logs in the user assuming their
     credentials pass muster.
  """
  if request.method == 'POST':
    username = request.form['username']
    password = request.form['password']
    cur = g.db.execute("select ident, password from users where name=?",
                       [username])
    user_pairs = cur.fetchall()
    cur.close()
    # print ("user-pairs")
    # print (user_pairs)
    # Should do some error handling here but for now:
    if not user_pairs:
      return goto_anonymous_login(error="No such user as: " + username)
    (ident, actual_password) = user_pairs[0]
    if password == actual_password:
      user = User(ident)
      user.name = username
      # user.password = password
      flasklogin.login_user(user)
      return redirect(url_for('show_entries'))
             # redirect(request.args.get("next"))
    else:
      return abort(401)
  else:
    return goto_anonymous_login(error=None)

@app.route('/register', methods=['GET', 'POST'])
def register():
  """ Similar to the login handler, handles both GET and POST for
      registering a user. The GET displays the register form whilst the
      POST registers the new user and logs them in, assuming that the
      registration was valid (eg. an unused username).
  """
  if request.method == 'POST':
    # TODO: We should check if the user name is available:
    name = request.form['username']
    cur = g.db.execute("select name from users where name=?", [name])
    user_pairs = cur.fetchall()
    cur.close()
    if user_pairs:
      error = "Username: " + name + " already exists :("
      return render_template('register.html', error=error,
                           user=flasklogin.AnonymousUser())
    # TODO: Of course we shouldn't store the plaintext password but
    # instead some hash of it.
    password = request.form['password']
    cur = g.db.execute('insert into users (name, password) values (?, ?)',
                       [name, password])
    new_ident = cur.lastrowid
    g.db.commit()
    user = User(new_ident)
    user.name = name
    flasklogin.login_user(user)
    return redirect(url_for('show_entries'))
  else: # We assume it was a GET
    return render_template('register.html', error=None,
                           user=flasklogin.AnonymousUser())
  
@app.route("/logout")
@flasklogin.login_required
def logout():
  """Handler for logging the current user out."""
  flasklogin.logout_user()
  flash('You were logged out')
  return redirect(url_for('show_entries'))

def create_arguments_parser():
  """Create the command-line arguments parser"""
  description = "A web based interface to biopepa-to-latex"
  parser = argparse.ArgumentParser(add_help=True,
                                   description=description)
  parser.add_argument('--init-db', action='store_true',
                      help="Initialise the data base on startup")
  parser.add_argument('--auto-reload', action='store_true',
                      help="Use auto-reload of code on edit in folder")
  parser.add_argument('--debug', action='store_true',
                      help="Run in debug mode")

  return parser


def run():
  """The main method, starts the application running"""
  parser = create_arguments_parser() 
  arguments = parser.parse_args()
  if arguments.init_db:
    datastore.init_db(app)
  else:
    app.run(debug=arguments.debug, use_reloader=arguments.auto_reload)

if __name__ == '__main__':
  run()
