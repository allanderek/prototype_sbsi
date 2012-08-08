""" A web application for translating Bio-PEPA files into LaTeX"""
import argparse
import sqlite3
from flask import Flask, request, session, g, redirect, url_for, \
     abort, render_template, flash
from contextlib import closing

import flask.ext.login as flasklogin

# configuration
DATABASE = '/tmp/flaskr.db'
DEBUG = True
SECRET_KEY = 'development key'
USERNAME = 'admin'
PASSWORD = 'default'

# create our little application :)
app = Flask(__name__)
app.config.from_object(__name__)

# Create the login code for our little application
login_manager = flasklogin.LoginManager()
login_manager.init_app(app)
login_manager.login_view = "login"


# silly user model
class User(flasklogin.UserMixin):
  def __init__(self, id):
    self.id = id
    self.name = None

  def __repr__(self):
    return "%d/%s/%s" % (self.id, self.name, self.password)

# callback to reload the user object        
@login_manager.user_loader
def load_user(userid):
  g.db = connect_db()
  cur = g.db.execute("select name from users where id=?", [userid])
  users = cur.fetchall()
  if not users:
    return flasklogin.AnonymousUser()

  user_tuple = users[0]
  user = User(userid)
  user.name = user_tuple[0]
  
  return user

def connect_db():
  return sqlite3.connect(app.config['DATABASE'])

def init_db():
  with closing(connect_db()) as db:
    with app.open_resource('schema.sql') as f:
      db.cursor().executescript(f.read())
    db.commit()


@app.before_request
def before_request():
  g.db = connect_db()

@app.teardown_request
def teardown_request(exception):
  g.db.close()


@app.route('/')
@flasklogin.login_required
def show_entries():
  # First get the current user:
  current_user = flasklogin.current_user 

  fields = "id, title, text, latex"
  cur = g.db.execute("select " + fields + " from entries order by id desc")
  entries = [ dict(id=row[0], title=row[1], text=row[2], latex=row[3])
              for row in cur.fetchall()]
  return render_template('show_entries.html',
                         entries=entries,
                         user=current_user)

@app.route('/add', methods=['POST'])
def add_entry():
  if not session.get('logged_in'):
    abort(401)
  g.db.execute('insert into entries (title, text) values (?, ?)',
               [request.form['title'], request.form['text']])
  g.db.commit()
  flash('New entry was successfully posted')
  return redirect(url_for('show_entries'))

import biopepa.biopepa_to_latex as biopepa_to_latex
def convert_model(source):
  return biopepa_to_latex.convert_source(source)

@app.route('/convert', methods=['POST'])
def convert_entry():
  if not session.get('logged_in'):
    abort(401)
  convert_id = request.form["convert_id"]
  cur = g.db.execute('select text from entries where id=?', [convert_id])
  texts = cur.fetchall()
  # Should do some error handling here but for now:
  text = texts[0]
  new_text = convert_model(text[0])
  # g.db.execute("update entries set latex=?", [new_text])
  g.db.execute('update entries set latex=? where id=?',
               [new_text, convert_id])
  g.db.commit()
  flash ("You wanted to convert: " + str(convert_id))
  return redirect(url_for('show_entries'))

@app.route('/delete', methods=['POST'])
def delete_entry():
  if not session.get('logged_in'):
    abort(401)
  convert_id = request.form["convert_id"]
  g.db.execute('delete from entries where id=?', [convert_id])
  g.db.commit()
  flash ("You just deleted the id: " + convert_id)
  return redirect(url_for('show_entries'))

@app.route('/no-store-biopepa-latex', methods=['GET', 'POST'])
def no_store_biopepa_latex():
  error = None
  model_source = None
  latex = None
  if request.method == 'POST':
    model_source = request.form['modelsource']
    latex = convert_model(model_source)
  current_user = flasklogin.current_user 
  return render_template('raw-to-latex.html', 
                         source=model_source,
                         latex=latex, error=error,
                         user=current_user)

@app.route('/login', methods=['GET', 'POST'])
def login():
  if request.method == 'POST':
    username = request.form['username']
    password = request.form['password']
    cur = g.db.execute("select id, password from users where name=?",
                       [username])
    user_pairs = cur.fetchall()
    cur.close()
    # print ("user-pairs")
    # print (user_pairs)
    # Should do some error handling here but for now:
    (id, actual_password) = user_pairs[0]
    if password == actual_password:
      user = User(id)
      user.name = username
      # user.password = password
      flasklogin.login_user(user)
      return redirect(url_for('show_entries'))
             # redirect(request.args.get("next"))
    else:
      return abort(401)
  else:
    return render_template('login.html', error=None,
                           user=flasklogin.AnonymousUser())

@app.route('/register', methods=['GET', 'POST'])
def register():
  if request.method == 'POST':
    # TODO: We should check if the user name is available:
    name = request.form['username']
    # TODO: Of course we shouldn't store the plaintext password but
    # instead some hash of it.
    password = request.form['password']
    g.db.execute('insert into users (name, password) values (?, ?)',
               [name, password])
    g.db.commit()
    return redirect(url_for('show_entries'))
  else: # We assume it was a GET
    return render_template('register.html', error=None,
                           user=flasklogin.AnonymousUser())
  
@app.route("/logout")
@flasklogin.login_required
def logout():
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

  return parser


def run():
  parser = create_arguments_parser() 
  arguments = parser.parse_args()
  if arguments.init_db:
    init_db()

  app.run()

if __name__ == '__main__':
  run()
