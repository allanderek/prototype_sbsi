import argparse
import sqlite3
from flask import Flask, request, session, g, redirect, url_for, \
     abort, render_template, flash
from contextlib import closing


# configuration
DATABASE = '/tmp/flaskr.db'
DEBUG = True
SECRET_KEY = 'development key'
USERNAME = 'admin'
PASSWORD = 'default'

# create our little application :)
app = Flask(__name__)
app.config.from_object(__name__)

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
def show_entries():
  fields = "id, title, text, latex"
  cur = g.db.execute("select " + fields + " from entries order by id desc")
  entries = [ dict(id=row[0], title=row[1], text=row[2], latex=row[3])
              for row in cur.fetchall()]
  return render_template('show_entries.html', entries=entries)

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

@app.route('/login', methods=['GET', 'POST'])
def login():
  error = None
  if request.method == 'POST':
    if request.form['username'] != app.config['USERNAME']:
      error = 'Invalid username'
    elif request.form['password'] != app.config['PASSWORD']:
      error = 'Invalid password'
    else:
      session['logged_in'] = True
      flash('You were logged in')
      return redirect(url_for('show_entries'))
  return render_template('login.html', error=error)

@app.route('/logout')
def logout():
    session.pop('logged_in', None)
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
