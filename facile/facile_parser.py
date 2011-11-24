"""
A module that implements a parser for the Bio-PEPA language
"""
import parcon
from parcon import Forward, InfixExpr, Translate, Optional, ZeroOrMore


model_parser = parcon.alphanum_word + parcon.End()

def parse_model(model_source):
 """Takes in the string which represents the source of the model.
    You can instead call 'parse_model_file' but this is useful if
    you have the source already, for example perhaps as part of a
    web or gui application
 """
 return model_parser.parse_string(model_source)

def parse_model_file(model_file):
  """Given a model file (handle, not filename), parse the contents
     of the file as a facile model
  """
  parse_result = model_parser.parse_string(model_file.read())
  return parse_result


