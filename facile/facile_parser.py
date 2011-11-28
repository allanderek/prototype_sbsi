"""
A module that implements a parser for the Bio-PEPA language
"""
import parcon
from parcon import Translate, SignificantLiteral

import biopepa.biopepa_parser as biopepa_parser
from biopepa.biopepa_parser import create_separated_by 
import create_sbml
import outline_sbml


name_syntax = parcon.alphanum_word

list_of_species_syntax = create_separated_by(name_syntax, "+")
reaction_arrow_syntax = parcon.First(SignificantLiteral ("->"), 
                                     SignificantLiteral ("<-"),
                                     SignificantLiteral ("<->"))

rate_law_syntax = name_syntax + "=" + biopepa_parser.expr + ";"

def create_reaction(parse_result):
  """The post-parse action for the reaction parser"""
  reaction = outline_sbml.Reaction("reaction") 
  reaction.reactants = [ outline_sbml.ReactionParticipant(n, 1) 
                            for n in parse_result[0] ]
  reaction.products =  [ outline_sbml.ReactionParticipant(n, 1) 
                            for n in parse_result[2] ]
  reaction.kinetic_law = parse_result[4]
  return reaction
  

reaction_syntax = (list_of_species_syntax +
                   reaction_arrow_syntax +
                   list_of_species_syntax +
                   ";"
                   + rate_law_syntax
                  )
reaction_parser = Translate(reaction_syntax, create_reaction)

# No reason why parcon.number couldn't be an expression
var_dec_syntax = "variable" + name_syntax + "=" + parcon.number
def create_var_dec(parse_result):
  """post-parsing method for variable declarations"""
  name = parse_result[0]
  expression = create_sbml.NumExpression(parse_result[1])
  var_dec = create_sbml.VariableDeclaration(name, expression)
  return var_dec
var_dec_parser = Translate(var_dec_syntax, create_var_dec)

reaction_or_var_dec_parser = parcon.First(var_dec_parser, reaction_parser)
equation_section_syntax = parcon.OneOrMore(reaction_or_var_dec_parser)


class InitialCondition(object):
  """Class representing an initial condition"""
  def __init__(self):
    self.name = None
    self.initial_amount = None
    self.units = None
    self.location = None

def create_initial_condition(parse_result):
  """post parsing method for initial conditions"""
  init_cond = InitialCondition()
  init_cond.name = parse_result[0]
  init_cond.initial_amount = parse_result[1]
  init_cond.units = parse_result[2]
  return init_cond
  
init_cond_units_syntax = parcon.First(SignificantLiteral ("N"),
                                      SignificantLiteral ("M"),
                                      SignificantLiteral ("uM")
                                     )
initial_condition_syntax = (name_syntax +
                            "=" +
                            parcon.number +
                            init_cond_units_syntax
                           )
initial_condition_parser = Translate(initial_condition_syntax,
                                     create_initial_condition)

init_cond_section_syntax = parcon.ZeroOrMore(initial_condition_parser)


class FacileModel(object):
  """A class for representing a parsed facile model"""
  def __init__(self):
    self.equations = None
    self.var_decs = None
    self.initial_conditions = None

model_syntax = (equation_section_syntax + 
                "INIT" +
                init_cond_section_syntax +
                parcon.End()
               )
def create_model(parse_result):
  """post parsing method for an entire facile model"""
  facile_model = FacileModel()
  eqn_section = parse_result[0]
  facile_model.equations = [ s for s in eqn_section 
                                 if isinstance(s, outline_sbml.Reaction) ]
  var_decs = [ s for s in eqn_section
                   if isinstance(s, create_sbml.VariableDeclaration)
             ]
  facile_model.var_decs = var_decs
  facile_model.initial_conditions = parse_result[1]
  return facile_model

model_parser = Translate(model_syntax, create_model)

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


