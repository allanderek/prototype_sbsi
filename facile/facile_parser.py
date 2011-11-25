"""
A module that implements a parser for the Bio-PEPA language
"""
import parcon
from parcon import Translate, SignificantLiteral

import biopepa.biopepa_parser as biopepa_parser
from biopepa.biopepa_parser import create_separated_by 
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
equation_section_syntax = parcon.OneOrMore(reaction_parser)


class InitialCondition(object):
  """Class representing an initial condition"""
  def __init__(self):
    self.name = None
    self.value = None
    self.units = None

def create_initial_condition(parse_result):
  """post parsing method for initial conditions"""
  init_cond = InitialCondition()
  init_cond.name = parse_result[0]
  init_cond.value = parse_result[1]
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


model_parser = (equation_section_syntax + 
                "INIT" +
                init_cond_section_syntax +
                parcon.End()
               )

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


