"""
A module that implements a parser for the Bio-PEPA language
"""
import parcon
from parcon import Translate, SignificantLiteral, OneOrMore, First

import biopepa.biopepa_parser as biopepa_parser
from biopepa.biopepa_parser import create_separated_by 
import sbml_ast

name_syntax = parcon.Word(parcon.alphanum_chars + "_", 
                          init_chars=parcon.alpha_chars)
empty_species_list = Translate(parcon.Literal("null"), lambda x : [])
list_of_species_syntax = First (empty_species_list, 
                                create_separated_by(name_syntax, "+"))

# Using Longest might not be the fastest way to do this
# another possibility is to use 'First' but ensure that the
# alternatives are ordered such that the longest will be the first
# one matched if possible.
reaction_arrow_syntax = parcon.Longest(SignificantLiteral ("->"), 
                                       SignificantLiteral ("=>"),
                                       SignificantLiteral ("<-"),
                                       SignificantLiteral ("<="),
                                       SignificantLiteral ("<->"),
                                       SignificantLiteral ("<=>"))


# In facile a rate law can be written with an optional defining name
# as in:
# A -> B; k1 = 0.1 ;
# Additionally the expression part (whether the name is there or not)
# can be surrounded in quotes, this is to indicate that the implicit
# 'fMA' function should not be applied to the rate law.
# So we parse the expression part as a rate law leaving the 'name' part
# of the rate law as undefined. A rate law has an optional name and if
# this is applied then we simply update the 'rate_law' returned from
# the expression parser to include the given name.
class RateLaw(object):
  """A class representing the abstract syntax of a rate law in facile"""
  def __init__(self):
    self.name = None
    self.value_expr = None
    self.no_implicit_fma = False

def rate_law_expr (expr):
  """Parse a rate law which is not surrouned by quotes and therefore
     does not implicitly require fMA added around it
  """
  rate_law = RateLaw()
  rate_law.value_expr = expr
  return rate_law
naked_expression = Translate(biopepa_parser.expr, rate_law_expr)
def rate_law_explicit(expr):
  """Parse an explicit rate law surrounded in quotes"""
  rate_law = rate_law_expr(expr) 
  rate_law.no_implicit_fma = True
  return rate_law
explicit_expression = Translate("\"" + biopepa_parser.expr + "\"",
                                rate_law_explicit)
rate_law_syntax = (parcon.Optional(name_syntax + "=") + 
                   First (naked_expression, explicit_expression) +
                   ";"
                  )
def add_rate_law_name(parse_result):
  """The post-parsing method for rate laws"""
  if isinstance(parse_result, RateLaw):
    rate_law = parse_result
  else:
    rate_law = parse_result[1]
    rate_law.name = parse_result[0]

  return rate_law

rate_law_parser = rate_law_syntax[add_rate_law_name]
reaction_core_syntax = (list_of_species_syntax +
                        reaction_arrow_syntax +
                        list_of_species_syntax +
                        ";"
                       )
def decide_how_many_rate_laws(parse_result):
  """The 'binding' method for the reaction parser, once we have parsed
     the reaction we know whether the reaction arrow is uni-directional
     or bi-directional and hence whether to expect 1 or 2 rate laws.
     This returns the appropriate parser
  """
  reaction = sbml_ast.Reaction("reaction") 
  reaction.reactants = [ sbml_ast.ReactionParticipant(n, 1) 
                            for n in parse_result[0] ]
  reaction.products =  [ sbml_ast.ReactionParticipant(n, 1) 
                            for n in parse_result[2] ]

  def add_single_rate_law(rate_law):
    """The post parsing method for a uni-directional reaction"""
    reaction.kinetic_law = rate_law
    return reaction
  def add_double_rate_law(rate_laws):
    """The post parsing function for a reversed reaction"""
    reaction.kinetic_law = rate_laws[0]
    reaction.reverse_kinetic_law = rate_laws[1]
    return reaction

  if parse_result[1] in [ "<->", "<=>" ]:
    return Translate(rate_law_parser + rate_law_parser,
                     add_double_rate_law)
  else:
    return Translate(rate_law_parser, add_single_rate_law)

reaction_syntax = parcon.Bind (reaction_core_syntax,
                               decide_how_many_rate_laws)

def create_reaction(parse_result):
  """The post-parse action for the reaction parser"""
  reaction = sbml_ast.Reaction("reaction") 
  reaction.reactants = [ sbml_ast.ReactionParticipant(n, 1) 
                            for n in parse_result[0] ]
  reaction.products =  [ sbml_ast.ReactionParticipant(n, 1) 
                            for n in parse_result[2] ]
  rate_law = parse_result[3]
  reaction.kinetic_law = rate_law.value_expr
  return reaction

reaction_parser = reaction_syntax 
                  # Translate(reaction_syntax, create_reaction)

var_dec_syntax = ("variable" + 
                  name_syntax + 
                  "=" + 
                  biopepa_parser.expr
                 )
def create_var_dec(parse_result):
  """post-parsing method for variable declarations"""
  name = parse_result[0]
  expression = parse_result[1]
  var_dec = sbml_ast.VariableDeclaration(name, expression)
  return var_dec
var_dec_parser = Translate(var_dec_syntax, create_var_dec)

reaction_or_var_dec_parser = First(var_dec_parser, reaction_parser)
equation_section_syntax = OneOrMore(reaction_or_var_dec_parser)


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
  
init_cond_units_syntax = First(SignificantLiteral ("N"),
                               SignificantLiteral ("M"),
                               SignificantLiteral ("uM")
                              )
initial_condition_syntax = (name_syntax +
                            "=" +
                            # No reason why this can't be an expression?
                            biopepa_parser.scientific_number +
                            init_cond_units_syntax +
                            ";"
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
                                 if isinstance(s, sbml_ast.Reaction) ]
  var_decs = [ s for s in eqn_section
                   if isinstance(s, sbml_ast.VariableDeclaration)
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
  newline_syntax = parcon.CharIn("\n\r")
  not_new_line = parcon.Except(parcon.AnyChar(), newline_syntax)
  line_comment = "#" + parcon.ZeroOrMore(not_new_line)
  comments_and_whitespace = OneOrMore(First(parcon.Whitespace(),
                                            line_comment))
  return model_parser.parse_string(model_source,
                                   whitespace=comments_and_whitespace)

def parse_model_file(model_file):
  """Given a model file (handle, not filename), parse the contents
     of the file as a facile model
  """
  return parse_model(model_file.read())


