"""A module implementing a transformation from facile to sbml"""

import sys
import argparse
import parcon

import utils
import facile_parser
import sbml_ast

def uniquely_name_reactions(reactions):
  """Give a unique name to each reaction, this is currently rather
     dumb and does not depend the current name for the reaction"""
  name_dictionary = dict()
  for reaction in reactions:
    if reaction.name in name_dictionary:
      number = name_dictionary[reaction.name]
      name_dictionary[reaction.name] = number + 1
      reaction.name = reaction.name + "_" + str(number)
    else:
      name_dictionary[reaction.name] = 1

def translate_facile_model(facile_model):
  """Translate the parsed facile model into an SBML xml document"""
  sbml_model = sbml_ast.SBMLModel()
  equations = facile_model.equations
  for equation in equations:
    equation.canonicalise_participants()
  uniquely_name_reactions(equations)
  reverse_equations = [ r.reverse_reaction() 
                          for r in equations
                            if r.reverse_kinetic_law ]
  all_equations = equations + reverse_equations
  # So first of all this should somehow check that the original
  # rate law was not given surrounded by quotes (probably a simple)
  # flag within 'reaction'. Additionally, again I sort of think this
  # belongs in facile_to_sbml.py
  def apply_fma(kinetic_law):
    """Apply the fMA method to the given rate law as is implied
       in the syntax of facile models"""
    return sbml_ast.ApplyExpression("fMA", [kinetic_law])
  # Note that this changes all equation kinetic laws from RateLaw
  # to simple Expression which is expected by the formatter
  # for the sbml ast.
  for equation in all_equations:
    kinetic_expr = equation.kinetic_law.value_expr
    if equation.kinetic_law.no_implicit_fma:
      equation.kinetic_law = kinetic_expr
    else:
      equation.kinetic_law = apply_fma(kinetic_expr)
    if equation.reverse_kinetic_law:
      rev_kin_expr = equation.reverse_kinetic_law.value_expr
      if equation.reverse_kinetic_law.no_implicit_fma:
        equation.reverse_kinetic_law = rev_kin_expr
      else:
        equation.reverse_kinetic_law = apply_fma(rev_kin_expr)

  sbml_model.reactions = all_equations

  sbml_model.component_defs = facile_model.initial_conditions
  sbml_model.var_decs = facile_model.var_decs
  document = sbml_model.create_sbml_document()
  return document

def translate_file(filename, arguments):
  """Translate a facile file into sbml"""
  model_file = open(filename, "r")
  try:
    facile_model = facile_parser.parse_model_file(model_file)
  except parcon.ParseException as parse_except:
    print parse_except
    sys.exit(1)
  model_file.close()

  document = translate_facile_model(facile_model)
  sbml_ast.output_to_sbml_file(filename, arguments, document)


def main():
  """A simple main function to parse in the arguments
     as facile model file
  """
  description = "Translate facile model file(s) into sbml"
  parser = argparse.ArgumentParser(add_help=True,
                                   description=description)
   # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="A facile model file to translate")
  utils.add_output_file_arg(parser)
  arguments = parser.parse_args()

  for filename in arguments.filenames:
    translate_file(filename, arguments)


if __name__ == '__main__':
  main()


