"""A module implementing a transformation from facile to sbml"""

import sys
import argparse
import parcon

import utils
import facile.facile_parser as facile_parser
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

def add_implicit_modifiers(reactions, facile_model):
  """Some models do not explicitly specify modifier species in reactions
     but instead simply use the species name in the kinetic law.
     For example where I would write
     A + B -> A + C; k1 = 0.1
     they might write instead
     B -> C; k1 = "A * B * 0.1"
     note that they are forced to write out the explicit rate equation
     since mass action does not work since A is not specified as a
     reactant/modifier.
     This function then inspects the kinetic law and adds any species
     mentioned as a modifier. Doing so allows the sbml model thus
     produced to be validated, however of course this is a questionable
     practice since it offers the modeller less checking.
  """
  # We would love to use the 'RateAnalyser from rate_checker_sbml here
  # however it only operates over sbml models because it literally
  # makes calls to element tags etc.
  # First we must determine all the species' names since we cannot
  # simply add any name since it might be a variable name.
  all_species_names = [ i.name for i in facile_model.initial_conditions ]
  # Now we get the names referenced by the kinetic law, this is
  # currently incomplete since it is only returning those names which
  # are literally used within the kinetic law, but we may have variables
  # set to something else, eg:
  # variable v = E + D
  # A -> B ; k1 = v;
  # would currently only return 'v' whereas we wish for it to return
  # E,D
  for reaction in reactions:
    used_names = reaction.kinetic_law.used_names()
    referenced_species = [ n for n in used_names 
                               if n in all_species_names ]
    for species in referenced_species:
      if (not reaction.is_reactant(species) and
          not reaction.is_modifier(species)):
        # Now this should NOT be added if the name itself is not even
        # a species
        reaction.add_modifier(species)


def translate_facile_model(facile_model, arguments):
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
    # Note that this could be 'tighter' basically this will still
    # add fMA to any reaction with no reactants/modifiers
    # We are somewhat relying on the Reaction's 'remove_rate_law_sugar'
    # method to detect this situation and do nothing.
    if equation.kinetic_law.no_implicit_fma:
      equation.kinetic_law = kinetic_expr
    else:
      equation.kinetic_law = apply_fma(kinetic_expr)
 
    # Similarly to the above, if the forward reaction has no products
    # and hence the reverse has no reactants then fMA is still added. 
    if equation.reverse_kinetic_law:
      rev_kin_expr = equation.reverse_kinetic_law.value_expr
      if equation.reverse_kinetic_law.no_implicit_fma:
        equation.reverse_kinetic_law = rev_kin_expr
      else:
        equation.reverse_kinetic_law = apply_fma(rev_kin_expr)

  # Now if the 'add-modifiers' flag is set, we must add any species
  # which are referenced in the kinetic law as a modifier. This is a
  # risky practice but will otherwise result in a non-validating sbml
  # model. 
  # Note doing this here, rather than in 'output_reaction' means that
  # the reaction is changed, should we ever decide to do something with
  # the model/reactions after outputting the sbml model.
  if arguments.add_modifiers:
    add_implicit_modifiers(all_equations, facile_model)

  sbml_model.reactions = all_equations

  def create_assign_init_cond(assign_rule):
    """From an assignment rule create an initial condition, this ensures
       that there will be a 'species' definition created for the
       assignment rule, which will come from a compound species/variable
       in facile, such as: ABC = AB + C;
    """
    init_cond = facile_parser.InitialCondition()
    init_cond.name = assign_rule.variable
    return init_cond
  assign_species = [ create_assign_init_cond(a) 
                       for a in facile_model.assign_rules ]
  sbml_model.component_defs = (facile_model.initial_conditions + 
                                assign_species )
  sbml_model.var_decs = facile_model.var_decs
  sbml_model.assign_rules = facile_model.assign_rules
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

  document = translate_facile_model(facile_model, arguments)
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
  parser.add_argument('--add-modifiers', action='store_true',
    help="Species mentioned in kinetic law added as modifiers")
  utils.add_output_file_arg(parser)
  arguments = parser.parse_args()

  for filename in arguments.filenames:
    translate_file(filename, arguments)


if __name__ == '__main__':
  main()


