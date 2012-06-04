"""A module implementing a transformation from biopepa to xml"""

import sys
import argparse

import utils
import biopepa.biopepa_parser as biopepa_parser
import sbml_ast

class Reaction(sbml_ast.Reaction):
  """A class representing reactions within a Bio-PEPA model. These
     reactions have been derived from the rate laws and the component
     definitions (and their associated reaction behaviours)
  """
  def __init__(self, name):
    super(Reaction, self).__init__(name)

  @property
  def rate(self):
    """The Bio-PEPA code essentially called the kinetic_law part
       the rate part so we set up a property to mimic that.
    """
    return self.kinetic_law

  @rate.setter
  def rate(self, rate):
    """The setter for the rate property"""
    self.kinetic_law = rate


  def add_behaviour(self, comp_name, behaviour):
    """Add the given behaviour to the reaction, essentially this means
       we are adding a component as a participant (reactant, product
       or modifier) to the reaction.
    """
    stoichimetries = behaviour.stoichiometry
    operators = behaviour.operator
    if len(stoichimetries) != len(operators):
      print ("Must have the same number of operators as stoichimetries")
      sys.exit(1)

    # Not sure if this is really the best way to do this, this means
    # we will be adding multiple behaviours for the same component.
    # eg, we'll have A + B -> B + B, where as arguably we should retain
    # the knowledge, perhaps we should just have the operator <> or ><
    for stoich, operator in zip(stoichimetries, operators):
      participant = sbml_ast.ReactionParticipant(comp_name)
      participant.stoich = stoich
      if operator == ">>" :
        self.products.append(participant)
      elif operator == "<<": 
        self.reactants.append(participant)
      elif operator == "(+)":
        self.modifiers.append(participant)
      elif operator == "(-)":
        self.modifiers.append(participant)
      elif operator == "(.)":
        self.modifiers.append(participant)
      else:
        print ("Unrecognised behaviour operator: " + operator)
        sys.exit(1)


def build_reaction_dictionary(components, rate_definitions):
  """From a list of components build a reaction dictionary which maps
     reaction names to reaction representations"""
  reaction_dictionary = dict()
  for rate_def in rate_definitions:
    name = rate_def.get_name()
    reaction = Reaction(name)
    reaction.rate = rate_def.get_rate()
    reaction_dictionary[name] = reaction

  def have_same_locations(reaction, reaction_behaviour):
    """returns true if the reaction has the same location as
       the reaction behaviour
    """
    return reaction.location == reaction_behaviour.location

  def new_reaction_for_behaviour(behaviour):
    """If there is no existing reaction to which to add a given
       behaviour then call this method to create a new such reaction
       and add it to the reaction dictionary
    """
    reaction = Reaction(b_name)
    reaction_dictionary[b_name] = reaction
    if behaviour.location:
      reaction.location = behaviour.location
    return reaction

  for component_def in components:
    for behaviour in component_def.get_behaviours():
      b_name = behaviour.get_name() 
      # So basically we add the current behaviour to a reaction
      # if we have already created a reaction with the given name
      # at the given location (if any) of the behaviour. Otherwise
      # we create a new reaction and add it to the dictionary, but either
      # way we add the current behaviour to some reaction.

      # Ach this totally won't work if we have the same reaction name
      # with different locations because we will overwrite the mapping
      # with the second (and subsequent) ones we encounter.

      if b_name in reaction_dictionary:
        reaction = reaction_dictionary[b_name]
        if not have_same_locations(reaction, behaviour):
          # So there exists a reaction with the correct name, but not
          # at the correct location
          reaction = new_reaction_for_behaviour(behaviour)
      else:
        # There is no reaction with the correct name, regardless of
        # the location of such a reaction.
        reaction = new_reaction_for_behaviour(behaviour)
      # Whether we have created a new reaction or we are adding
      # behaviour to an existing one we must add the current
      # behaviour to the reaction.
      reaction.add_behaviour(component_def.get_name(), behaviour)
  return reaction_dictionary


def translate_biopepa_model(biopepa_model):
  """Translate the parsed Bio-PEPA model into an SBML xml document"""
  
  react_dict = build_reaction_dictionary(biopepa_model.component_defs,
                                         biopepa_model.rate_defs)
  reactions = react_dict.values()

  sbml_model = sbml_ast.SBMLModel()
  sbml_model.reactions = reactions
  sbml_model.var_decs = biopepa_model.var_decs
  sbml_model.component_defs = biopepa_model.system_equation

  return sbml_model.create_sbml_document()

def convert_file(filename, arguments):
  """Parse in a Bio-PEPA file, translate to SBML and create an
     SBML file which should be the translated Bio-PEPA model.
  """
  model_file = open(filename, "r")
  parse_result = biopepa_parser.parse_model_file_exit_on_error(model_file)
  model_file.close()

  document = translate_biopepa_model(parse_result)
  sbml_ast.output_to_sbml_file(filename, arguments, document)

def main():
  """A simple main function to parse in the arguments as
     Bio-PEPA model files"""
  description = "Translates Bio-PEPA model files to SBML"
  parser = argparse.ArgumentParser(add_help=True,
                                   description=description)
   # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="A Bio-PEPA model file to translate")
  utils.add_output_file_arg(parser)
  arguments = parser.parse_args()

  for filename in arguments.filenames:
    convert_file(filename, arguments)


if __name__ == '__main__':
  main()


