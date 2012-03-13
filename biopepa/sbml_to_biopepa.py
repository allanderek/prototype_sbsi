"""A module to convert SBML models into Bio-PEPA models where possible"""
import xml.dom.minidom
import argparse
import sys

import utils
from outline_sbml import format_math_element
import outline_sbml
import biopepa_parser as biopepa


def output_system_equation(output_file, species, init_assigns):
  """Output the species and their initial concentrations as
     the system equation
  """

  # This prefix is a little trick to put <*> in between
  # all the component initialisers. We could do this via
  # a "<*>".join(...) but this way is a bit more efficient.
  prefix = ""
  for component in species:
    init_string = component.initial_amount
    if not init_string:
      # This is NOT correct, but I'm not quite sure how to deal with
      # concentrations as of yet.
      init_string = component.initial_conc
    for init_assign in init_assigns:
      if init_assign.variable == component.name:
        init_string = format_math_element(init_assign.expression)
        break
    output_file.write(prefix)
    prefix = "<*> "
    output_file.write(component.name)
    output_file.write("[" + init_string + "]")
    output_file.write("\n")


def calculate_component_defs(reactions):
  """From the set of reactions calculate what the set of Bio-PEPA
     component definitions should be
  """
  def_map = dict()
  for reaction in reactions:
    def add_behaviour(participant, operator):
      name = participant.name 
      if name in def_map:
        comp_def = def_map[name]
      else:
        comp_def = biopepa.ComponentDefinition(name, [])
        def_map[name] = comp_def
    
      behaviour = biopepa.Behaviour(reaction.name, operator)
      behaviour.stoichiometry = participant.stoich 
      comp_def.behaviours.append(behaviour) 

    for participant in reaction.reactants:
      add_behaviour(participant, "<<")
    for participant in reaction.products:
      add_behaviour(participant, ">>")
    for participant in reaction.modifiers:
      add_behaviour(participant, "(.)")


  return def_map.values()

def convert_file(filename, arguments):
  """Given an sbml file convert it into a corresponding Bio-PEPA model"""
  dom = xml.dom.minidom.parse(filename)
  model = dom.getElementsByTagName("model")[0]

  if arguments.output_file:
    eqn_filename = arguments.output_file
  else:
    eqn_filename = utils.change_filename_ext(filename, ".biopepa")
  if eqn_filename == "stdout":  
    output_file = sys.stdout
  else:
    output_file = open(eqn_filename, "w")

  output_file.write("\n\n// Parameters\n")
  # We also need to do the same for variables and parameters
  parameters = outline_sbml.get_list_of_parameters(model)
  for param in parameters:
    output_file.write(param.name +
                      " = " +
                      param.value +
                      " ;"
                      )
    output_file.write("\n")

  # Reactions
  reactions = outline_sbml.get_list_of_reactions(model)

  output_file.write("\n\n// Rate Definitions\n")
  for reaction in reactions:
    output_file.write(reaction.name + " = [ ")
    output_file.write(format_math_element(reaction.kinetic_law))
    output_file.write(reaction.name + " ] ;\n")

  output_file.write("\n\n// Component Definitions\n")

  component_defs = calculate_component_defs(reactions)
  for component_def in component_defs:
    output_file.write(component_def.show_definition())
    output_file.write("\n")
 
  output_file.write("\n// System Equation\n")
  species = outline_sbml.get_list_of_species(model)
  init_assigns = outline_sbml.get_list_of_init_assigns(model)

  output_system_equation(output_file, species, init_assigns)

  if eqn_filename != "stdout":
    output_file.close()


def run():
  """ Command-line argument processing and then real work.
  """
  description = "Convert SBML models into Bio-PEPA model files"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="an SBML file to translate")
  utils.add_output_file_arg(parser)
  arguments = parser.parse_args()

  for filename in arguments.filenames:
    convert_file(filename, arguments)

if __name__ == "__main__":
  run()
