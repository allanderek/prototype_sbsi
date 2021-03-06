"""A module to convert SBML models into Bio-PEPA models where possible"""
import argparse
import sys

import utils
from outline_sbml import format_math_element
import outline_sbml
import biopepa.biopepa_parser as biopepa


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

    # Essentially we could check if it is an assignment rule, when
    # a species is defined by an assignment rule essentially it is
    # a variable, rather than actual species in which case we do not
    # want it in the system equation, we should of course though have
    # it in the biopepa file as a dynamic variable
    if init_string != None and init_string != "":
      output_file.write(prefix)
      prefix = "<*> "
      output_file.write(component.ident)
      output_file.write("[" + init_string + "]")
      output_file.write("\n")


def calculate_component_defs(reactions):
  """From the set of reactions calculate what the set of Bio-PEPA
     component definitions should be
  """
  def_map = dict()
  for reaction in reactions:
    # This is wrong we should check if the same named
    # component is listed as both a reactant and product / modifier.
    # And essentially just do the correct thing, which would be to have
    # a proper list of operators and stoichiometries rather than the
    # singleton lists used here.
    def add_behaviour(participant, operator):
      """Creates a behaviour from the participant and adds it to the
         corresponding component definition in the def map. If there is no
         current definition of this process in the def map then a new one
         created and added.
      """
      name = participant.name 
      if name in def_map:
        comp_def = def_map[name]
      else:
        comp_def = biopepa.ComponentDefinition(name, [])
        def_map[name] = comp_def
    
      behaviour = biopepa.Behaviour(reaction.name, [ operator ])
      behaviour.stoichiometry = [ participant.stoich ]
      comp_def.behaviours.append(behaviour) 

    for participant in reaction.reactants:
      add_behaviour(participant, "<<")
    for participant in reaction.products:
      add_behaviour(participant, ">>")
    for participant in reaction.modifiers:
      add_behaviour(participant, "(.)")


  return def_map.values()

def translate_file(filename, arguments):
  """Given an sbml file convert it into a corresponding Bio-PEPA model"""
  model = outline_sbml.get_model_from_sbml_file(filename)

  if arguments.output_file:
    biopepa_filename = arguments.output_file
  else:
    biopepa_filename = utils.change_filename_ext(filename, ".biopepa")
  if biopepa_filename == "stdout":  
    output_file = sys.stdout
  else:
    output_file = open(biopepa_filename, "w")

  output_file.write("\n\n// Parameters\n")
  # We also need to do the same for variables and parameters
  parameters = outline_sbml.get_list_of_parameters(model)
  for param in parameters:
    output_file.write(" ".join([param.name, "=", param.value, ";\n"]))
    output_file.write("\n")

  # Reactions
  reactions = outline_sbml.get_list_of_reactions(model)

  output_file.write("\n\n// Rate Definitions\n")
  for reaction in reactions:
    output_file.write(reaction.name + " = [ ")
    output_file.write(format_math_element(reaction.kinetic_law))
    output_file.write(" ] ;\n")

  output_file.write("\n\n// Component Definitions\n")

  component_defs = calculate_component_defs(reactions)
  for component_def in component_defs:
    output_file.write(component_def.show_definition())
    output_file.write("\n")
 
  output_file.write("\n// System Equation\n")
  species = outline_sbml.get_list_of_species(model)
  init_assigns = outline_sbml.get_list_of_init_assigns(model)

  output_system_equation(output_file, species, init_assigns)

  if biopepa_filename != "stdout":
    output_file.close()


def run():
  """ Command-line argument processing and then real work.
  """
  description = "Convert SBML models into Bio-PEPA model files"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="an SBML file to translate to Bio-PEPA")
  utils.add_output_file_arg(parser)
  arguments = parser.parse_args()

  for filename in arguments.filenames:
    translate_file(filename, arguments)

if __name__ == "__main__":
  run()
