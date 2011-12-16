"""A module to convert SBML models into facile models where possible"""
import xml.dom.minidom
import argparse
import sys

import utils
from outline_sbml import format_math_element
import outline_sbml


def output_reaction(output_file, reaction):
  """Output the facile format of an sbml reaction"""
  def names_of_participants(participants):
    """ We would like to get the names of the participants, ie
        reactants, modifiers and products by simply mapping the
        participants to their names. But some of them have different
        stoichiometeries, so this turns a list of participants into
        a list of names, for example if A has a stoich of three then
        it generates A,A,A, such that it can be formatted as: A + A + A.
    """ 
    return [ item for p in participants
                  for item in [ p.name ] * p.stoich ]
  reactant_names = names_of_participants(reaction.reactants)
  # modifiers shouldn't have stoichiometries
  modifier_names = [ m.name for m in reaction.modifiers ]
  product_names = names_of_participants(reaction.products)
  left = " + ".join(reactant_names + modifier_names)
  right = " + ".join(product_names + modifier_names)
  output_file.write(left) 
  output_file.write(" -> ")
  output_file.write(right)
  output_file.write(" ; ")

  output_file.write(format_math_element(reaction.kinetic_law))
  output_file.write("\n")


def convert_file(filename, arguments):
  """Given an sbml file convert it into a corresponding facile model"""
  dom = xml.dom.minidom.parse(filename)
  model = dom.getElementsByTagName("model")[0]

  if arguments.output_file:
    eqn_filename = arguments.output_file
  else:
    eqn_filename = utils.change_filename_ext(filename, ".eqn")
  if eqn_filename == "stdout":  
    output_file = sys.stdout
  else:
    output_file = open(eqn_filename, "w")

  # We also need to do the same for variables and parameters
  parameters = outline_sbml.get_list_of_parameters(model)
  for param in parameters:
    output_file.write("variable " +
                      param.name +
                      " = " +
                      param.value +
                      " ;"
                      )
    output_file.write("\n")
  # Reactions
  reactions = outline_sbml.get_list_of_reactions(model)
  for reaction in reactions:
    output_reaction(output_file, reaction)
 
  output_file.write("\nINIT\n")
  species = outline_sbml.get_list_of_species(model)
  init_assigns = outline_sbml.get_list_of_init_assigns(model)
  for component in species:
    init_string = component.initial_amount
    for init_assign in init_assigns:
      if init_assign.variable == component.name:
        init_string = format_math_element(init_assign.expression)
    output_file.write(component.name +
                      " = " +
                      init_string +
                      " N;\n")

  if eqn_filename != "stdout":
    output_file.close()


def run():
  """Perform the banalities of command-line argument processing
     and then do the actual work on each argument file
  """
  description = "Convert SBML models into facile model files"""
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
