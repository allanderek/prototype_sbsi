"""A module to convert SBML models into facile models where possible"""
import xml.dom.minidom
import argparse
import sys

from outline_sbml import format_math_element
import outline_sbml


def convert_file(filename):
  """Given an sbml file convert it into a corresponding facile model"""
  dom = xml.dom.minidom.parse(filename)
  model = dom.getElementsByTagName("model")[0]

  output_file = sys.stdout

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
    # We should look at the stoichiometry and possibly do P + P
    # if the stoich is two for example.
    reactant_names = [ r.name for r in reaction.reactants ]
    modifier_names = [ m.name for m in reaction.modifiers ]
    product_names = [ p.name for p in reaction.products ]
    left = " + ".join(reactant_names + modifier_names)
    right = " + ".join(product_names)
    output_file.write(left) 
    output_file.write(" -> ")
    output_file.write(right)
    output_file.write(" ; ")

    output_file.write(format_math_element(reaction.kinetic_law))
    output_file.write("\n")
  
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


def run():
  """Perform the banalities of command-line argument processing
     and then do the actual work on each argument file
  """
  description = "Convert SBML models into facile model files"""
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="an SBML file to translate")

  arguments = parser.parse_args()

  for filename in arguments.filenames:
    convert_file(filename)

if __name__ == "__main__":
  run()
