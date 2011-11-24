""" A simple script to replace the values of parameters within an
    sbml model file"""
import xml.dom.minidom
import argparse

import parameters
import utils


def parameterise_model(model, dictionary):
  """Given the model and dictionary of new parameter values,
     parameterise the model"""
  lists_of_parameters = model.getElementsByTagName("listOfParameters")
  for list_of_params in lists_of_parameters:
    parameter_elements = list_of_params.getElementsByTagName("parameter")
    for parameter_element in parameter_elements:
      param_name = parameter_element.getAttribute("id")
      if param_name in dictionary:
        new_value = dictionary[param_name]
        parameter_element.setAttribute("value", str(new_value))
            

def parameterise_model_file (filename, dictionary):
  """Parse in a file as an SBML model, and re-parameterise it based
     on the given dictionary"""
  dom = xml.dom.minidom.parse(filename)
  for model in dom.getElementsByTagName("model"):
    parameterise_model(model, dictionary)
  # Print out the modified dom model.
  print(dom.toxml("UTF-8"))

  
def run():
  """The main entry point, parameterise SBML model files given
     on the command line
  """
  description = "Parameterise an SBML model based on a given param file"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="input files: parameters and sbml model files")

  arguments = parser.parse_args()

  sbml_extentions = [ ".xml", ".sbml" ]
  param_files = [ x for x in arguments.filenames
                        if not utils.has_extension(x, sbml_extentions) ]

  sbml_files = [ x for x in arguments.filenames
                   if utils.has_extension(x, sbml_extentions) ]

  dictionary = dict()
  for param_file in param_files:
    parameters.parse_param_file(param_file, dictionary=dictionary)
  for sbml_file in sbml_files:
    parameterise_model_file(sbml_file, dictionary)

if __name__ == "__main__":
  run()
