""" A simple script to replace the values of parameters within an
    sbml model file"""
import xml.dom.minidom
import os
import argparse
import parameters


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
  

def has_extension(filename, extensions):
  """Returns true if the given filename has one of the given
     list of extensions"""
  file_extension = os.path.splitext(filename)[1]
  return file_extension in extensions

def run():
  """Perform the banalities of command-line argument processing and
     and then get under way in parameterising the model"""
  description = "Parameterise an SBML model based on a given param file"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="input files, parameters and sbml model files")
  # parser.add_argument('--pretty', action='store_true',
  #                     help="Pretty print the xml")

  arguments = parser.parse_args()

  sbml_extentions = [ ".xml", ".sbml" ]
  param_files = [ x for x in arguments.filenames
                        if not has_extension(x, sbml_extentions) ]

  sbml_files = [ x for x in arguments.filenames
                   if has_extension(x, sbml_extentions) ]

  dictionary = dict()
  for param_file in param_files:
    parameters.parse_param_file(param_file, dictionary=dictionary)
  for sbml_file in sbml_files:
    parameterise_model_file(sbml_file, dictionary)

if __name__ == "__main__":
  run()
