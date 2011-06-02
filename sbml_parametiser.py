""" A simple script to replace the values of parameters within an
    sbml model file"""
import sys
import xml.dom.minidom
import os


def parameterise_model(model, dictionary):
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
  

def parse_param_file(param_filename, dictionary):
  param_file = open(param_filename, "r")

  for line in param_file:
    (name, value_string) = line.split(":", 1)
    value = float(value_string.lstrip().rstrip())
    # We should also check if the name already exists in the dictionary
    dictionary[name] = value

  param_file.close()

def has_extension(filename, extension):
  file_extension = os.path.splitext(filename)[1]
  return file_extension == extension

def run():
  """Perform the banalities of command-line argument processing and
     and then get under way in parameterising the model"""
  # The command line arguments not including this script itself
  arguments    = sys.argv 
  # file names are arguments that don't affect the configuration such
  # as limit=10 or x=k1
  filenames = [ x for x in arguments if '=' not in x and
                  not has_extension(x, ".py") ]
  
  param_files = [ x for x in filenames if not has_extension(x, ".xml") ] 
  sbml_files =  [ x for x in filenames if has_extension(x, ".xml") ]


  dictionary = dict()
  for param_file in param_files:
    parse_param_file(param_file, dictionary)
  for sbml_file in sbml_files:
    parameterise_model_file(sbml_file, dictionary)

if __name__ == "__main__":
  run()
