"""A simple script to guess the initparams file for an SBML model
   to be used in an sbsi optimisation"""
import sys
import os
import xml.dom.minidom
import outline_sbml
import argparse

def change_filename_ext(filename, new_ext):
  """Returns a new file name based on the first but with the extension
     replaced by the given extension. The new_ext should include the
     '.' separator if that is desired"""
  basename = os.path.splitext(filename)[0]
  return basename + new_ext
  

class Parameter:
  def __init__(self, identifier, name, value):
    self.identifier = identifier
    self.name = name
    self.value = value 

  def format_init(self, factor):
    """Format the parameter specifically for an sbsi initparams file.
       Note that this entails calculating the high and low values"""
    high = self.value * factor
    low  = self.value / factor
    line = (self.name + "\t" + str(low) + "\t" +
                               str(high) + "\t" +
                               str(self.value))
    return line        

def get_parameter_of_element(param_element):
  """a function to return a parameter object from a
     parameter sbml element"""
  ident = param_element.getAttribute("id")
  name = param_element.getAttribute("name")
  value = float(param_element.getAttribute("value") )
  parameter = Parameter(ident, name, value)
  
  return parameter

def get_parameters_from_model(model):
  params = outline_sbml.get_elements_from_lists_of_list("listOfParameters",
                "parameter",
                get_parameter_of_element,
                model)
  return params

class AssignmentRule:
  def __init__(self, variable, math):
    self.variable = variable
    self.math = math

def get_assignment_rule_from_element(arule_element):
  variable = arule_element.getAttribute("variable")
  math = arule_element.getElementsByTagName("math")
  return AssignmentRule(variable, math)

def get_assignment_rules_from_model(model):
  arules = outline_sbml.get_elements_from_lists_of_list("listOfRules",
               "assignmentRule",
               get_assignment_rule_from_element,
               model)
  return arules

def init_params_model_file(filename, factor):
  """Parse in a file as an SBML model, and extract the parameter
     information from it"""
  dom = xml.dom.minidom.parse(filename)
  model = dom.getElementsByTagName("model")[0]
  params = get_parameters_from_model(model)
  arules = get_assignment_rules_from_model(model)
  arule_names = [ arule.variable for arule in arules ]

  output_filename = change_filename_ext(filename, ".initparams")
  output_file = open (output_filename, "w") 
  for param in params:
    if param.name not in arule_names: 
      output_file.write (param.format_init(factor) + "\n")
  output_file.close()
 
def run():
  """Perform the banalities of command-line argument processing
     and try to guess the initparams for the given SBML model"""
  description = "Analyse SBML files for invariants"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="an sbml file to check invariants for")
 
  arguments = parser.parse_args()

  # factor should be a command-line option
  factor = 10.0

  for filename in arguments.filenames:
    init_params_model_file(filename, factor)

if __name__ == "__main__":
  run()