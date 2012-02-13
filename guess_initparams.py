"""A simple script to guess the initparams file for an SBML model
   to be used in an sbsi optimisation"""
import xml.dom.minidom
import outline_sbml
import argparse
import utils
 

class Parameter:
  """A simple class to represent a parameter description"""
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
  value_str = param_element.getAttribute("value") 
  if value_str:
    value = float(value_str)
    parameter = Parameter(ident, name, value)
    return parameter
  return None

def get_parameters_from_model(model):
  """Get all the SBML parameters from the model"""
  params = outline_sbml.get_elements_from_lists_of_list("listOfParameters",
                "parameter",
                get_parameter_of_element,
                model)
  return [ p for p in params if p is not None ]

class AssignmentRule:
  """A simple class to represent an assignment rule"""
  def __init__(self, variable, math):
    self.variable = variable
    self.math = math

def get_assign_rule_from_element(arule_element):
  """Return an AssignmentRule instance from the xml element
     representing the assignment rule"""
  variable = arule_element.getAttribute("variable")
  math = arule_element.getElementsByTagName("math")
  return AssignmentRule(variable, math)

def get_assignment_rules_from_model(model):
  """Get all the assignment rules from the model"""
  arules = outline_sbml.get_elements_from_lists_of_list("listOfRules",
               "assignmentRule",
               get_assign_rule_from_element,
               model)
  return arules

def init_params_model_file(filename):
  """Parse in a file as an SBML model, and extract the probably
     optimisable parameters from it. Essentially then that is all
     the parameters which are not associated with an assignment rule"""
  dom = xml.dom.minidom.parse(filename)
  model = dom.getElementsByTagName("model")[0]
  params = get_parameters_from_model(model)
  arules = get_assignment_rules_from_model(model)
  arule_names = [ arule.variable for arule in arules ]

  return [ p for p in params if p.name not in arule_names ]


def create_init_params_file(filename, factor):
  """Get the parameters of an sbml file which we think are likely to
     be 'optimisable' and based on that create an initparams file"""
  params = init_params_model_file(filename)
  output_filename = utils.change_filename_ext(filename, ".initparams")
  output_file = open (output_filename, "w") 
  for param in params:
    output_file.write (param.format_init(factor) + "\n")
  output_file.close()
 
def run():
  """Perform the banalities of command-line argument processing
     and try to guess the initparams for the given SBML model"""
  description = "Analyse SBML files for invariants"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="an sbml file to guess init params for")
  parser.add_argument('--factor', action='store',
                      type=float, default=10.0,
                      help="Set the factor either side for min-max")

  arguments = parser.parse_args()

  for filename in arguments.filenames:
    create_init_params_file(filename, arguments.factor)

if __name__ == "__main__":
  run()
