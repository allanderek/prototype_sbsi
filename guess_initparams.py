"""A simple script to guess the initparams file for an SBML model
   to be used in an sbsi optimisation"""
import xml.dom.minidom
import outline_sbml
import argparse
import utils
 

class Parameter:
  """A simple class to represent a parameter description"""
  def __init__(self, identifier):
    self.identifier = identifier
    self.name = identifier
    self.value = None 
    self.high = None
    self.low = None
    self.initial_step = None
    self.use = True

  def format_init(self, factor):
    """Format the parameter specifically for an sbsi initparams file.
       Note that this entails calculating the high and low values"""
    if not self.high:
      high = self.value * factor
    else:
      high = self.high
    if not self.low:
      low  = self.value / factor
    else:
      low = self.low
    if not self.initial_step:
      initial_step = (high - low) / 100
    else:
      initial_step = self.initial_step
    line = (self.name + "\t" + str(low) + "\t" +
                               str(high) + "\t" +
                               str(initial_step))
    return line        

def get_parameter_of_element(param_element):
  """a function to return a parameter object from a
     parameter sbml element"""
  ident = param_element.getAttribute("id")
  name = param_element.getAttribute("name")
  value_str = param_element.getAttribute("value") 
  if value_str:
    value = float(value_str)
    parameter = Parameter(ident)
    parameter.value = value 
    if name:
      parameter.name = name
    return parameter
  return None

# TODO: this and 'get_parameter_of_element' are somewhat similarly defined
# in outline_sbml, merge the two definitions and have only one.
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


def float_from_attribute(element, attribute):
  """Returns the float value from a named attribute of the given
     element. None if the attribute doesn't exist or is empty.
  """
  value_str = element.getAttribute(attribute)
  if not value_str:
    return None
  return float(value_str)


def get_param_of_history(history_element):
  """Returns a parameter description from a parameterHistory element"""
  identifier = history_element.getAttribute("id")
  parameter = Parameter(identifier)
  name = history_element.getAttribute("name")
  if name:
    parameter.name = name

  parameter.low = float_from_attribute(history_element, "min")
  parameter.high = float_from_attribute(history_element, "max")
  parameter.initial_step = float_from_attribute(history_element, 
                                                "initialstep")
  use = history_element.getAttribute("use")
  if use == "No":
    parameter.use = False
  else:
    parameter.use = True

  return parameter


def params_from_param_histories(model):
  """Get all the parameter histories in the model file and parse
     these into parameter descriptions suitable for creating an
     init params file
  """
  annotations = model.getElementsByTagName("annotation")
  params = []
  for annotation in annotations:
    list_name = "listOfParameterHistories"
    annot_params = outline_sbml.get_elements_from_lists_of_list(list_name,
                      "parameterHistory",
                      get_param_of_history,
                      annotation)
    params += [ p for p in annot_params if p.use ]
                      
   
  return params 

def init_params_model_file(filename, arguments):
  """Parse in a file as an SBML model, and extract the probably
     optimisable parameters from it. Essentially then that is all
     the parameters which are not associated with an assignment rule"""
  dom = xml.dom.minidom.parse(filename)
  model = dom.getElementsByTagName("model")[0]

  # So first we attempt to get the parameters from parameter histories
  # assuming that hasn't been disabled by the user
  if not arguments.ignore_param_histories:
    params = params_from_param_histories(model)
  # If we got any parameters from the parameter histories then go
  # ahead and return those, otherwise we continue and just get them
  # from the parameter definitions
  if params:
    return params
 
  params = get_parameters_from_model(model)

  # We ignore parameters that are specifically assigned to via an
  # initial assignment since it is likely that these will not be
  # optimised for, in particular they are likely to be some combination
  # of parameters which are to be optimised for.
  arules = get_assignment_rules_from_model(model)
  arule_names = [ arule.variable for arule in arules ]
  return [ p for p in params if p.name not in arule_names ]



def create_init_params_file(filename, arguments):
  """Get the parameters of an sbml file which we think are likely to
     be 'optimisable' and based on that create an initparams file"""
  output_filename = utils.change_filename_ext(filename, ".initparams")
  output_file = open (output_filename, "w") 
  params = init_params_model_file(filename, arguments)
  for param in params:
    output_file.write (param.format_init(arguments.factor) + "\n")
  output_file.close()

def get_argument_parser():
  """Return the arguments parser for guess_initparams"""
  description = "Analyse SBML files for invariants"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="an sbml file to guess init params for")
  parser.add_argument('--factor', action='store',
                      type=float, default=10.0,
                      help="Set the factor either side for min-max")
  parser.add_argument('--ignore-param-histories', action='store_true',
                      help="Ignore any parameter histories in the SBML")
  return parser


def run():
  """Perform the banalities of command-line argument processing
     and try to guess the initparams for the given SBML model"""
  parser = get_argument_parser()
  arguments = parser.parse_args()

  for filename in arguments.filenames:
    create_init_params_file(filename, arguments)

if __name__ == "__main__":
  run()
