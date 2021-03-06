"""a simple script to give a simple outline of an sbml model"""
import xml.dom.minidom
import argparse

import sbml_ast
import expressions


def get_element_children(node):
  """A utility function for dealing with xml. Often you wish to deal
     with the children of an element, but only those that are actually
     nodes, not the text children. For example
     <math>
       <cn>0.1</cn>
     </math>
     The math element has three children, the node child <cn> plus the text
     childrend representing the whitespace either side of the <cn> child.
     So this method picks out only those children which are nodes.
  """
  if node.hasChildNodes():
    return [ x for x in node.childNodes
                 if x.nodeType == x.ELEMENT_NODE ]
  else:
    return []

def name_of_species_reference(spec_ref):
  """Return the name of the species referred to within a
     speciesReference sbml element. This is also used for 
     the modifierSpeciesReference elements
  """
  name = spec_ref.getAttribute("species")
  return name

def react_partic_of_species_ref(spec_ref):
  """Return the 'ReactionParticpant' of a speciesReference sbml
     element. This is also used for 
     the modifierSpeciesReference elements
  """
  name = spec_ref.getAttribute("species")
  stoich = spec_ref.getAttribute("stoichiometry")
  try:
    stoich = int(stoich)
  except ValueError:
    stoich = None
  if stoich == None:
    return sbml_ast.ReactionParticipant(name)
  else :
    return sbml_ast.ReactionParticipant(name, stoich=int(stoich))
  return name

def get_reaction_of_element(reaction_element):
  """a function to return a reaction object from a reaction sbml element"""
  name = reaction_element.getAttribute("id")
  reaction = sbml_ast.Reaction(name)
  location = reaction_element.getAttribute("compartment")
  if location:
    reaction.location = location
  reactants = get_elements_from_lists_of_list("listOfReactants",
                                              "speciesReference",
                                              react_partic_of_species_ref,
                                              reaction_element)
  for reactant in reactants:
    reaction.add_reactant(reactant)

  # of course we do the same for products
  products = get_elements_from_lists_of_list("listOfProducts",
                                             "speciesReference",
                                             react_partic_of_species_ref,
                                             reaction_element)
  for product in products:
    reaction.add_product(product)


  # And also for the modifiers
  modifiers = get_elements_from_lists_of_list("listOfModifiers",
                                             "modifierSpeciesReference",
                                             react_partic_of_species_ref,
                                             reaction_element)
  for modifier in modifiers:
    reaction.add_modifier(modifier)

  # This whole thing is not very defensively programmed
  # in that we just keep taking the first of any appropriate
  # sub nodes.
  kinetic_laws = reaction_element.getElementsByTagName("kineticLaw")
  if kinetic_laws:
    kinetic_law = kinetic_laws[0]
    math_elements = kinetic_law.getElementsByTagName("math")
    if math_elements:
      math_element = math_elements[0]
      reaction.kinetic_law = math_element

  return reaction
   

def get_fundef_of_element(fundef_element):
  """Return a FunctionDefinition object from a functionDefinition
     element
  """
  fundef = sbml_ast.FunctionDefinition()
  name = fundef_element.getAttribute("name")
  fundef.name = name

  lambda_elements = fundef_element.getElementsByTagName("lambda")
  if lambda_elements:
    lambda_element = lambda_elements[0]

    bvar_elements = []
    other_children = []
    for child_node in lambda_element.childNodes:
      if child_node.nodeType == child_node.ELEMENT_NODE:
        # is this case-sensitive?
        if child_node.tagName == "bvar":
          bvar_elements.append(child_node)
        else:
          other_children.append(child_node)

    for bvar_element in bvar_elements:
      ci_elements = bvar_element.getElementsByTagName("ci")
      if ci_elements:
        name = ci_elements[0].firstChild.data.lstrip().rstrip() 
        fundef.parameters.append(name)
    if other_children:
      fundef.body = parse_expression(other_children[0])

  return fundef

  

def get_elements_from_lists_of_list(list_element_name,
                                    element_name,
                                    extract_fun,
                                    parent_node):
  """returns a list of objects from an element which contains a list
     of elements. Such as 'listOfReactions'. This takes in the name
     of the list element, the name of the elements which represent the
     items of the list and a function which creates an object for each
     item in the list from the item element in the list"""
  lists_of_element = parent_node.getElementsByTagName(list_element_name)
  result_objects = []
  for list_of_element in lists_of_element:
    elements = list_of_element.getElementsByTagName(element_name)
    for element in elements:
      result_obj = extract_fun(element)
      result_objects.append(result_obj)
  return result_objects


def get_list_of_fundefs(model):
  """Returns a list or FunctionDefinition objects from an sbml model"""
  fundefs = get_elements_from_lists_of_list("listOfFunctionDefinitions",
                                            "functionDefinition",
                                            get_fundef_of_element,
                                            model)
  return fundefs 

def get_list_of_reactions(model, ignore_sources=False, ignore_sinks=False):
  """Returns a list of reaction objects from an sbml model"""
  reactions = get_elements_from_lists_of_list("listOfReactions",
                                              "reaction",
                                              get_reaction_of_element,
                                              model)
  if ignore_sources:
    reactions = [ r for r in reactions if not r.is_source() ] 
  if ignore_sinks:
    reactions = [ r for r in reactions if not r.is_sink() ]
 
  return reactions


 

def get_species_of_element(species_element):
  """Return a species object from a species element, in other words
     parse a species element"""
  return sbml_ast.Species(species_element)

def get_list_of_species(model):
  """Return the list of species from an sbml model"""
  return get_elements_from_lists_of_list("listOfSpecies",
                                         "species",
                                         get_species_of_element,
                                         model)

def get_parameter_of_element(param_element):
  """Return a Parameter object from a parameter element,
     in other words parse a parameter element
  """
  return sbml_ast.Parameter(param_element)

def get_list_of_parameters(model):
  """Return the list of parameters from an sbml model"""
  return get_elements_from_lists_of_list("listOfParameters",
                                         "parameter",
                                         sbml_ast.Parameter,
                                         model)
def get_assignment_rule_of_element(assign_rule_element):
  """Returns the AssignmentRule representation of an 
     AssignmentRule sbml element"""
  variable_name = assign_rule_element.getAttribute("variable")
  math_elements = assign_rule_element.getElementsByTagName("math")
  expression = None
  if math_elements:
    expression = math_elements[0]

  return sbml_ast.AssignmentRule(variable_name, expression)
 
def get_list_of_assignment_rules(model):
  """Return the list of assignment rules from an sbml model"""
  return get_elements_from_lists_of_list("listOfRules",
                                         "assignmentRule",
                                         get_assignment_rule_of_element,
                                         model)

def get_init_assign_of_element(init_assign_element):
  """Returns the list InitialAssignment representation of an
     InitialAssignment rule sbml element
  """
  symbol_name = init_assign_element.getAttribute("symbol")
  math_elements = init_assign_element.getElementsByTagName("math")
  expression = None
  if math_elements:
    expression = parse_expression(math_elements[0])

  return sbml_ast.InitialAssignment(symbol_name, expression)
  

def get_list_of_init_assigns(model):
  """Return the list of initial assignments from an sbml model"""
  return get_elements_from_lists_of_list("listOfInitialAssignments",
                                         "initialAssignment",
                                         get_init_assign_of_element,
                                         model)

def print_amount(num, singular, plural):
  """Utility function to print an expression for the number of
     something. Picks the singular or plural noun appropriately"""
  if num == 0:
    print("no " + plural)
  elif num == 1:
    print("1 " + singular)
  else:
    print(str(num) + " " + plural)


class ExprVisitor(object):
  """A parent class for classes which descend through SBML
     expression elements storing a result along the way.
  """
  def __init__(self):
    self.result = None

  def get_results(self):
    """Return the string result of visiting the expression"""
    return self.result

  def generic_visit(self, element):
    """The main entry point, decides on the kind of element we have
       and calls the appropriate visitor function from below"""
    if element.nodeType == element.ELEMENT_NODE :
      tag_name = element.tagName
      if tag_name == "apply":
        self.visit_apply(element)
      elif tag_name == "ci":
        self.visit_ci(element)
      elif tag_name == "cn":
        self.visit_cn(element)
      elif tag_name == "math":
        self.visit_maths(element)
      else:
        raise ValueError("unknown-tag: " + tag_name)
    else: 
      return ""

  def visit_maths(self, maths):
    """Visit a 'maths' element"""
    element_children = get_element_children(maths)
    no_chilren = len(element_children)
    if no_chilren == 1: 
      self.generic_visit(element_children[0])
    else:
      raise ValueError("math element with: " +
                        str(no_chilren) + " children")

  ###################################
  # These are the unimplemented methods that you would be likely
  # to override for your expression visitor.
  def visit_ci(self, element):
    """Visit a 'ci' element"""
    raise NotImplementedError("visit_ci element for expression visitor")

  def visit_cn(self, element):
    """Visit a 'cn' element"""
    raise NotImplementedError("visit_cn element for expression visitor")

  def visit_apply(self, element):
    """Visit an 'apply' element"""
    raise NotImplementedError("visit_apply element for expression visitor")


 
def parse_expression(expr_element):
  """Parse an sbml expression element into an sbml_ast.Expression
     subclass
  """
  expr_builder = ExprBuilder()
  expr_builder.generic_visit(expr_element)
  return expr_builder.get_results()

class ExprBuilder(ExprVisitor):
  """An expression visitor class which will build an abstract syntax
     representation of the expression defined by an sbml expression
     element (typically under a math element).
  """
  def __init__(self):
    super(ExprBuilder, self).__init__()
    self.result = None

  def visit_ci(self, element):
    """Visit a 'ci' element"""
    name = element.firstChild.data.lstrip().rstrip()
    self.result = expressions.NameExpression(name)

  def visit_cn(self, element):
    """Visit a 'cn' element"""
    number_string =  element.firstChild.data.lstrip().rstrip()
    self.result = expressions.NumExpression(float(number_string))

  def visit_apply(self, element):
    """Visit an 'apply' element"""
    children = [ x for x in element.childNodes 
                   if x.nodeType == x.ELEMENT_NODE
               ]
    function = children[0]
    function_name = function.tagName
 
    arguments = []
    for child in children[1:]:
      self.generic_visit(child)
      arguments.append(self.result)
    apply_expr = expressions.ApplyExpression(function_name, arguments) 
    self.result = apply_expr

class ExprFormatter(ExprVisitor):
  """A class which descends through SBML math expressions storing
     a formatted string representing the math expression"""
  def __init__(self):
    super(ExprFormatter, self).__init__()
    self.result = ""

  def print_str(self, string):
    """A (private) utility function for printing to the result string"""
    self.result += string

  def visit_ci(self, element):
    """Visit a 'ci' element"""
    self.print_str(element.firstChild.data)

  def visit_cn(self, element):
    """Visit a 'cn' element"""
    self.print_str(element.firstChild.data)

  def visit_apply(self, element):
    """Visit an 'apply' element"""
    
    children = [ x for x in element.childNodes 
                   if x.nodeType == x.ELEMENT_NODE
               ]
    function = children[0]
    function_name = function.tagName
    arg_strings = []
    for child in children[1:]:
      current_result = self.result
      self.generic_visit(child)
      arg_strings.append(self.result)
      self.result = current_result

    formatted = expressions.show_apply_expression(function_name, arg_strings)
    self.print_str(formatted)

def evaluate_function_application(function_name, arguments):
  """Takes in the function name and the arguments which we will assume
     are numbers and returns the result of the application.
  """
  if function_name == "plus":
    return sum(arguments)
  elif function_name == "minus":
    answer = arguments[0]
    for arg in arguments[1:]:
      answer -= arg
    return answer
  elif function_name == "times":
    answer = arguments[0]
    for arg in arguments[1:]:
      answer *= arg
    return answer
  elif function_name == "div":
    answer = arguments[0]
    for arg in arguments[1:]:
      answer /= arg
    return answer
  elif function_name == "power":
    # power is interesting because it associates to the right
    exponent = 1
    # counts downwards from the last index to the 0.
    # As an example, consider power(3,2,3), the answer should be
    # 3 ** (2 ** 3) = 3 ** 8 = 6561, not (3 ** 2) ** 3 = 9 ** 3 = 81
    # going through our loop here we have
    # exp = 1
    # exp = 3 ** exp = 3
    # exp = 2 ** exp = 2 ** 3 = 8
    # exp = 3 ** exp = 3 ** 8 = 6561
    for i in range(len(arguments) - 1, -1, -1):
      exponent = arguments[i] ** exponent
    return exponent 
  else:
    raise ValueError("Unknown function: " + function_name)

 
      

class ExprEvaluator(ExprVisitor):
  """A class which descends through SBML math expressions to evaluate
     the current value of the expression. It relies on a dictionary
     mapping names to current values to evaluate symbolic names, this
     should be set prior to visiting the expression element with
     expr_eval.name_mapping = ...
  """
  def __init__(self):
    super(ExprEvaluator, self).__init__()
    self.result = None
    self.name_mapping = None


  def visit_ci(self, element):
    """Visit a 'ci' element"""
    name = element.firstChild.data.lstrip().rstrip() 
    self.result = self.name_mapping[name]

  def visit_cn(self, element):
    """Visit a 'cn' element"""
    self.result = float(element.firstChild.data)

  def visit_apply(self, element):
    """Visit an 'apply' element"""
    children = [ x for x in element.childNodes 
                   if x.nodeType == x.ELEMENT_NODE
               ]
    function = children[0]
    function_name = function.tagName

    # Evaluate all the arguments to the function
    def evaluate_argument(arg):
      """Quick helper function to evaluate an argument"""
      self.generic_visit(arg)
      return self.result

    arguments = [ evaluate_argument(a) for a in children[1:] ]
    self.result = evaluate_function_application(function_name, arguments)
 
def evaluate_expression(population_dictionary, expression):
  """ Expressions are evaluated via an ExprEvaluator. This is just a
      utility function which creates the evaluator, initialises it, calls
      the visit function, obtains the result and returns it.
  """
  expr_evaluator = ExprEvaluator()
  expr_evaluator.name_mapping = population_dictionary
  expr_evaluator.generic_visit(expression)
  value = expr_evaluator.get_results()
  return value 

def format_expression(expr):
  """ A utility function which creates an expression formattor and then
      applies that to an expression returning the results.
  """
  expr_visitor = ExprFormatter()
  expr_visitor.generic_visit(expr)
  return expr_visitor.get_results()
 
def format_math_element(maths):
  """Format an math element as an expression
  """
  expr_visitor = ExprFormatter()
  expr_visitor.visit_maths(maths)
  return expr_visitor.get_results()
 

def format_rate_rule(raterule):
  """Return a string representing a 'rateRule' element"""
  name = raterule.getAttribute("variable")
  maths = raterule.getElementsByTagName ("math")[0]
  return name + " = " + format_math_element(maths)


def outline_rate_rules(model):
  """Print out an outline for the rate rules of an SBML model"""
  rate_rules = get_elements_from_lists_of_list("listOfRules",
                                               "rateRule",
                                               format_rate_rule,
                                               model)
  no_rate_rules = len(rate_rules)
  print_amount(no_rate_rules, "rate rule", "rate rules")
  for rate_rule in rate_rules:
    print("  " + rate_rule)

def outline_model(model, arguments):
  """format and print out an outline for the given model"""
  parameters = get_list_of_parameters(model)
  
  print_amount(len(parameters), "parameter", "parameters")
  for param in parameters:
    print ("  " + param.format_param())


  fundefs = get_list_of_fundefs(model)

  print_amount(len(fundefs), "function definition", "function definitions")
  for fundef in fundefs:
    print (fundef.format_fundef())

  reactions = get_list_of_reactions(model,
                  ignore_sources = arguments.ignore_sources,
                  ignore_sinks = arguments.ignore_sinks)


  print_amount(len(reactions), "reaction", "reactions")
  for reaction in reactions:
    print("  " + reaction.format_reaction())

  species_list = get_list_of_species(model)
  # Ha okay we don't really need 'print_amount' since the plural
  # and singular of species is the same.
  print_amount(len(species_list), "species", "species")
  for species in species_list:
    print ("  " + species.format_species())
    for reaction in reactions:
      if reaction.involves(species):
        print ("       " + reaction.format_reaction())

  outline_rate_rules(model)

def get_model_from_sbml_file(filename):
  """Parses the sbml file as an xml file and then returns the appropriate
     model element
  """
  dom = xml.dom.minidom.parse(filename)
  model = dom.getElementsByTagName("model")[0]
  return model
 
def outline_sbml_file(filename, arguments):
  """parse in a file as an sbml model, extract the outline information
     and then format that information and print it out"""
  model = get_model_from_sbml_file(filename)
  outline_model(model, arguments)
  
def run():
  """perform the banalities of command-line argument processing and
     then go ahead and calculate the outline for each model file"""
  description = "Print out an outline of an SBML file"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="an sbml file to outline")
  parser.add_argument("--ignore-sources",
                      action="store_true", default=False,
    help="Ingore source reactions when outlining the model")
  parser.add_argument("--ignore-sinks",
                      action="store_true", default=False,
    help="Ignore sink reactions when outlining the model")
  arguments = parser.parse_args()

  for filename in arguments.filenames:
    outline_sbml_file(filename, arguments)

if __name__ == "__main__":
  run()
