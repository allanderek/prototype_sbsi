"""a simple script to give a simple outline of an sbml model"""
import xml.dom.minidom
import argparse

import sbml_ast

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
  if not stoich:
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
    expression = math_elements[0]

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


class ExprVisitor: 
  """A class which descends through SBML math expressions storing
     a formatted string representing the math expression"""
  def __init__(self):
    self.result = ""

  def get_results(self):
    """Return the string result of visiting the expression"""
    return self.result

  def print_str(self, string):
    """A (private) utility function for printing to the result string"""
    self.result += string

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
      else:
        self.print_str("unknown-tag: " + tag_name)
    else: 
      return ""

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
    function_dict = { "plus" : "+", 
                      "minus" : "-",
                      "divide" : "/",
                      "times" : "*",
                      "power" : "^",
                     }
    if function_name in function_dict :
      self.print_str("(")
      self.generic_visit(children[1])
      self.print_str(" " + function_dict[function_name] + " ")
      self.generic_visit(children[2])
      self.print_str(")") 
    else:
      self.print_str (function_name + "(")
      for child in children[1:]:
        self.generic_visit(child)
        self.print_str (", ")
      self.print_str (")")

  def visit_maths(self, maths):
    """Visit a 'maths' element"""
    for child in maths.childNodes:
      if maths.nodeType == maths.ELEMENT_NODE:
        self.generic_visit(child)

def format_math_element(maths):
  """Format an math element as an expression
  """
  expr_visitor = ExprVisitor()
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

def outline_sbml_file(filename, arguments):
  """parse in a file as an sbml model, extract the outline information
     and then format that information and print it out"""
  dom = xml.dom.minidom.parse(filename)
  model = dom.getElementsByTagName("model")[0]
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
