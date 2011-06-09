"""A simple script to give a simple outline of an sbml model"""
import sys
import xml.dom.minidom

class Reaction:
  def __init__(self):
    """Initialise a reaction which has no reactants or products
       these may in turn be added with 'add_reactant' and 'add_product'"""
    self.reactants = []
    self.products = []
  
  def add_reactant(self, r):
    """Add a reactant into the reaction definition"""
    self.reactants.append(r)

  def add_product(self, p):
    """Add a product to the reaction"""
    self.products.append(p)

  def format_reaction(self):
    results = ""

    # So the first reactant has nothing attached to the front of it
    prefix = ""
    for reactant in self.reactants:
      results += prefix
      results += reactant
      # Hence, every reactant other than the first will have ", "
      # prefixed to the front of it, separating it from the previous one
      prefix = ", "
    results += " --> "
    prefix = ""
    for product in self.products:
      results += prefix
      results += product
      prefix = ", "
  
    return results


def obfuscate_named_attribute(element, attribute):
  """A helper function to obfuscate a given named attribute if
     that attribute is present in the given element"""
  name = element.getAttribute(attribute)
  if attribute and name:
    new_name = symbol_factory.get_new_symbol(name)
    element.setAttribute(attribute, new_name)

def obfuscate_compartment_type(compartment_type):
  """Simple as is function to obfuscate a compartment type definition"""
  obfuscate_named_attribute(compartment_type, "name")
  obfuscate_named_attribute(compartment_type, "id")

def obfuscate_compartment(compartment):
  """Simple obfuscate the names in a compartment definition"""
  obfuscate_named_attribute(compartment, "name")
  obfuscate_named_attribute(compartment, "id")
  
def obfuscate_species(species):
  """Simple obfuscation for a species declaration"""
  obfuscate_named_attribute(species, "name")
  obfuscate_named_attribute(species, "id")
  obfuscate_named_attribute(species, "compartment")

def obfuscate_parameter(parameter):
  """Simple obfuscation of a parameter definition"""
  obfuscate_named_attribute(parameter, "name")
  obfuscate_named_attribute(parameter, "id")

def obfuscate_maths(maths):
  """A function to obfuscate math nodes. Essentially all we're doing is
     descending down through the sub-elements until we come across
     a 'ci' node at which point we symbol-obfuscate its contents"""
  if maths.nodeType == maths.ELEMENT_NODE:
    if maths.tagName == "ci":
      text = maths.firstChild
      old_string = text.data
      new_string = symbol_factory.get_new_symbol(old_string)
      text.data = new_string
    else:
      for child in maths.childNodes:
        obfuscate_maths(child)

def obfuscate_raterule(raterule):
  """A simple as is possible to obfuscate a rate rule"""
  obfuscate_named_attribute(raterule, "variable")
  maths = raterule.getElementsByTagName("math")[0]
  obfuscate_maths(maths)

def obfuscate_init_assign(init_assign):
  """A simple function to obfuscate an initial assignment"""
  obfuscate_named_attribute(init_assign, "symbol")
  maths = init_assign.getElementsByTagName("math")[0]
  obfuscate_maths(maths)

def obfuscate_species_reference(spec_ref):
  """A simple function to obfuscate a <speciesReference> tag"""
  obfuscate_named_attribute(spec_ref, "species")

def obfuscate_kinetic_law(kin_law):
  """A simple function to obfuscate a <kineticLaw> tag"""
  # There should only be one 'maths' element but it doesn't hurt
  # to call it as a possible many maths sub-elements
  obfuscate_list_of(kin_law, "math", obfuscate_maths)

def obfuscate_reaction(reaction):
  """A function to obfuscate a reaction"""
  obfuscate_named_attribute(reaction, "id")
  obfuscate_list_of_list_of (reaction,
                             "listOfReactants",
                             "speciesReference",
                             obfuscate_species_reference)
  obfuscate_list_of_list_of (reaction,
                             "listOfProducts",
                             "speciesReference",
                             obfuscate_species_reference)
  obfuscate_list_of(reaction, "kineticLaw", obfuscate_kinetic_law)

def name_of_species_reference(spec_ref):
  name = spec_ref.getAttribute("species")
  return name

def get_reaction_of_element(reaction_element):
  """A function to return a reaction object from a reaction sbml element"""
  reaction = Reaction()
  reactants = get_elements_from_lists_of_list("listOfReactants",
                                              "speciesReference",
                                              name_of_species_reference,
                                              reaction_element)
  for reactant in reactants:
    reaction.add_reactant(reactant)

  # Of course we do the same for products
  products = get_elements_from_lists_of_list("listOfProducts",
                                             "speciesReference",
                                             name_of_species_reference,
                                             reaction_element)
  for product in products:
    reaction.add_product(product)

  return reaction
   

def get_elements_from_lists_of_list(list_element_name,
                                    element_name,
                                    extract_fun,
                                    parent_node):
  lists_of_element = parent_node.getElementsByTagName(list_element_name)
  result_objects = []
  for list_of_element in lists_of_element:
    elements = list_of_element.getElementsByTagName(element_name)
    for element in elements:
      result_obj = extract_fun(element)
      result_objects.append(result_obj)
  return result_objects

def get_list_of_reactions(model):
  return get_elements_from_lists_of_list("listOfReactions",
                                         "reaction",
                                         get_reaction_of_element,
                                         model)

def print_amount(num, singular, plural):
  if num == 0:
    print("No " + plural)
  elif num == 1:
    print("1 " + singular)
  else:
    print(str(num) + " " + plural)


class ExprVistor: 
  def __init__(self):
    self.result = ""

  def get_results(self):
    return self.result

  def print_str(self, string):
    self.result += string

  def generic_visit(self, element):
    if element.nodeType == element.ELEMENT_NODE :
      tag_name = element.tagName
      if tag_name == "apply":
        self.visit_Apply(element)
      elif tag_name == "ci":
        self.visit_ci(element)
      else:
        self.print_str("unknown-tag: " + tag_name)
    else: 
      return ""

  def visit_ci(self, element):
    self.print_str(element.firstChild.data)

  def visit_Apply(self, element):
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

  def visit_Maths(self, maths):
    for child in maths.childNodes:
      if maths.nodeType == maths.ELEMENT_NODE:
        self.generic_visit(child)

def format_rate_rule(raterule):
  maths_visitor = ExprVistor()
  maths = raterule.getElementsByTagName("math")[0]
  maths_visitor.visit_Maths(maths)
  name = raterule.getAttribute("variable")
  return name + " = " + maths_visitor.get_results()


def outline_rate_rules(model):
  rate_rules = get_elements_from_lists_of_list("listOfRules",
                                               "rateRule",
                                               format_rate_rule,
                                               model)
  no_rate_rules = len(rate_rules)
  print_amount(len(rate_rules), "rate rule", "rate rules")
  for rate_rule in rate_rules:
    print("  " + rate_rule)

def outline_model(model):
  """Format and print out an outline for the given model"""
  reactions = get_list_of_reactions(model)
  print_amount(len(reactions), "reaction", "reactions")
  for reaction in reactions:
    print("  " + reaction.format_reaction())

  outline_rate_rules(model)

def outline_sbml_file(filename):
  """Parse in a file as an SBML model, extract the outline information
     and then format that information and print it out"""
  dom = xml.dom.minidom.parse(filename)
  outline_model(dom.getElementsByTagName("model")[0])
  


def run():
  """Perform the banalities of command-line argument processing and
     then get under way with the proper grunt work of obfuscating the model"""
  # The command line arguments not including this script itself
  arguments    = sys.argv 
  # file names are arguments that don't affect the configuration such
  # as limit=10 or x=k1
  filenames = [ x for x in arguments if '=' not in x and
                  not x.endswith(".py") ]
  
  for filename in filenames:
    outline_sbml_file(filename)


if __name__ == "__main__":
  run()
