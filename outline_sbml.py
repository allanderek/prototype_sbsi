"""a simple script to give a simple outline of an sbml model"""
import xml.dom.minidom
import argparse

class ReactionParticipant:
  """A simple class to represent a reaction participant"""
  def __init__(self, name, stoich=1):
    self.name = name
    self.stoich = stoich

  def get_name(self):
    """Return the name of the reaction participant"""
    return self.name
  def get_stoichiometry(self):
    """Return the stoichiometry of this reaction participant
       within the reaction"""
    return self.stoich

class Species:
  """A simple class to represent a species within an sbml model"""
  def __init__(self, name, compartment):
    self.name = name
    self.compartment = compartment
  def get_name(self):
    """ return the name of the species"""
    return self.name

class Reaction:
  """A class which represents a reaction"""
  def __init__(self, name):
    """initialise a reaction which has no reactants or products
       these may in turn be added with 'add_reactant' and 'add_product'"""
    self.reactants = []
    self.products = []
    self.name = name

  def get_name(self):
    """return the name of the reaction"""
    return self.name
  
  def add_reactant(self, reactant):
    """add a reactant into the reaction definition"""
    self.reactants.append(reactant)

  def add_product(self, product):
    """add a product to the reaction"""
    self.products.append(product)

  def get_reactants(self):
    """return the list of reactant names"""
    return self.reactants

  def get_products(self):
    """return the list of product names"""
    return self.products

  def is_sink(self):
    """returns true if the reaction is a sink,
       in that it has at least one reactant and no products"""
    return self.reactants and not self.products

  def is_source(self):
    """returns true if the reaction is a source reaction,
       in that it has at least one product and no reactants"""
    return self.products and not self.reactants

  def format_reaction(self):
    """Return a string representing the reaction in a format 
       suitable for human consumption"""
    results = self.name + ": "

    # so the first reactant has nothing attached to the front of it
    prefix = ""
    for reactant in self.reactants:
      results += prefix
      results += reactant.get_name()
      # hence, every reactant other than the first will have ", "
      # prefixed to the front of it, separating it from the previous one
      prefix = ", "
    results += " --> "
    prefix = ""
    for product in self.products:
      results += prefix
      results += product.get_name()
      prefix = ", "
  
    return results


def name_of_species_reference(spec_ref):
  """Return the name of the species referred to within a
     speciesReference sbml element"""
  name = spec_ref.getAttribute("species")
  return name

def get_reaction_of_element(reaction_element):
  """a function to return a reaction object from a reaction sbml element"""
  name = reaction_element.getAttribute("id")
  reaction = Reaction(name)
  reactants = get_elements_from_lists_of_list("listOfReactants",
                                              "speciesReference",
                                              name_of_species_reference,
                                              reaction_element)
  # todo we should add the modifiers

  for reactant in reactants:
    reaction.add_reactant(ReactionParticipant(reactant))

  # of course we do the same for products
  products = get_elements_from_lists_of_list("listOfProducts",
                                             "speciesReference",
                                             name_of_species_reference,
                                             reaction_element)
  for product in products:
    reaction.add_product(ReactionParticipant(product))

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

def get_list_of_reactions(model):
  """Returns a list of reaction objects from an sbml model"""
  return get_elements_from_lists_of_list("listOfReactions",
                                         "reaction",
                                         get_reaction_of_element,
                                         model)

def get_species_of_element(species_element):
  """Return a species object from a species element, in other words
     parse a species element"""
  name = species_element.getAttribute("id")
  compartment = species_element.getAttribute("compartment")
  species = Species(name, compartment)
  return species

def get_list_of_species(model):
  """Return the list of species from an sbml model"""
  return get_elements_from_lists_of_list("listOfSpecies",
                                         "species",
                                         get_species_of_element,
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
      else:
        self.print_str("unknown-tag: " + tag_name)
    else: 
      return ""

  def visit_ci(self, element):
    """Visit a 'ci' element"""
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

def format_rate_rule(raterule):
  """Return a string representing a 'rateRule' element"""
  maths_visitor = ExprVisitor()
  maths = raterule.getElementsByTagName ("math")[0]
  maths_visitor.visit_maths(maths)
  name = raterule.getAttribute("variable")
  return name + " = " + maths_visitor.get_results()


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
  reactions = get_list_of_reactions(model)
  print_amount(len(reactions), "reaction", "reactions")
  for reaction in reactions:
    if ( (arguments.ignore_sources and reaction.is_sources()) or
         (arguments.ignore_sinks and reaction.is_sink()) ):
      pass
    else:
      print("  " + reaction.format_reaction())

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
