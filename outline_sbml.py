"""a simple script to give a simple outline of an sbml model"""
import sys
import xml.dom.minidom
import argparse

class ReactionParticipant:
  def __init__(self, name, stoich=1):
    self.name = name
    self.stoich = stoich

  def get_name(self):
    return self.name
  def get_stoichiometry(self):
    return self.stoich

class Species:
  def __init__(self, name, compartment):
    self.name = name
    self.compartment = compartment
  def get_name(self):
    return self.name

class Reaction:
  def __init__(self, name):
    """initialise a reaction which has no reactants or products
       these may in turn be added with 'add_reactant' and 'add_product'"""
    self.reactants = []
    self.products = []
    self.name = name

  def get_name(self):
    return self.name
  
  def add_reactant(self, r):
    """add a reactant into the reaction definition"""
    self.reactants.append(r)

  def add_product(self, p):
    """add a product to the reaction"""
    self.products.append(p)

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
  # todo we should add the 

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

def get_species_of_element(species_element):
  name = species_element.getAttribute("id")
  compartment = species_element.getAttribute("compartment")
  species = Species(name, compartment)
  return species

def get_list_of_species(model):
  return get_elements_from_lists_of_list("listOfSpecies",
                                         "species",
                                         get_species_of_element,
                                         model)

def print_amount(num, singular, plural):
  if num == 0:
    print("no " + plural)
  elif num == 1:
    print("1 " + singular)
  else:
    print(str(num) + " " + plural)


class ExprVisitor: 
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
        self.visit_apply(element)
      elif tag_name == "ci":
        self.visit_ci(element)
      else:
        self.print_str("unknown-tag: " + tag_name)
    else: 
      return ""

  def visit_ci(self, element):
    self.print_str(element.firstChild.data)

  def visit_apply(self, element):
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
    for child in maths.childNodes:
      if maths.nodeType == maths.ELEMENT_NODE:
        self.generic_visit(child)

def format_rate_rule(raterule):
  maths_visitor = ExprVisitor()
  maths = raterule.getElementsByTagName ("math")[0]
  maths_visitor.visit_maths(maths)
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

def outline_model(model, ignore_sources, ignore_sinks):
  """format and print out an outline for the given model"""
  reactions = get_list_of_reactions(model)
  print_amount(len(reactions), "reaction", "reactions")
  for reaction in reactions:
    if ( (ignore_sources and reaction.is_sources()) or
         (ignore_sinks and reaction.is_sink()) ):
      pass
    else:
      print("  " + reaction.format_reaction())

  outline_rate_rules(model)

def outline_sbml_file(filename, ignore_sources, ignore_sinks):
  """parse in a file as an sbml model, extract the outline information
     and then format that information and print it out"""
  dom = xml.dom.minidom.parse(filename)
  model = dom.getElementsByTagName("model")[0]
  outline_model(model, ignore_sources, ignore_sinks)
  
def run():
  """perform the banalities of command-line argument processing and
     then go ahead and calculate the outline for each model file"""
  description = "Print out an outline of an SBML file"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="an sbml file to check invariants for")
  parser.add_argument("--ignore-sources",
                      action="store_true", default=False)
  parser.add_argument("--ignore-sinks",
                      action="store_true", default=False)
 
  arguments = parser.parse_args()

  for filename in arguments.filenames:
    outline_sbml_file(filename,
                      arguments.ignore_sources,
                      arguments.ignore_sinks)


if __name__ == "__main__":
  run()
