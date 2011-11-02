"""a simple script to give a simple outline of an sbml model"""
import xml.dom.minidom
import argparse

class ReactionParticipant:
  """A simple class to represent a reaction participant"""
  def __init__(self, name, stoich=1):
    self.name = name
    self.stoich = stoich

  def __cmp__(self, other):
    """Implementation of comparison for general use"""
    if self.name == other.name and self.stoich == other.stoich:
      return 0
    elif self.name == other.name:
      return cmp(self.stoich, other.stoich)
    else:
      return cmp(self.name, other.name)

  def get_name(self):
    """Return the name of the reaction participant"""
    return self.name
  def get_stoichiometry(self):
    """Return the stoichiometry of this reaction participant
       within the reaction"""
    return self.stoich

  def format_participant(self):
    """Returns a reasonable format of this reaction participant,
       basically just the name if the stoichiometry is 1, or
       (name, stoich) otherwise
    """
    if self.stoich == 1:
      return self.get_name()
    else:
      return "(" + self.get_name() + ", " + str(self.stoich) + ")"

class Species:
  """A simple class to represent a species within an sbml model"""
  def __init__(self, name, compartment):
    self.name = name
    self.compartment = compartment
  def get_name(self):
    """ return the name of the species"""
    return self.name

  def format_species(self):
    """Returns a reasonable format of this species definition"""
    return self.name + " in " + self.compartment

class Reaction:
  """A class which represents a reaction"""
  def __init__(self, name):
    """initialise a reaction which has no reactants or products
       these may in turn be added with 'add_reactant' and 'add_product'
       Similarly there are initially no modifiers which can be added with
       'add_modifier'
    """
    self.reactants = []
    self.products = []
    self.modifiers = []
    self.name = name
    self.kinetic_law = None

  def get_name(self):
    """return the name of the reaction"""
    return self.name
  
  def add_reactant(self, reactant):
    """add a reactant into the reaction definition"""
    self.reactants.append(reactant)

  def add_product(self, product):
    """add a product to the reaction"""
    self.products.append(product)

  def add_modifier(self, modifier):
    """Add a modifier to the reaction"""
    self.modifiers.append(modifier)

  def get_modifiers(self):
    """return the list of modifier names"""
    return self.modifiers

  def get_reactants(self):
    """return the list of reactant names"""
    return self.reactants

  def get_products(self):
    """return the list of product names"""
    return self.products

  def get_kinetic_law(self):
    """Return the kinetic law for the rate of this reaction if any"""
    return self.kinetic_law

  def set_kinetic_law(self, kinetic_law):
    """Set the kinetic law, for the rate of this reaction"""
    self.kinetic_law = kinetic_law

  def is_sink(self):
    """returns true if the reaction is a sink,
       in that it has at least one reactant and no products"""
    return self.reactants and not self.products

  def is_source(self):
    """returns true if the reaction is a source reaction,
       in that it has at least one product and no reactants"""
    return self.products and not self.reactants

  def is_reverse(self, reversed):
    """Returns true if the given reaction is the reverse of
       the current reaction
    """
    def equal_lists(left, right):
      """Checks if two lists are 'equal', equal if they are considered
         to be sets
      """
      # Could check the lengths here, but I think we want
      # 'l,l,r' to be equal to 'l,r', so we're ignoring duplicates.
      # Not sure if that's correct to do so though.
      for l in left:
        if l not in right:
          return False
      for r in right :
        if r not in left :
          return False
      return True
    if (equal_lists(self.reactants, reversed.products) and
        equal_lists(self.products, reversed.reactants)):
       return True
    else:
       return False

  def involves(self, species):
    """Returns true if the given species is involved with this reaction"""
    name = species.name
    reactant_names = [ r.get_name() for r in self.reactants ] 
    product_names = [ p.get_name() for p in self.products ]
    modifier_names = [ m.get_name() for m in self.modifiers ]
    if (name in reactant_names or
        name in product_names  or
        name in modifier_names):
      return True
    return False

  def format_reaction(self):
    """Return a string representing the reaction in a format 
       suitable for human consumption"""
    results = self.name + ": "

    # so the first reactant has nothing attached to the front of it
    reactants_and_modifiers = self.reactants + self.modifiers
    results += ", ".join([ r.format_participant() 
                           for r in reactants_and_modifiers])
    results += " --> "
    results += ", ".join([ p.format_participant() 
                           for p in self.products])
 
    return results


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
    return ReactionParticipant(name)
  else :
    return ReactionParticipant(name, stoich=float(stoich))
  return name

def get_reaction_of_element(reaction_element):
  """a function to return a reaction object from a reaction sbml element"""
  name = reaction_element.getAttribute("id")
  reaction = Reaction(name)
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

  kinetic_laws = reaction_element.getElementsByTagName("kineticLaw")
  if kinetic_laws:
    reaction.set_kinetic_law(kinetic_laws[0])

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


class Assignment:
  def __init__(self, variable, expression):
    self.variable = variable
    self.expression = expression

  def get_variable_name(self):
    return self.variable

  def get_assigned_expr(self):
    return self.expression

class InitialAssignment(Assignment):
  pass

class AssignmentRule(Assignment):
  pass

def get_assignment_rule_of_element(assign_rule_element):
  variable_name = assign_rule_element.getAttribute("variable")
  math_elements = assign_rule_element.getElementsByTagName("math")
  expression = None
  if math_elements:
    expression = math_elements[0]

  return AssignmentRule(variable_name, expression)
 
def get_list_of_assignment_rules(model):
  """Return the list of assignment rules from an sbml model"""
  return get_elements_from_lists_of_list("listOfRules",
                                         "assignmentRule",
                                         get_assignment_rule_of_element,
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
