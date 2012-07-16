"""
A module to ease creation of sbml documents
"""

import sys
import xml.dom.minidom as minidom
import copy

import utils

default_location_name = "default_compartment"


class SBMLConfiguration(object):
  """A class representing the SBML configuration which may be used to
     configure the input/output of SBML models
  """
  def __init__(self):
    self.sbml_level = 3
    self.sbml_level_version = 1
    self.copasi_spec_ref_workaround = False


def output_to_sbml_file(filename, arguments, document):
  """A utility function for consumers of this module to calculate
     the output file for given a set of parsed arguments which includes
     the possible argument "--output-file". We then output the given
     document to that file. This is basically a suitable place to put
     this functionality that is common to nearly all x_to_sbml translators
     that we have.
     Note that the first argument is the filename of the file from which
     we have translated into sbml and not the sbml filename itself.
  """
  sbml_filename = utils.get_output_filename(filename, arguments, ".sbml")
  if sbml_filename == "stdout":
    sbml_file = sys.stdout
  else:
    sbml_file =  open(sbml_filename, "w")

  document.writexml(sbml_file, 
                    encoding="UTF-8",
                    indent="",
                    addindent="  ",
                    newl="\n")
  if sbml_filename != "stdout":
    sbml_file.close()

 
class Expression:
  """The base class for all classes which represent the AST of some
     kind of expression"""
  def __init__(self):
    pass

  def show_expr(self):
    """This is a virtual method stub, this method should be overridden
       by any class inheriting from this class"""
    raise NotImplementedError("Expression is really an abstract class")

  def used_names(self):
    """This is a virtual method stud, this method should be overridden
       by any class inheriting from this class. The overriding method
       should return a set of names used within the expression
    """
    raise NotImplementedError("Expression is really an abstract class")

  def convert_to_sbml(self):
    """This is a virtual method stub, this method should be overridden
       by any class inheriting from this class"""
    raise NotImplementedError("Expression is really an abstract class")

  def get_value(self, environment=None):
    """Returns the underlying value of this expression. For complex
       expressions a dictionary mapping names to values may be supplied.
       We return None, the value cannot be derived, generally this will
       mean that it uses a name which is not defined by the provided
       variabile dictionary (or one was not provided).
    """
    # pylint: disable=W0613
    # pylint: disable=R0201
    return None

  def munge_names(self, function):
    """Munges the names used within the expression using the function
       supplied. This is a virtual method stud, see below on our comment
       of the remove_rate_law_sugar method. Essentially I think I should
       be able to do something much nicer, using a visitor pattern.
       Again here this is a bit more than a stub since all the simple
       expressions which cannot contain any names do not need to override
       this stub implementation.
    """
    # pylint: disable=W0613
    # pylint: disable=R0201
    return None


  def remove_rate_law_sugar(self, reaction=None):
    """This is a virtual method stub, this method should be overridden
       by any class inheriting from this class. In fact we should be
       doing this with something like a visitor pattern, but I have not
       yet fully groked visitor patterns for python.
       Well it's a bit more than a stub, all the very simple expressions
       which don't have sub-expressions do not need to override this."""
    # pylint: disable=W0613
    return self 


class NumExpression(Expression):
  """A class to represent the AST of an number literal expression"""
  def __init__(self, number):
    Expression.__init__(self)
    self.number = number
  def convert_to_sbml(self):
    """Convert the number expression to SBML code for the expression"""
    return "<cn>" + str(self.number) + "</cn>"
  def show_expr(self):
    """Display the underlying number of the numerical expression"""
    return str(self.number)

  def get_value(self, environment=None):
    """Returns the underlying value of this expression"""
    return self.number

  def used_names(self):
    """Return the set of used names, clearly here there are none"""
    return []
  def create_sbml_element(self, document):
    """Create an sbml xml element for the sbml code for the expression"""
    cn_element = document.createElement("cn")
    cn_text = document.createTextNode(str(self.number))
    cn_element.appendChild(cn_text)
    return cn_element


class NameExpression(Expression):
  """A class to represent the AST of a variable (name) expression"""
  def __init__(self, name):
    Expression.__init__(self) 
    self.name = name

  def show_expr(self):
    """Format as a string the name expression"""
    return self.name

  def get_value(self, environment=None):
    """Evalutes this expression based on the given variable_dictionary,
       returns None if the variable diction is absent or does not define
       the name used in this expression.
    """
    if environment:
      return environment[self.name]
    else:
      return None
      

  def used_names(self):
    """Return the set of names used within this expression"""
    return set([self.name])

  def munge_names(self, function):
    self.name = function(self.name)

  def convert_to_sbml(self):
    """Convert the variable(name) expression to SBML code for
       the expression"""
    return "<ci>" + self.name + "</ci>"

  def create_sbml_element(self, document):
    """Create an sbml xml element for the sbml code for the expression"""
    ci_element = document.createElement("ci")
    ci_text = document.createTextNode(self.name)
    ci_element.appendChild(ci_text)
    return ci_element


def show_apply_expression(function_name, children):
  """Formats an apply expression as a string and returns that string.
     Checks for common arithmetic operators and outputs the appropriate
     infix expression in the case that it finds one.
  """
  function_dict = { "plus" : "+", 
                    "minus" : "-",
                    "divide" : "/",
                    "times" : "*",
                    "power" : "^",
                  }
  # The check on the length of children is just in case someone
  # has managed to say apply 'times' to no arguments which would
  # otherwise cause an error when we attempt to print the first one.
  # It's unclear what we should do in that case, but for now I fall
  # through to the generic case and basically you'll end up with
  # just the 'times' (named as 'times' not as *) printed out.
 
  result = ""
 
  if function_name in function_dict and len(children) > 1 :
    result += "("
    # Could just put the spaces in the dictionary above?
    operator = " " + function_dict[function_name] + " "
    result += operator.join(children)
    result += ")"
  else:
    result += function_name + "("
    result += ", ".join(children) 
    result += ")"
 
  return result



class ApplyExpression(Expression):
  """A class to represent the AST of an apply expression, applying a
     named function to a list of argument expressions"""
  def __init__(self, name, args):
    Expression.__init__(self)
    self.name = name
    self.args = args

  def show_expr(self):
    """Format as a string the application expression"""
    arg_strings = [ arg.show_expr() for arg in self.args ]
    return show_apply_expression(self.name, arg_strings)

  def used_names(self):
    """Return the set of names used within this apply expression"""
    result_set = set()
    for expr in self.args:
      result_set = result_set.union(expr.used_names())
    return result_set
    
  def munge_names(self, function):
    """Must munge all the names, we do not munge the name of the
       function of the apply expression however.
    """
    for child in self.args:
      child.munge_names(function)

  def convert_to_sbml(self):
    """return a string representing the sbml of an math apply expression"""
    result = "<apply>\n"
    result += "  <" + self.name + "/>\n"
    for argument in self.args:
      result += "  " + argument.convert_to_sbml() + "\n"
    result += "</apply>"
    return result

  def create_sbml_element(self, document):
    """Create an sbml xml element for the sbml code for the expression"""
    apply_el = document.createElement("apply")
    operator_el = document.createElement(self.name)
    apply_el.appendChild(operator_el)
    for argument in self.args:
      arg_el = argument.create_sbml_element(document)
      apply_el.appendChild(arg_el)
    
    return apply_el

  def remove_rate_law_sugar(self, reaction=None):
    # First apply this to all of the argument expressions.
    new_args = [ arg.remove_rate_law_sugar(reaction) for arg in self.args ]
    self.args = new_args
 
    if reaction != None and self.name == "fMA":
      # Should do some more error checking, eg if there is exactly
      # one argument.
      mass_action_reactants = reaction.get_mass_action_participants()
      extra_args = [ NameExpression(reactant.get_name())
                       for reactant in mass_action_reactants ]
      all_args = new_args + extra_args
      # fMA should have exactly one argument, the additional arguments
      # are the populations of all the reactants/modifiers of the reaction.
      # It could be that there are no such, in otherwords we have a
      # source reaction
      if len(all_args) > 1:
        new_expr = ApplyExpression("times", new_args + extra_args)
        return new_expr
      else:
        # If there is only the original argument then just return that
        # even without the surrounding 'fMA' application.
        return all_args[0]
    # I'm not really comfortable having this here, it's really only for
    # COPASI and should be a separate and configurable expression
    # transformer.
    elif self.name == "sqrt":
      if len(new_args) != 1:
        # This is a terrible way to error as well!
        print ("sqrt function takes only a single argument")
        sys.exit(1)
      argument = new_args[0]
      half_expr = NumExpression(0.5)
      new_expr = ApplyExpression("power", [argument, half_expr])
      return new_expr
    else:
      new_expr = ApplyExpression(self.name, new_args)
      return new_expr


class FunctionDefinition(object):
  """A class to represent a function definition in SBML"""
  def __init__(self):
    self.name = None
    self.parameters = []
    self.body = None

  def format_fundef(self):
    """A method to format the function definition as a string"""
    params_string = "(" + ", ".join(self.parameters) + ")"
    body_string = ""
    if self.body != None:
      body_string = self.body.show_expr()
    result = " ".join ([ "fundef", params_string, "=", body_string ])
    return result

class VariableDeclaration(object):
  """A class to represent the AST of a variable declaration in Bio-PEPA"""
  def __init__(self, variable, expression):
    self.variable = variable
    self.expression = expression
    self.init_assign = None

  def get_name(self):
    """Return the name of the variable being declared"""
    return self.variable

  def create_parameter_element(self, document):
    """Creates a parameter sbml element for the parameter corresponding
       to this variable declaration"""
    parameter = document.createElement("parameter")
    parameter.setAttribute("id", self.variable)
    parameter.setAttribute("name", self.variable)
    # true means that the parameter's value can only be set by an
    # initial assignment (or here as the 'value' attribute) so I think it
    # is reasonable here to set it to true.
    parameter.setAttribute("constant", "true")
    number = self.expression.get_value()

    # I would quite like it to be configurable such that we could
    # ALWAYS create an initial assignment, not just if there is no
    # number value.
    if number != None:
      parameter.setAttribute("value", str(number))
    else:
      self.init_assign = create_initial_assignment(document,
                                                   self.variable,
                                                   self.expression)
          
    return parameter
 
  # def create_initial_assignment(self, document):
  #   """Creates an sbml element for an initial assignment for the
  #      parameter corresponding to this variable declaration"""
  #   return create_initial_assignment(document, 
  #                                    self.variable,
  #                                    self.expression)

def create_math_element(document):
  """Creates a math element into which is generally stored some kind
     of expression, for example it a math element is the child of a
     kineticLaw element of a Reaction element.
  """
  math_element = document.createElement("math")
  mathxmlns = "http://www.w3.org/1998/Math/MathML"
  math_element.setAttribute("xmlns", mathxmlns)
  mathxmlnssbml = "http://www.sbml.org/sbml/level3/"
  math_element.setAttribute("xmlns:sbml", mathxmlnssbml)
  return math_element
  

def create_initial_assignment(document, name, expression):
  """Create an initial assignment element which assigns the given
     name to the given expression
  """
  init_assign = document.createElement("initialAssignment")
  init_assign.setAttribute("symbol", name)

  math_element = create_math_element(document)
  init_assign.appendChild(math_element)

  expression = expression.remove_rate_law_sugar()
  expr_element = expression.create_sbml_element(document)
  math_element.appendChild(expr_element)

  return init_assign


class IdNamedElement(object):
  """Many of the classes below are a parsed representation of a
     kind of sbml element. Many such elements have a required id
     attribute and an optional name element. This is a base class
     which will define and parse those fields such that we needn't
     do it ourselves for each kind of element.
  """
  def __init__(self, element):
    self.ident = element.getAttribute("id")
    name = element.getAttribute("name") 
    if name:
      self.name = name 
    else:
      self.name = self.ident

class Species(IdNamedElement):
  """A simple class to represent a species within an sbml model"""
  def __init__(self, element):
    super(Species, self).__init__(element)
    self.compartment = element.getAttribute("compartment")
    self.initial_amount = element.getAttribute("initialAmount")
    self.initial_conc = element.getAttribute("initialConcentration")
    self.substance_units = element.getAttribute("substanceUnits")
    self.has_only_substance_units = element.getAttribute(
          "hasOnlySubstanceUnits")
    self.boundary_condition = element.getAttribute("boundaryCondition")
    self.constant = element.getAttribute("constant")
    self.conversion_factor = element.getAttribute("conversionFactor")

  def get_name(self):
    """ return the name of the species"""
    return self.name

  def format_species(self):
    """Returns a reasonable format of this species definition"""
    return self.name + " in " + self.compartment


class Parameter(IdNamedElement):
  """A simple class to represent a parameter definition
     within an sbml model
  """
  def __init__(self, element):
    super(Parameter, self).__init__(element)
    self.value = element.getAttribute("value")
    self.units = element.getAttribute("units")
    self.boolean = element.getAttribute("boolean")

  def format_param(self):
    """Format this parameter as an assignment if there is a value and
       just as the name otherwise"""
    result = self.name
    if self.value != None:
      result += " = " + str(self.value)
    return result

class ReactionParticipant:
  """A simple class to represent a reaction participant"""
  def __init__(self, name, stoich=1):
    self.name = name
    self.stoich = stoich
    # If more members are added here we should consider how this
    # affects Reaction.canonicalise_participants.

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


class Reaction(object):
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
    #  Because of create_element, whatever you set this to,
    # should have a 'remove_rate_law_sugar' method and the result
    # of that call should be an object (possibly the same one) with
    # a method for creating an sbml element to represent the expression
    self.kinetic_law = None
    self.reverse_kinetic_law = None
    self.location = None

  def add_reactant(self, reactant):
    """add a reactant into the reaction definition"""
    if isinstance(reactant, ReactionParticipant):
      self.reactants.append(reactant)
    else:
      self.reactants.append(ReactionParticipant(reactant, stoich=1))

  def add_product(self, product):
    """add a product to the reaction"""
    if isinstance(product, ReactionParticipant):
      self.products.append(product)
    else:
      self.products.append(ReactionParticipant(product, stoich=1))

  def add_modifier(self, modifier):
    """Add a modifier to the reaction"""
    if isinstance(modifier, ReactionParticipant):
      self.modifiers.append(modifier)
    else:
      self.modifiers.append(ReactionParticipant(modifier, stoich=1))

  def canonicalise_participants(self):
    """We look at the reaction participants and decide whether we could
       more succintly represent them. For example if 'A' appears as
       both a single reactant and product then it can be instead
       represented as a modifier. Additionally if 'A' appears twice as
       a reactant then we can simply increase its stoichiometry
    """
    # Note that as we are creating new 'ReactionParticipants' then this
    # is a little helicoptor code, in that if change ReactantParticipant
    # then we would need to change this as well. 
    name_dict = dict()
    for participant in self.reactants:
      if participant.name in name_dict:
        name_dict[participant.name] += participant.stoich
      else:
        name_dict[participant.name] = participant.stoich
    for participant in self.products:
      if participant.name in name_dict:
        name_dict[participant.name] -= participant.stoich
      else:
        name_dict[participant.name] = -(participant.stoich)
    # For modifiers then we essentially just make sure that the name
    # is in the dictionary such that if it is unchanged it will be added
    # in the final part of this method.
    for participant in self.modifiers:
      if participant.name not in name_dict:
        name_dict[participant.name] = 0
    

    # Note that there is no stoichiometry on modifier species so for now
    # can essentially get away with ignoring the stoichiometry on
    # modifier species references.
    self.reactants = []
    self.products = []
    self.modifiers = []
    for name, value in name_dict.iteritems():
      if value == 0:
        participant = ReactionParticipant(name, value)
        self.modifiers.append(participant)
      elif value < 0:
        participant = ReactionParticipant(name, - value)
        self.products.append(participant)
      else:
        participant = ReactionParticipant(name, value)
        self.reactants.append(participant)


  def get_mass_action_participants(self):
    """Return the left hand side participants which contribute to the
       rate in a mass action rate method. For example A + B --> C 
       would return A and B. Essentially this is so that fMA(r) could
       be translated to r * A * B.
    """
    # I think we return all the modifiers but perhaps not the inhibitors?
    return self.reactants + self.modifiers

  def is_sink(self):
    """returns true if the reaction is a sink,
       in that it has at least one reactant and no products"""
    return bool (self.reactants and not self.products)

  def is_source(self):
    """returns true if the reaction is a source reaction,
       in that it has at least one product and no reactants"""
    return bool(self.products and not self.reactants)

  def is_reverse(self, the_inverse):
    """Returns true if the given reaction is the reverse of
       the current reaction
    """
    if (utils.equal_lists(self.reactants, the_inverse.products) and
        utils.equal_lists(self.products, the_inverse.reactants)):
      return True
    else:
      return False

  def reverse_reaction(self):
    """Creates a copy of this reaction but reversed"""
    reversed_reaction = copy.copy(self)
    reversed_reaction.name = reversed_reaction.name + "_rev"
    # In theory actually we should copy these, rather than simply
    # assign, on the basis that if we changed one do we want the
    # other to also change?
    reversed_reaction.reactants = self.products
    reversed_reaction.products = self.reactants
    reversed_reaction.kinetic_law = self.reverse_kinetic_law
    reversed_reaction.reverse_kinetic_law = self.kinetic_law

    return reversed_reaction

  def is_in_species_list(self, species_list, species):
    """A helper method for is_reactant, is_modifier and is_product.
       Essentially checks if the given species_list contains the given
       species and that that species is in the correction location
    """
    # So if the species given is just a name, just check for that
    if isinstance(species, str):
      for species_item in species_list:
        if species_item.name == species:
          return True
      return False
    # If it's not just a name, assume it is a species with a compartment
    # and hence check whether it is in the correct location
    for species_item in species_list:
      if (species.name == species_item.name and
          species.compartment == self.location):
        return True
    return False

  def is_reactant(self, species):
    """Returns true if the given species is a reactant of this reaction
    """
    return self.is_in_species_list(self.reactants, species)

  def is_modifier(self, species):
    """Returns true if the given species is a product of this reaction
    """
    return self.is_in_species_list(self.modifiers, species)  

  def is_product(self, species):
    """Returns true if the given species is a product of this reaction
    """
    return self.is_in_species_list(self.products, species)


  def involves(self, species):
    """Returns true if the given species is involved with this reaction"""
    return (self.is_reactant(species) or
            self.is_modifier(species) or
            self.is_product(species))

  def format_reaction(self):
    """Return a string representing the reaction in a format 
       suitable for human consumption"""
    results = self.name 
    if self.location:
      results += "@" + self.location

    results += ": "
    reactant_names = [ r.format_participant() for r in self.reactants ]
    modifier_names = [ "$" + r.format_participant() 
                        for r in self.modifiers]
    results += ", ".join(reactant_names + modifier_names)
                          
    results += " --> "

    product_names = [ r.format_participant() for r in self.products ]
    results += ", ".join(product_names)

    if self.is_source():
      results += "  (is source)"
    if self.is_sink():
      results += "  (is sink)"
 
    return results

  def create_element(self, document, sbml_configuration):
    """Create an xml element representing this reaction"""
    reaction_element = document.createElement("reaction")
    reaction_element.setAttribute("id", self.name)
    reaction_element.setAttribute("reversible", "false")
    reaction_element.setAttribute("fast", "false")
    if self.location:
      reaction_element.setAttribute("compartment", self.location)
    else:
      reaction_element.setAttribute("compartment", default_location_name)
    if self.reactants:
      list_of_reactants = document.createElement("listOfReactants")
      reaction_element.appendChild(list_of_reactants)
      for reactant in self.reactants:
        spec_ref = document.createElement("speciesReference")
        # spec_ref.setAttribute("id", reactant.name)
        spec_ref.setAttribute("species", reactant.name)
        # for biopepa models the stoichiometry values cannot
        # change during the simulation so this 'constant' attribute is
        # always true.
        if not sbml_configuration.copasi_spec_ref_workaround:
          spec_ref.setAttribute("constant", "true")
        spec_ref.setAttribute("stoichiometry", str(reactant.stoich))
        list_of_reactants.appendChild(spec_ref)
    if self.products:
      list_of_products = document.createElement("listOfProducts")
      reaction_element.appendChild(list_of_products)
      for product in self.products:
        spec_ref = document.createElement("speciesReference")
        # spec_ref.setAttribute("id", product.name)
        spec_ref.setAttribute("species", product.name)
        # for biopepa models the stoichiometry values cannot
        # change during the simulation so this 'constant' attribute is
        # always true.
        if not sbml_configuration.copasi_spec_ref_workaround:
          spec_ref.setAttribute("constant", "true")
        spec_ref.setAttribute("stoichiometry", str(product.stoich))
        list_of_products.appendChild(spec_ref)
    if self.modifiers:
      list_of_modifiers = document.createElement("listOfModifiers")
      reaction_element.appendChild(list_of_modifiers)
      for modifier in self.modifiers:
        spec_ref = document.createElement("modifierSpeciesReference")
        # spec_ref.setAttribute("id", modifier.name)
        # Note there is no stoichiometry attribute on
        # modifierSpeciesReference elements in SBML.
        spec_ref.setAttribute("species", modifier.name)
        list_of_modifiers.appendChild(spec_ref)
    if self.kinetic_law:
      kinetic_law = document.createElement("kineticLaw")
      reaction_element.appendChild(kinetic_law)
      math_element = create_math_element(document)
      kinetic_law.appendChild(math_element)
      simplified_rate = self.kinetic_law.remove_rate_law_sugar(self)
      expr_element = simplified_rate.create_sbml_element(document)
      math_element.appendChild(expr_element)
    return reaction_element


class Assignment:
  """A class representing an assignment of an sbml expression to
     a variable name"""
  def __init__(self, variable, expression):
    self.variable = variable
    self.expression = expression

  def get_variable_name(self):
    """return the name of the variable being assigned to"""
    return self.variable

  def get_assigned_expr(self):
    """return the expression part of the assignment, that is
       the right hand side"""
    return self.expression

class InitialAssignment(Assignment):
  """A class representing an sbml initial assignment"""
  pass

class AssignmentRule(Assignment):
  """A class representing an sbml assignment rule"""
  def create_element(self, document):
    """Create an assignmentRule element representing this assignment rule
    """
    assign_rule = document.createElement("assignmentRule")
    assign_rule.setAttribute("variable", self.variable)

    math_element = create_math_element(document)
    assign_rule.appendChild(math_element)

    expr_element = self.expression.create_sbml_element(document)
    math_element.appendChild(expr_element)

    return assign_rule



def find_element_after(top_element, candidate_names):
  """Occasionally we need to insert an element in a given position
     but the xml dom library only provides us with insertBefore.
     Therefore we need to find the element that occurs directly after
     the position we want. Sometimes those that follow are optional.
     For example we wish to insert the listOfInitialAssignments, before
     listOfRules, but there may not be any rules and hence no
     listOfRules element, so we need to search for the next candidate
     extra. So this method takes in an element and a list of candidiate
     child elements and returns the first one that exists.
  """
  for candidate in candidate_names:
    elements = top_element.getElementsByTagName(candidate)
    if elements:
      return elements[0]

  # If none of the candidates are found, then there is no choice but
  # to return None. 
  return None

class SBMLModel(object):
  """A class representing an SBML model which could be written
     out as an sbml document. This is primarily for writing the document
     out rather than parsing it in. In particular we have no notion of
     such things as assignments or initial assignments, these are
     garnered automatically from the species, compartments and variable
     declarations.
  """
  def __init__(self):
    self.compartments = None
    self.reactions = None
    self.component_defs = None
    self.var_decs = None
    self.assign_rules = None
    self.init_assign_elements = []
    self.sbml_configuration = SBMLConfiguration()
    
  def create_compartment_elements(self, document, model_element):
    """Create the elements to describe the compartments within
       the sbml model
    """
    list_of_compartments = document.createElement("listOfCompartments")
    model_element.appendChild(list_of_compartments)

    if not self.compartments:
      default_compartment = document.createElement("compartment")
      default_compartment.setAttribute("id", default_location_name)
      default_compartment.setAttribute("constant", "true")
      default_compartment.setAttribute("size", "1.0")
      list_of_compartments.appendChild(default_compartment)
    else:
      for compartment in self.compartments:
        compartment_el = document.createElement("compartment")
        compartment_el.setAttribute("id", compartment.ident)
        compartment_el.setAttribute("constant", "true")
        if compartment.size:
          compartment_el.setAttribute("size", compartment.size)
        else:
          compartment_el.setAttribute("size", "1.0")
        list_of_compartments.appendChild(compartment_el)
   
  def create_species_elements(self, document, model_element):
    """Creates the listOfSpecies element and all of the associated
       species elements under it for a set of component definitions"""
    list_of_species = document.createElement("listOfSpecies")
    model_element.appendChild(list_of_species)
    if self.component_defs:
      for comp_def in self.component_defs:
        species_element = document.createElement("species")
        name = comp_def.name
        species_element.setAttribute("id", name)
        species_element.setAttribute("name", name)
        compartment = comp_def.location
        expr = comp_def.initial_expression
        if not compartment:
          compartment = default_location_name
        species_element.setAttribute("compartment", compartment)
        if comp_def.initial_amount:
          species_element.setAttribute("initialAmount", 
                                       comp_def.initial_amount)
        else:
          species_element.setAttribute("initialAmount", "0")
          if expr:
            init_assign = create_initial_assignment(document, name, expr)
            self.init_assign_elements.append(init_assign)
          
        species_element.setAttribute("hasOnlySubstanceUnits", "true")
        species_element.setAttribute("constant", "false")
        species_element.setAttribute("boundaryCondition", "false")
        list_of_species.appendChild(species_element) 

  def convert_variable_declarations(self, document, top_element):
    """Convert a set of variable declarations into SBML elements and
       add them to an sbml document
    """
    if self.var_decs:
      list_of_params = document.createElement("listOfParameters")
      top_element.appendChild(list_of_params)
 
      # Variable declarations are often used inside things such as
      # initial assignments for species, therefore we take care to
      # put the variable assigns first. However we don't just insert
      # each on at the head of the initial assignments since that would
      # reverse the order of the variable declarations themselves which
      # may be in a specific order. 
      init_assigns = [] 
      for var_dec in self.var_decs:
        param_element = var_dec.create_parameter_element(document)
        list_of_params.appendChild(param_element)

        if var_dec.init_assign != None:
          init_assigns.append(var_dec.init_assign)
     
      # So here we add the variable declaration initial assignments
      # in order to the front of the overall initial assignments 
      self.init_assign_elements = init_assigns + self.init_assign_elements


  def create_sbml_document(self):
    """Create an sbml model document for this SBML model"""
                          
    xml_implementation = minidom.getDOMImplementation()
    document = xml_implementation.createDocument(None, "sbml", None)
    top_element = document.documentElement 
    xmlns = "http://www.sbml.org/sbml/level3/version1/core" 
    top_element.setAttribute("xmlns", xmlns)
    top_element.setAttribute("level", "3")
    top_element.setAttribute("version", "1")

    model_element = document.createElement("model")
    top_element.appendChild(model_element)

    self.create_compartment_elements(document, model_element)
    self.create_species_elements(document, model_element)
    self.convert_variable_declarations(document, model_element)

    if self.assign_rules:
      assign_rules = document.createElement("listOfRules")
      for assign_rule in self.assign_rules:
        assign_rule_element = assign_rule.create_element(document)
        assign_rules.appendChild(assign_rule_element)
      model_element.appendChild(assign_rules)

    if self.reactions:
      list_of_reactions = document.createElement("listOfReactions")
      for reaction in self.reactions:
        reaction_element = reaction.create_element(document,
                                                   self.sbml_configuration)
        list_of_reactions.appendChild(reaction_element)
      model_element.appendChild(list_of_reactions)

    if self.init_assign_elements: 
      init_assigns = document.createElement("listOfInitialAssignments")
      for init_assign in self.init_assign_elements:
        init_assigns.appendChild(init_assign)
      model_element.appendChild(init_assigns) 
    
      after_init_assigns = find_element_after(model_element,
                                              [ "listOfRules",
                                                "listOfConstraints",
                                                "listOfReactions",
                                                "listOfEvents"
                                              ])
      model_element.insertBefore(init_assigns, after_init_assigns)

    return document

