"""
A module to ease creation of sbml documents
"""

import sys
import xml.dom.minidom as minidom
import copy

import utils

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

  if arguments.output_file:
    sbml_filename = arguments.output_file
  else:
    sbml_filename = utils.change_filename_ext(filename, ".sbml")
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


  def convert_to_sbml(self):
    """This is a virtual method stub, this method should be overridden
       by any class inheriting from this class"""
    raise NotImplementedError("Expression is really an abstract class")

  def remove_rate_law_sugar(self, reaction):
    """This is a virtual method stub, this method should be overridden
       by any class inheriting from this class. In fact we should be
       doing this with something like a visitor pattern, but I have not
       yet fully groked visitor patterns for python.
       Well it's a bit more than a stub, all the very simple expressions
       which don't have sub-expressions do not need to override this."""
    # pylint: disable-msg=W0613
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


class ApplyExpression(Expression):
  """A class to represent the AST of an apply expression, applying a
     named function to a list of argument expressions"""
  def __init__(self, name, args):
    Expression.__init__(self)
    self.name = name
    self.args = args

  def show_expr(self):
    """Format as a string the application expression"""
    result = self.name + "("
    prefix = ""
    for arg in self.args:
      result += prefix
      result += arg.show_expr()
      prefix = ", "
    result += ")"
    return result

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

  def remove_rate_law_sugar(self, reaction):
    # First apply this to all of the argument expressions.
    new_args = [ arg.remove_rate_law_sugar(reaction) for arg in self.args ]
    self.args = new_args
 
    if self.name == "fMA":
      # Should do some more error checking, eg if there is exactly
      # one argument.
      mass_action_reactants = reaction.get_mass_action_participants()
      extra_args = [ NameExpression(reactant.get_name())
                       for reactant in mass_action_reactants ]
      new_expr = ApplyExpression("times", new_args + extra_args)
      return new_expr
    else:
      new_expr = ApplyExpression(self.name, new_args)
      return new_expr


class VariableDeclaration:
  """A class to represent the AST of a variable declaration in Bio-PEPA"""
  def __init__(self, variable, expression):
    self.variable = variable
    self.expression = expression

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
    if isinstance(self.expression, NumExpression):
      parameter.setAttribute("value", self.expression.show_expr())
    return parameter
 
  def create_initial_assignment(self, document):
    """Creates an sbml element for an initial assignment for the
       parameter corresponding to this variable declaration"""
    init_assign = document.createElement("initialAssignment")
    init_assign.setAttribute("symbol", self.variable)

    math_element = document.createElement("math")
    mathxmlns = "http://www.w3.org/1998/Math/MathML"
    math_element.setAttribute("xmlns", mathxmlns)
    mathxmlnssbml = "http://www.sbml.org/sbml/level3/"
    math_element.setAttribute("xmlns:sbml", mathxmlnssbml)
    init_assign.appendChild(math_element)

    expr_element = self.expression.create_sbml_element(document)
    math_element.appendChild(expr_element)

    return init_assign

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
    self.reactants.append(reactant)

  def add_product(self, product):
    """add a product to the reaction"""
    self.products.append(product)

  def add_modifier(self, modifier):
    """Add a modifier to the reaction"""
    self.modifiers.append(modifier)

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
    return self.reactants and not self.products

  def is_source(self):
    """returns true if the reaction is a source reaction,
       in that it has at least one product and no reactants"""
    return self.products and not self.reactants

  def is_reverse(self, the_inverse):
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
      for l_item in left:
        if l_item not in right:
          return False
      for r_item in right :
        if r_item not in left :
          return False
      return True
    if (equal_lists(self.reactants, the_inverse.products) and
        equal_lists(self.products, the_inverse.reactants)):
      return True
    else:
      return False

  def reverse_reaction(self):
    """Creates a copy of this reaction but reversed"""
    reversed_reaction = copy.copy(self)
    # In theory actually we should copy these, rather than simply
    # assign, on the basis that if we changed one do we want the
    # other to also change?
    reversed_reaction.reactants = self.products
    reversed_reaction.products = self.reactants
    reversed_reaction.kinetic_law = self.reverse_kinetic_law
    reversed_reaction.reverse_kinetic_law = self.kinetic_law

    return reversed_reaction

  def involves(self, species):
    """Returns true if the given species is involved with this reaction"""
    name = species.name
    reactant_names = [ r.get_name() for r in self.reactants ] 
    product_names = [ p.get_name() for p in self.products ]
    modifier_names = [ m.get_name() for m in self.modifiers ]
    if ( (name in reactant_names or
          name in product_names  or
          name in modifier_names) and
         (species.compartment == self.location) ):
      return True
    return False

  def format_reaction(self):
    """Return a string representing the reaction in a format 
       suitable for human consumption"""
    results = self.name 
    if self.location:
      results += "@" + self.location

    results += ": "
    # so the first reactant has nothing attached to the front of it
    reactants_and_modifiers = self.reactants + self.modifiers
    results += ", ".join([ r.format_participant() 
                           for r in reactants_and_modifiers])
    results += " --> "
    results += ", ".join([ p.format_participant() 
                           for p in self.products])
 
    return results

  def create_element(self, document):
    """Create an xml element representing this reaction"""
    reaction_element = document.createElement("reaction")
    reaction_element.setAttribute("id", self.name)
    reaction_element.setAttribute("reversible", "false")
    reaction_element.setAttribute("fast", "false")
    if self.location:
      reaction_element.setAttribute("compartment", self.location)
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
      math_element = document.createElement("math")
      kinetic_law.appendChild(math_element)
      mathxmlns = "http://www.w3.org/1998/Math/MathML"
      math_element.setAttribute("xmlns", mathxmlns)
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
  pass


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
    
  def create_compartment_elements(self, document, model_element):
    """Create the elements to describe the compartments within
       the sbml model
    """
    list_of_compartments = document.createElement("listOfCompartments")
    model_element.appendChild(list_of_compartments)

    if not self.compartments:
      default_compartment = document.createElement("compartment")
      default_compartment.setAttribute("id", "default_compartment")
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
        if not compartment:
          compartment = "default_compartment"
        species_element.setAttribute("compartment", compartment)
        if comp_def.initial_amount:
          species_element.setAttribute("initialAmount", 
                                       comp_def.initial_amount)
        else:
          species_element.setAttribute("initialAmount", "0")
          
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
   
      init_assigns = document.createElement("listOfInitialAssignments")
      top_element.appendChild(init_assigns)
      for var_dec in self.var_decs:
        param_element = var_dec.create_parameter_element(document)
        list_of_params.appendChild(param_element)

        init_assign = var_dec.create_initial_assignment(document)
        init_assigns.appendChild(init_assign)


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

    list_of_reactions = document.createElement("listOfReactions")
    model_element.appendChild(list_of_reactions)

    if self.reactions:
      for reaction in self.reactions:
        reaction_element = reaction.create_element(document)
        list_of_reactions.appendChild(reaction_element)

    return document
