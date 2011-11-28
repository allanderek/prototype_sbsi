"""
A module to ease creation of sbml documents
"""

import sys
import xml.dom.minidom as minidom

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


class SBML_Model(object):
  """A class representing an SBML model which could be written
     out as an sbml document
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

    if not compartments:
      default_compartment = document.createElement("compartment")
      default_compartment.setAttribute("id", "default_compartment")
      default_compartment.setAttribute("constant", "true")
      default_compartment.setAttribute("size", "1.0")
      list_of_compartments.appendChild(default_compartment)
    else:
      for compartment in self.compartments:
        compartment_el = document.createElement("compartment")
        compartment_el.setAttribute("id", compartment.id)
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

