"""
A module that implements a parser for the Bio-PEPA language
"""
import argparse
import parcon
from parcon import Forward, InfixExpr, Translate, Optional, ZeroOrMore
import xml.dom.minidom as minidom


# A simply utility for creating parsers which accept a list of
# somethings, separated by something elses.
# Note here that the separator_parser should return None
def create_separated_by(element_parser, separator_parser):
  """A utility to create a parser for a list of elements which are
     separated by a given parser. This is useful for writing things
     such as a comma-separated list of arguments
  """
  def create_list_result(parse_result):
    """Simple generator function to create a list from a parser result
       of the list_syntax
    """
    first_item = parse_result[0]
    rest_items = parse_result[1]
    rest_items.insert(0, first_item)
    return rest_items
  # Could call 'parcon.Discard' on the separator parser to ensure that
  # it returns None.
  following_parser = separator_parser + element_parser  
  list_syntax = element_parser + Optional(ZeroOrMore(following_parser))
  return Translate(list_syntax, create_list_result)


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

def create_number_expression(number_str):
  """Simple utility to create a number expression from the parse result
     of parsing a simple number"""
  return NumExpression(float(number_str))
 

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

 
def make_name_expression(parse_result):
  """Simple post parsing creation method for NameExpression"""
  return NameExpression(parse_result)

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

   
def make_apply_expression(parse_result):
  """simple post-parsing creation method for apply expressions"""
  if len(parse_result) == 1:
    # We assume then that there are no arguments, and that it was
    # an application expression such as: f()
    return ApplyExpression(parse_result[0], [])
  # Otherwise the second part of the parse result should be a list
  # of arguments
  return ApplyExpression(parse_result[0], parse_result[1])
def make_plus_expression(left, right):
  """simple post-parsing creation method for add expressions"""
  return ApplyExpression("plus", [left, right])
def make_minus_expression(left, right):
  """simple post-parsing creation method for subtract expressions"""
  return ApplyExpression("minus", [left, right])
def make_times_expression(left, right):
  """simple post-parsing creation method for multiply expressions"""
  return ApplyExpression("times", [left, right])
def make_divide_expression(left, right):
  """simple post-parsing creation method for divide expressions"""
  return ApplyExpression("divide", [left, right])


expr = Forward()

argument_list = Optional(create_separated_by(expr, ","))


# This way seems to work, provided in 'term' we have apply_expr
# before name_expr. We could also try something like:
# apply_expr = parcon.alpha_word + Optional("(" + argument_list + "))
# provided we can then interpret the result correctly, in particular
# we must make sure that "h()" is different from "h", the former being
# an application of the function h to no arguments and the latter being
# just a name expression referencing the 'h' variable.
variable_name = parcon.alpha_word
name_expr = Translate (parcon.alpha_word, make_name_expression)
apply_expr = Translate (parcon.alpha_word + "(" + argument_list + ")",
                        make_apply_expression)

term = ( parcon.number[create_number_expression] 
         | "(" + expr + ")" 
         | apply_expr
         | name_expr
       )
term = InfixExpr(term, [("*", make_times_expression), 
                        ("/", make_divide_expression)])
term = InfixExpr(term, [("+", make_plus_expression), 
                        ("-", make_minus_expression)])
expr.set(term(name="expr"))


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

def create_variable_dec(parse_result):
  """The parse action for the variable_definition parser"""
  return VariableDeclaration (parse_result[0], parse_result[1])

variable_definition = Translate(variable_name + "=" + expr + ";",
                                create_variable_dec)


class RateDefinition:
  """The representation of a rate definition"""
  def __init__(self, name, rate):
    self.name = name
    self.rate = rate

  def get_name(self):
    """Return the name of the rate being defined (or more rather the
       name of the associated reaction)"""
    return self.name

  def get_rate(self):
    """Return the rate expression part of the rate definition"""
    return self.rate

def create_rate_def(parse_result):
  """Post parsing method for rate definitions"""
  return RateDefinition(parse_result[0], parse_result[1])

rate_def_syntax = variable_name + "=" + "[" + expr + "]" + ";"
rate_def_parser = Translate(rate_def_syntax, create_rate_def)

class ComponentDefinition:
  """The representation of a component definition"""
  def __init__(self, name, behaviours):
    self.name = name
    self.behaviours = behaviours

  def get_name(self):
    """Return the name of the component being defined"""
    return self.name

  def get_behaviours(self):
    """Return the behaviours of this component"""
    return self.behaviours

  def show_definition(self):
    """Prints out the component definition in Bio-PEPA format"""
    behaviour_strings = [ x.show_behaviour() for x in self.behaviours ]
    return self.name + " = " + " + ".join(behaviour_strings) + " ;"

class Behaviour:
  """A class representing the behaviour of a component within a reaction.
     That is more intuitively part of a component's definition"""
  def __init__(self, reaction, operator):
    self.reaction_name = reaction
    self.operator = operator
    self.stoichiometry = 1

  def get_name(self):
    """Return the name of the reaction to which the behaviour is
       referring"""
    return self.reaction_name

  def set_stoiciometry(self, stoich):
    """Set the stoichiometry"""
    self.stoichiometry = stoich

  def show_behaviour(self):
    """Print out the behaviour as a string in Bio-PEPA format"""
    return self.reaction_name + " " + self.operator

behaviour_op = parcon.First (parcon.SignificantLiteral(">>"), 
                             parcon.SignificantLiteral("<<"),
                             parcon.SignificantLiteral("(+)"),
                             parcon.SignificantLiteral("(-)"),
                             parcon.SignificantLiteral("(.)"))
behaviour_op_list = parcon.OneOrMore(behaviour_op)


def make_behaviour(parse_result):
  """The parse action for the behaviour parser"""
  print ("--- parse result ---")
  print (parse_result)
  reaction_name = parse_result[0]
  stoich = parse_result[1]
  operator = parse_result[2]
  behaviour = Behaviour(reaction_name, operator)
  behaviour.set_stoiciometry(stoich)
  return behaviour

stoich_list = create_separated_by(parcon.number, ",")
explicit_stoich_syntax = "(" + variable_name + "," + stoich_list + ")"
implicit_stoich_syntax = Translate(variable_name, lambda x : (x,[1]))
name_and_stoich = parcon.First(implicit_stoich_syntax,
                               explicit_stoich_syntax)
                               

behaviour_syntax = Translate(name_and_stoich + behaviour_op_list,
                             make_behaviour)


behaviour_list_parser = create_separated_by(behaviour_syntax, "+")

def make_component_def(parse_result):
  """The parse action for the component definition parser"""
  return ComponentDefinition(parse_result[0], parse_result[1])
component_definition_syntax = (variable_name + "=" + 
                               behaviour_list_parser + ";")
component_definition = Translate(component_definition_syntax,
                                 make_component_def)

any_definition = parcon.First(variable_definition,
                              rate_def_parser,
                              component_definition
                              )
definition_list = parcon.OneOrMore(any_definition)

class ComponentPopulation:
  """A class to hold the representation of a system_equation component.
     This is essentially a name and initial population
  """
  def __init__(self, name, population_expr):
    self.name = name
    self.population_expr = population_expr
    self.compartment = "default_compartment"

  def get_name(self):
    return self.name
  def get_compartment(self):
    return self.compartment

def create_component_population(parse_result):
  """A post-parsing function to create a component population"""
  return ComponentPopulation(parse_result[0], parse_result[1])

component_population = Translate(variable_name + "[" + expr + "]",
                                 create_component_population)
system_equation_parser = create_separated_by(component_population, "<*>")

class BioPEPAModel:
  """A simple class to hold the representation of a Bio-PEPA model"""
  def __init__(self, definitions, system_equation):
    self.definitions = definitions
    self.system_equation = system_equation

def make_model(parse_result):
  """Simple post-parsing creation function for the model parser"""
  return BioPEPAModel(parse_result[0], parse_result[1])
model_syntax = definition_list + system_equation_parser + parcon.End()
model_parser = Translate(model_syntax, make_model)


def parse_model(model_source):
  """A significant entry point to this module. Takes in the string
     of the source of a Bio-PEPA model and return the parse result
     of parsing the model
  """
  return model_parser.parse_string(model_source)

class ReactionParticipant:
  """A simple class to represent a reaction participant"""
  def __init__(self, name):
    self.name = name
    self.stoichiometry = 1

  def get_name(self):
    """Return the name of the reaction particpant"""
    return self.name

  def set_stoichiometry(self, stoichiometry):
    """Set the stoichiometry of the reaction participant"""
    self.stoichiometry = stoichiometry

class Reaction:
  """ A class representing a Bio-PEPA reaction"""
  def __init__(self, name):
    self.name = name
    # The user should not be forced into providing a rate
    # they may only wish to perform invariant analysis for example.
    self.rate = None
    self.reactants = []
    self.products = []
    self.modifiers = []

  def define_rate(self, rate):
    """The default rate is undefined since we may be attempting to
       perform a rateless analysis over the model and we needn't then
       insist that the user writes down rates. This method sets the
       rate to the given expression
    """
    self.rate = rate

  def get_mass_action_participants(self):
    """Return the left hand side participants which contribute to the
       rate in a mass action rate method. For example A + B --> C 
       would return A and B. Essentially this is so that fMA(r) could
       be translated to r * A * B.
    """
    # I think we return all the modifiers but perhaps not the inhibitors?
    return self.reactants + self.modifiers

  def add_behaviour(self, comp_name, behaviour):
    """Add the given behaviour to the reaction, essentially this means
       we are adding a component as a participant (reactant, product
       or modifier) to the reaction.
    """
    stoichimetries = behaviour.stoichiometry
    operators = behaviour.operator
    if len(stoichimetries) != len(operators):
      print ("Must have the same number of operators as stoichimetries")
      sys.exit(1)

    # Not sure if this is really the best way to do this, this means
    # we will be adding multiple behaviours for the same component.
    # eg, we'll have A + B -> B + B, where as arguably we should retain
    # the knowledge, perhaps we should just have the operator <> or ><
    for stoich, operator in zip(stoichimetries, operators):
      participant = ReactionParticipant(comp_name)
      participant.set_stoichiometry(stoich)
      if operator == ">>" :
        self.products.append(participant)
      elif operator == "<<": 
        self.reactants.append(participant)
      elif operator == "(+)":
        self.modifiers.append(participant)
      elif operator == "(-)":
        self.modifiers.append(participant)
      elif operator == "(.)":
        self.modifiers.append(participant)
      else:
        print ("Unrecognised behaviour operator: " + operator)
        sys.exit(1)

  def create_element(self, document):
    """Create an xml element representing this reaction"""
    reaction_element = document.createElement("reaction")
    reaction_element.setAttribute("id", self.name)
    reaction_element.setAttribute("reversible", "false")
    reaction_element.setAttribute("fast", "false")
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
        spec_ref.setAttribute("stoichiometry", str(reactant.stoichiometry))
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
        spec_ref.setAttribute("stoichiometry", str(product.stoichiometry))
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
    if self.rate:
      kinetic_law = document.createElement("kineticLaw")
      reaction_element.appendChild(kinetic_law)
      math_element = document.createElement("math")
      kinetic_law.appendChild(math_element)
      mathxmlns = "http://www.w3.org/1998/Math/MathML"
      math_element.setAttribute("xmlns", mathxmlns)
      simplified_rate = self.rate.remove_rate_law_sugar(self)
      expr_element = simplified_rate.create_sbml_element(document)
      math_element.appendChild(expr_element)
    return reaction_element

def create_compartment_elements(document, model_element):
  """This is likely to change signature when we allow for user defined
     compartments but for now we have only a single default compartment
     to add to the model"""
  list_of_compartments = document.createElement("listOfCompartments")
  model_element.appendChild(list_of_compartments)

  default_compartment = document.createElement("compartment")
  default_compartment.setAttribute("id", "default_compartment")
  default_compartment.setAttribute("constant", "true")
  list_of_compartments.appendChild(default_compartment)

def create_species_elements(document, model_element, component_defs):
  """Creates the listOfSpecies element and all of the associated
     species elements under it for a set of component definitions"""
  list_of_species = document.createElement("listOfSpecies")
  model_element.appendChild(list_of_species)
  for comp_def in component_defs:
    species_element = document.createElement("species")
    name = comp_def.get_name()
    species_element.setAttribute("id", name)
    species_element.setAttribute("name", name)
    species_element.setAttribute("compartment", comp_def.get_compartment())
    species_element.setAttribute("hasOnlySubstanceUnits", "true")
    species_element.setAttribute("constant", "false")
    species_element.setAttribute("boundaryCondition", "false")
    list_of_species.appendChild(species_element) 

def convert_variable_declarations(document, top_element, var_decs):
  """Convert a set of variable declarations into SBML elements and
     add them to an sbml document
  """
  if var_decs:
    list_of_params = document.createElement("listOfParameters")
    top_element.appendChild(list_of_params)
 
    init_assigns = document.createElement("listOfInitialAssignments")
    top_element.appendChild(init_assigns)
    for var_dec in var_decs:
      param_element = var_dec.create_parameter_element(document)
      list_of_params.appendChild(param_element)

      init_assign = var_dec.create_initial_assignment(document)
      init_assigns.appendChild(init_assign)


def build_reaction_dictionary(components, rate_definitions):
  """From a list of components build a reaction dictionary which maps
     reaction names to reaction representations"""
  reaction_dictionary = dict()
  for rate_def in rate_definitions:
    name = rate_def.get_name()
    reaction = Reaction(name)
    reaction.define_rate(rate_def.get_rate())
    reaction_dictionary[name] = reaction


  for component_def in components:
    for behaviour in component_def.get_behaviours():
      b_name = behaviour.get_name() 
      if b_name in reaction_dictionary:
        reaction = reaction_dictionary[b_name]
      else:
        reaction = Reaction(b_name)
        reaction_dictionary[b_name] = reaction
      reaction.add_behaviour(component_def.get_name(), behaviour)
  return reaction_dictionary



def process_file(filename):
  """A simple method to process a Bio-PEPA model file and print
     out the parse result, mostly for debugging purposes"""
  # I'm not exactly sure how to do this, there is no parseFile
  # in parcon, I think we just need to open the file and pass in
  # the file handle but I haven't tried that yet. 
  model_file = open(filename, "r")
  parse_result = model_parser.parse_string(model_file.read())
  model_file.close()

  # print (parse_result.asXML())

  var_decs = [ x for x in parse_result.definitions
                     if isinstance(x, VariableDeclaration) ]
  rate_defs = [ x for x in parse_result.definitions
                      if isinstance(x, RateDefinition) ]
  components = [ x for x in parse_result.definitions
                       if isinstance(x, ComponentDefinition) ]

  reaction_dictionary = build_reaction_dictionary(components, rate_defs)

  xml_implementation = minidom.getDOMImplementation()
  document = xml_implementation.createDocument(None, "sbml", None)
  top_element = document.documentElement 
  xmlns = "http://www.sbml.org/sbml/level3/version1/core" 
  top_element.setAttribute("xmlns", xmlns)
  top_element.setAttribute("level", "3")
  top_element.setAttribute("version", "1")

  model_element = document.createElement("model")
  top_element.appendChild(model_element)

  create_compartment_elements(document, model_element)

  create_species_elements(document, model_element,
                          parse_result.system_equation)

  convert_variable_declarations(document, model_element, var_decs)

  list_of_reactions = document.createElement("listOfReactions")
  model_element.appendChild(list_of_reactions)

  for reaction in reaction_dictionary.values():
    reaction_element = reaction.create_element(document)
    list_of_reactions.appendChild(reaction_element)
     

  # top_element.writexml(sys.stdout, indent="  ")
  # print ("")
  formatted = document.toprettyxml(indent="  ", encoding="UTF-8")
  print (formatted)

def main():
  """A simple main function to parse in the arguments as
     Bio-PEPA model files"""
  description = "Parse a Bio-PEPA model file(s)"
  parser = argparse.ArgumentParser(add_help=True,
                                   description=description)
   # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="A Bio-PEPA model file to parse")
  arguments = parser.parse_args()
  for filename in arguments.filenames:
    process_file(filename)


if __name__ == '__main__':
  main()


