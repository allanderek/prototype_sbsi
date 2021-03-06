"""
A module that implements a parser for the Bio-PEPA language
"""
import sys
import parcon
from parcon import Forward, InfixExpr, Translate, Optional, ZeroOrMore

import sbml_ast
import expressions

# A simply utility for creating parsers which accept a list of
# somethings, separated by something elses.
# Note here that the separator_parser should return None
# Note also that I believe the latest version of parcon has this
# (or similar functionality) included.
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


def create_number_expression(number_str):
  """Simple utility to create a number expression from the parse result
     of parsing a simple number"""
  return expressions.NumExpression(float(number_str))
 
def make_name_expression(parse_result):
  """Simple post parsing creation method for NameExpression"""
  return expressions.NameExpression(parse_result)

   
def make_apply_expression(parse_result):
  """simple post-parsing creation method for apply expressions"""
  if len(parse_result) == 1:
    # We assume then that there are no arguments, and that it was
    # an application expression such as: f()
    return expressions.ApplyExpression(parse_result[0], [])
  # Otherwise the second part of the parse result should be a list
  # of arguments
  return expressions.ApplyExpression(parse_result[0], parse_result[1])

def make_plus_expression(left, right):
  """simple post-parsing creation method for add expressions"""
  return expressions.ApplyExpression("plus", [left, right])
def make_minus_expression(left, right):
  """simple post-parsing creation method for subtract expressions"""
  return expressions.ApplyExpression("minus", [left, right])
def make_times_expression(left, right):
  """simple post-parsing creation method for multiply expressions"""
  return expressions.ApplyExpression("times", [left, right])
def make_divide_expression(left, right):
  """simple post-parsing creation method for divide expressions"""
  return expressions.ApplyExpression("divide", [left, right])
def make_power_expression(left, right):
  """simple post-parsing creation method for power expressions"""
  return expressions.ApplyExpression("power", [left, right])


expr = Forward()

argument_list = Optional(create_separated_by(expr, ","))


# This way seems to work, provided in 'term' we have apply_expr
# before name_expr. We could also try something like:
# apply_expr = parcon.alpha_word + Optional("(" + argument_list + "))
# provided we can then interpret the result correctly, in particular
# we must make sure that "h()" is different from "h", the former being
# an application of the function h to no arguments and the latter being
# just a name expression referencing the 'h' variable.
#
# Note also that parcon.alphanum_word may be incorrect here, it should
# be alphchar followed by alphanum_word, or something like that.
variable_name = parcon.Word(parcon.alphanum_chars + "_", 
                            init_chars=parcon.alpha_chars)
location_suffix = parcon.SignificantLiteral("@") + variable_name
located_name = (variable_name + Optional(location_suffix))["".join]
name_expr = Translate (located_name, make_name_expression)
apply_expr = Translate (variable_name + "(" + argument_list + ")",
                        make_apply_expression)

my_number_syntax = parcon.Exact(Optional(parcon.CharIn("+-")) + 
                                Optional(parcon.ZeroOrMore(parcon.digit) + 
                                         parcon.SignificantLiteral(".")) +
                                parcon.OneOrMore(parcon.digit))
my_number_parser = my_number_syntax[parcon.flatten]["".join](name="number")
signed_number = my_number_parser
exponent_end = (parcon.CharIn("Ee") + signed_number)["".join]
scientific_number = (signed_number + Optional(exponent_end))["".join]

term = ( scientific_number[create_number_expression] 
         | "(" + expr + ")" 
         | apply_expr
         | name_expr
       )

term = InfixExpr(term, [("^", make_power_expression)])
term = InfixExpr(term, [("*", make_times_expression), 
                        ("/", make_divide_expression)])
term = InfixExpr(term, [("+", make_plus_expression), 
                        ("-", make_minus_expression)])
expr.set(term(name="expr"))


class VariableDeclaration(sbml_ast.VariableDeclaration):
  """A class to represent Bio-PEPA variable declarations"""
  def __init__(self, variable, expression):
    super(VariableDeclaration, self).__init__(variable, expression)

  def format_declaration(self):
    """A function to return a Bio-PEPA formatted string representation
       of this variable declaration
    """
    return self.variable + " = " + self.expression.show_expr() + " ;"

def create_variable_dec(parse_result):
  """The parse action for the variable_definition parser"""
  return VariableDeclaration (parse_result[0], parse_result[1])

variable_definition = Translate(variable_name + "=" + expr + ";",
                                create_variable_dec)


class LocationDefinition(object):
  """The representation of a location definition"""
  def __init__(self, name, size_expr):
    self.name = name
    self.size_expr = size_expr

def create_location_definition(parse_result):
  """Post parsing method for location definitions"""
  return LocationDefinition(parse_result[0], parse_result[1])

# location M : size = 1;
location_name_syntax = variable_name
location_definition_syntax = ( "location" + 
                               location_name_syntax + 
                               ":" +
                               "size" +
                               "=" +
                               expr +
                               ";" )
location_definition_parser = Translate(location_definition_syntax,
                                       create_location_definition)


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

  def format_rate_def(self):
    """Return a string of a Bio-PEPA formatted rate definition"""
    return self.name + " = [ " + self.rate.show_expr() + " ] ;"

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
    behaviour_strings = [ x.show_behaviour(self.name) 
                            for x in self.behaviours ]
    return self.name + " = " + " + ".join(behaviour_strings) + " ;"

class Behaviour:
  """A class representing the behaviour of a component within a reaction.
     That is more intuitively part of a component's definition"""
  def __init__(self, reaction, operator):
    self.reaction_name = reaction
    self.operator = operator
    self.stoichiometry = ['1']
    self.location = None

  def get_name(self):
    """Return the name of the reaction to which the behaviour is
       referring"""
    return self.reaction_name

  def set_stoiciometry(self, stoich):
    """Set the stoichiometry"""
    self.stoichiometry = stoich

  def show_behaviour(self, component_name):
    """Print out the behaviour as a string in Bio-PEPA format"""
    if self.stoichiometry == ['1']:
      result = self.reaction_name + " " + self.operator[0]
    else:
      stoichs = ", ".join([str(s) for s in self.stoichiometry])
      opers = " ".join(self.operator)
      result =  ("(" + self.reaction_name + ", " + 
             stoichs + ") " + opers)

    if self.location != None:
      result += " " + component_name + "@" + self.location

    return result

  def set_location (self, location):
    """Sets the location of this behaviour"""
    self.location = location

behaviour_op = parcon.First (parcon.SignificantLiteral(">>"), 
                             parcon.SignificantLiteral("<<"),
                             parcon.SignificantLiteral("(+)"),
                             parcon.SignificantLiteral("(-)"),
                             parcon.SignificantLiteral("(.)"))
behaviour_op_list = parcon.OneOrMore(behaviour_op)


def make_behaviour(parse_result):
  """The parse action for the behaviour parser"""
  reaction_name = parse_result[0]
  stoich = parse_result[1]
  operator = parse_result[2]
  behaviour = Behaviour(reaction_name, operator)
  behaviour.set_stoiciometry(stoich)
  if len(parse_result) > 3:
    # We should check somehow that the component is the correct name
    component_placement = parse_result[3]
    behaviour.set_location(component_placement.location)
  return behaviour

stoich_list = create_separated_by(parcon.number, ",")
explicit_stoich_syntax = "(" + variable_name + "," + stoich_list + ")"
implicit_stoich_syntax = Translate(variable_name, lambda x : (x, ['1']))
name_and_stoich = parcon.First(implicit_stoich_syntax,
                               explicit_stoich_syntax)
                               
component_name_syntax = variable_name
component_name_placement = located_name

class ComponentLocation(object):
  """A class representing a component as well as a location, which
     might be None. Intended to represent the parsed "A@M" syntax
  """
  def __init__(self, name, location):
    self.name = name
    self.location = location

def create_component_with_location(parse_result):
  """Post-parsing method for a component with a location,
     eg "A@M".
  """
  names = parse_result.split("@") 
  component = names[0]
  if len(names) == 1:
    location = None
  else:
    # I must admit I'd be a little happier to check that it's not
    # longer than two long, in which case we had multiple @ symbols
    # and I wouldn't be sure about what to do.
    location = names[1]
  return ComponentLocation(component, location)

component_with_location = Translate(component_name_placement,
                                    create_component_with_location)
optional_component_placement = Optional (component_with_location)
                             

behaviour_syntax = Translate(name_and_stoich +
                             behaviour_op_list + 
                             optional_component_placement,
                             make_behaviour)


behaviour_list_parser = create_separated_by(behaviour_syntax, "+")

def make_component_def(parse_result):
  """The parse action for the component definition parser"""
  # print ("component def parse result")
  # print (parse_result)
  return ComponentDefinition(parse_result[0], parse_result[1])
component_definition_syntax = (variable_name + "=" + 
                               behaviour_list_parser + ";")
component_definition = Translate(component_definition_syntax,
                                 make_component_def)

any_definition = parcon.First(variable_definition,
                              rate_def_parser,
                              component_definition,
                              location_definition_parser
                              )
definition_list = parcon.OneOrMore(any_definition)

class ComponentPopulation:
  """A class to hold the representation of a system_equation component.
     This is essentially a name and initial population
  """
  def __init__(self, name, initial_expression):
    self.name = name
    # The 'create_sbml' module requires this to have an 'initial_amount'
    # field. In theory we could populate this if the expression is a
    # simple number.
    self.initial_amount = None
    self.initial_expression = initial_expression
    self.location = "default location"

  def format(self):
    """Return a string of a Bio-PEPA formatted component population"""
    result = self.name
    if self.location != None and self.location != "default location":
      result += "@" + self.location
    result += "[" + self.initial_expression.show_expr() + "]"
    return result
  

def create_component_population(parse_result):
  """A post-parsing function to create a component population"""
  component_location = parse_result[0]
  component_name = component_location.name
  population = parse_result[1]
  result = ComponentPopulation(component_name, population)
  result.location = component_location.location
  return result

component_population = Translate(component_with_location + 
                                 "[" + expr + "]",
                                 create_component_population)
system_equation_parser = create_separated_by(component_population, "<*>")

class BioPEPAModel:
  """A simple class to hold the representation of a Bio-PEPA model"""
  def __init__(self):
    self.var_decs = None
    self.rate_defs = None
    self.component_defs = None
    self.system_equation = None

  def output_model(self, outfile):
    """Output the model to the given outfile"""

    if self.var_decs != None and self.var_decs:
      outfile.write("// Variable declarations\n")
      for var_dec in self.var_decs:
        outfile.write(var_dec.format_declaration())
        outfile.write("\n")
      outfile.write("\n")
 
    # Admittedly it would be a strange Bio-PEPA outfile indeed which
    # contained no rate definitions and the same goes for component
    # definitions, but there is nothing wrong with a bit of defensive
    # programming, what if someone is extracting a module or something?
    if self.rate_defs != None and self.rate_defs:
      outfile.write("// Rate definitions\n")
      for rate_def in self.rate_defs:
        outfile.write(rate_def.format_rate_def())
        outfile.write("\n")
      outfile.write("\n")

    if self.component_defs != None and self.component_defs:
      outfile.write("// Component defintions\n")
      for component_def in self.component_defs:
        outfile.write(component_def.show_definition())
        outfile.write("\n")
      outfile.write("\n")

    comp_pops = [ comp_pop.format() for comp_pop in self.system_equation ]
    system_equation = " <*>\n".join(comp_pops)

    outfile.write ("// The system equation\n")
    outfile.write (system_equation)
    outfile.write ("\n")

    
    

def make_model(parse_result):
  """Simple post-parsing creation function for the model parser"""
  biopepa_model = BioPEPAModel()
  definitions = parse_result[0]

  var_decs = [ x for x in definitions
                     if isinstance(x, VariableDeclaration) ]
  rate_defs = [ x for x in definitions
                      if isinstance(x, RateDefinition) ]
  components = [ x for x in definitions
                       if isinstance(x, ComponentDefinition) ]
  biopepa_model.var_decs = var_decs
  biopepa_model.rate_defs = rate_defs
  biopepa_model.component_defs = components

  biopepa_model.system_equation = parse_result[1] 
  return biopepa_model

model_syntax = definition_list + system_equation_parser + parcon.End()
model_parser = Translate(model_syntax, make_model)


def parse_model(model_source):
  """A significant entry point to this module. Takes in the string
     of the source of a Bio-PEPA model and return the parse result
     of parsing the model
  """
  # This code is nearly identical to code in the facile parser, we
  # should factor out the ability to easily add a line comment
  # facility to any parser.
  newline_syntax = parcon.CharIn("\n\r")
  not_new_line = parcon.Except(parcon.AnyChar(), newline_syntax)
  line_comment = "//" + parcon.ZeroOrMore(not_new_line)
  either_or = parcon.First(parcon.Whitespace(), line_comment) 
  comments_and_whitespace = parcon.OneOrMore(either_or)
  return model_parser.parse_string(model_source,
                                   whitespace=comments_and_whitespace)

def parse_model_file(model_file):
  """Parse an entire Bio-PEPA model file"""
  # parse_result = model_parser.parse_string(model_file.read())
  parse_result = parse_model(model_file.read())
  return parse_result

def parse_model_file_exit_on_error(model_file):
  """Calls the above parse_model_file and if there is a parse error
     prints the error report to the screen and then exits. Generally we
     expect callers of parse_model_file will do something more
     sophisticated but this works well for a quick script.
  """
  try:
    parse_result = parse_model_file(model_file)
  except parcon.ParseException as parse_except:
    print (parse_except)
    sys.exit(1)
  return parse_result
 
