import argparse
import parcon
from parcon import Forward, InfixExpr, Translate, Optional, ZeroOrMore


# A simply utility for creating parsers which accept a list of
# somethings, separated by something elses.
# Note here that the separator_parser should return None
def create_separated_by(element_parser, separator_parser):
  def create_list_result(parse_result):
    first_item = parse_result[0]
    rest_items = parse_result[1]
    rest_items.insert(0, first_item)
    return rest_items
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


def make_apply_expression(parse_result):
  return ApplyExpression(parse_result[0], parse_result[1])

def make_plus_expression(left, right):
  return ApplyExpression("plus", [left, right])
def make_minus_expression(left, right):
  return ApplyExpression("minus", [left, right])
def make_times_expression(left, right):
  return ApplyExpression("times", [left, right])
def make_divide_expression(left, right):
  return ApplyExpression("divide", [left, right])


expr = Forward()

expr_list = Optional(expr + Optional(ZeroOrMore(", " + expr)))

def make_argument_list(exprs):
  if not exprs:
    return []
  else:
    return [ exprs[0] ] + exprs[1]

argument_list = Translate(expr_list, make_argument_list)


# This way seems to work, provided in 'term' we have apply_expr
# before name_expr. We could also try it something like:
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
expr << term(name="expr")


class VariableDeclaration:
  """A class to represent the AST of a variable declaration in Bio-PEPA"""
  def __init__(self, variable, expression):
    self.variable = variable
    self.expression = expression

  def convert_to_sbml(self):
    """Convert the variable declaration into an SBML parameter definition.
       If the expression is a simple number then this is added as a
       'value' attribute to the parameter element"""
    value_att = ""
    if isinstance(self.expression, NumExpression):
      value_att = "value=\"" + self.expression.show_number() + "\""
    return "<parameter id=\"" + self.variable + "\" " + value_att + "/>"

  def show_initial_assignment(self):
    """Print out the SBML initial assigment for the variable declaration.
       This does so blindly whether or not an initial assignment is
       required"""
    # Here we are giving each variable an initial assignment
    # We could instead decide that it doesn't need one because it
    # has a value attribute associated with its parameter def.
    result = "<initialAssignment symbol=\"" + self.variable + "\">"
    result += "<math xmlns=\"http://www.w3.org/1998/Math/MathML\"\n"
    result += "      xmlns:sbml=\"http://www.sbml.org/sbml/level3/"
    result += "version1/core\">"

    result += self.expression.convert_to_sbml()

    result += "\n</math>\n</initialAssignment>"
    
    return result

def create_variable_dec(parse_result):
  """The parse action for the variable_definition parser"""
  return VariableDeclaration (parse_result[0], parse_result[1])

variable_definition = Translate(variable_name + "=" + expr + ";",
                                create_variable_dec)

class ComponentDefinition:
  """The representation of a component definition"""
  def __init__(self, name, behaviours):
    self.name = name
    self.behaviours = behaviours

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

  def show_behaviour(self):
    """Print out the behaviour as a string in Bio-PEPA format"""
    return self.reaction_name + " " + self.operator

behaviour_op = parcon.First (parcon.SignificantLiteral(">>"), 
                             parcon.SignificantLiteral("<<"),
                             parcon.SignificantLiteral("(+)"),
                             parcon.SignificantLiteral("(-)"),
                             parcon.SignificantLiteral("(.)"))
def make_behaviour(parse_result):
  """The parse action for the behaviour parser"""
  return Behaviour(parse_result[0], parse_result[1])
behaviour = Translate(variable_name + behaviour_op,
                      make_behaviour)


behaviour_list_parser = create_separated_by(behaviour, "+")

def make_component_def(parse_result):
  """The parse action for the component definition parser"""
  return ComponentDefinition(parse_result[0], parse_result[1])
component_definition_syntax = (variable_name + "=" + 
                               behaviour_list_parser + ";")
component_definition = Translate(component_definition_syntax,
                                 make_component_def)

any_definition = parcon.First(variable_definition,
                              component_definition
                              )
definition_list = parcon.OneOrMore(any_definition)

class ComponentPopulation:
  def __init__(self, name, expr):
    self.name = name
    self.population_expr = expr
def create_component_population(parse_result):
  return ComponentPopulation(parse_result[0], parse_result[1])

component_population = Translate(variable_name + "[" + expr + "]",
                                 create_component_population)
system_equation = create_separated_by(component_population, "<*>")

# This is clearly pretty temporary and just returning the
# definition list.
def make_model(parse_result):
  return parse_result[0]
model_syntax = definition_list + system_equation + parcon.End()
model_parser = Translate(model_syntax, make_model)


def parse_model(model_source):
  """A significant entry point to this module. Takes in the string
     of the source of a Bio-PEPA model and return the parse result
     of parsing the model
  """
  return model_parser.parse_string(model_source)



def process_file(filename):
  """A simple method to process a Bio-PEPA model file and print
     out the parse result, mostly for debugging purposes"""
  # I'm not exactly sure how to do this, there is no parseFile
  # in parcon, I think we just need to open the file and pass in
  # the file handle but I haven't tried that yet. 
  parse_result = model_parser.parseFile(filename)

  # print (parse_result.asXML())

  for mvd in parse_result:
    if isinstance(mvd, VariableDeclaration):
      print(mvd.convert_to_sbml())
  for mvd in parse_result:
    if isinstance(mvd, VariableDeclaration):
      print(mvd.show_initial_assignment())

  components = [ x for x in parse_result 
                       if isinstance(x, ComponentDefinition) ]

  for component_def in components:
    print(component_def.show_definition())
  # When you have something for a whole biopepafile you can use
  # parser.parseFile(source_file)

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


