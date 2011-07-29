"""A module for parsing in Bio-PEPA files together with the necessary
   data type definitions to hold the AST of a Bio-PEPA model.
"""
import argparse
import pyparsing
from pyparsing import ( Literal, Word, OneOrMore, Optional, 
                        Combine, nums, Or
                      )

def update_forwarded_parser(forward, parser):
  """A simple utility function to help with forwarded parsers.
     These are created with pyparsing.Forward and are used when
     you need to create a recursive parser definition.
     expr = pyparsing.Forward()
     infix_expression = Or([ expr + operator + expr, number])
     Now we need to update 'expr' to be 'infix_expression. The usual
     way to do this is,
     expr << infix_expression
     but pylint complains about the statement not assigning to anything
     so insteady we use:
     expr.__lshift__(infix_expression)
     but this looks a bit ugly so this utility function helps reduce
     the ugliness a little
  """
  forward.__lshift__(parser)  

capital_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
small_letters = capital_letters.lower()
digits = "0123456789"
all_letters = capital_letters + small_letters

variable_name = Word(all_letters, all_letters + digits)
integer_number = Word(digits)
sci_notation_start = pyparsing.CaselessLiteral( "E" )
# Must be quite careful about this, it seems to lenient when we allow
# also binary operations, for example you can have weird things like:
# 1 + -2. The operatorPrecendence allows for sign operators to be used
# so we don't need the sign here.
float_number = Combine( Word( nums ) + 
                        Optional( "." + Optional( Word( nums ) ) ) +
                        Optional( sci_notation_start + 
                                  Word( "+-"+nums, nums ) ) )
float_number.setParseAction(float)

number_expression = float_number
def make_num_expression(numbers):
  """The parse action for literal number expression parser"""
  return NumExpression(numbers[0])
number_expression.setParseAction(make_num_expression)

# I could just use variable_name here, but we do not want to
# set the parse action of the variable name parser hence we use 'copy'.
# name_expression = variable_name.copy()
# def make_name_expression(names):
#   """The parse action for a name expression parser"""
#   return VarExpression(names[0])
# name_expression.setParseAction(make_name_expression)

expression_parser = pyparsing.Forward()


expression_list = pyparsing.delimitedList(expression_parser,
                                          delim=",",
                                          combine=False)
argument_list = Optional(expression_list)
apply_expression = variable_name + Optional("(" + argument_list + ")")
def create_apply_expression(tokens):
  """The parse action for an apply expression, specifically an
     expression which is a call of some named function to a list of
     arguments."""
  name = tokens[0]

  # We might not have a apply expression here at all, if there were
  # no following arguments (including the empty argument list defined
  # as ()) then it is a name expression. Note that f is really different
  # from f(). One is a name expression the other is the application of
  # a function to zero arguments.
  if len(tokens) < 2:
    return VarExpression(name)
  # tokens[1] is opening bracket, at len-1 is closing bracket
  args = tokens[2:len(tokens) - 1]
  return ApplyExpression(name, args)
apply_expression.setParseAction(create_apply_expression)


expop = Literal('^')
signop = pyparsing.oneOf('+ -')
multop = pyparsing.oneOf('* /')
plusop = pyparsing.oneOf('+ -')

def create_operator_expression(tokens):
  """The parse action for an operator expression to return the AST
     representation given the list of tokens parsed for the expression"""
  if len(tokens) > 3:
    print ("No idea how to deal with greater than 3 token expressions")
  elif len(tokens) == 1 and isinstance(tokens[0], pyparsing.ParseResults):
    return create_operator_expression(tokens[0])
  elif len(tokens) == 1:
    return tokens[0]
  elif len(tokens) == 2:
    # Assume it is a unary operator expression
    return ApplyExpression(tokens[0], tokens[1:])
  elif len(tokens) == 3:
    # assume it is a binary operator expression
    return ApplyExpression(tokens[1], [tokens[0], tokens[2]])
  else:
    print ("No idea how to and token expressions with no tokens")


operand = Or([number_expression, apply_expression])
# Careful the number of operators here has a significant impact
# on how quickly nested function applications are parsed, I'm not quite
# sure why. Definitely don't need factorial. Seems like 3 is the magical
# number.
operator_expression = pyparsing.operatorPrecedence( operand,
    [ # ("!", 1, pyparsing.opAssoc.LEFT, create_operator_expression),
      (signop, 1, pyparsing.opAssoc.RIGHT, create_operator_expression),
      # ("^", 2, pyparsing.opAssoc.RIGHT, create_operator_expression),
      (multop, 2, pyparsing.opAssoc.LEFT, create_operator_expression),
      (plusop, 2, pyparsing.opAssoc.LEFT, create_operator_expression),
    ])

operator_expression.setParseAction(create_operator_expression)

# Now that we have defined our expression parsers we can update the
# definition of an expression parser
update_forwarded_parser(expression_parser, operator_expression)
                        # Or([apply_expression, operator_expression]))

variable_definition = variable_name + "=" + expression_parser + ";"

class Expression:
  """The base class for all classes which represent the AST of some
     kind of expression"""
  def __init__(self):
    pass

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
  def show_number(self):
    """Display the underlying number of the numerical expression"""
    return str(self.number)

class VarExpression(Expression):
  """A class to represent the AST of a variable (name) expression"""
  def __init__(self, name):
    Expression.__init__(self) 
    self.name = name

  def convert_to_sbml(self):
    """Convert the variable(name) expression to SBML code for
       the expression"""
    return "<ci>" + self.name + "</ci>"



def rename_operator(name):
  """A simple utility operator to convert an operator in to the named
     function to perform the same operation as an apply expression"""
  if name == "+":
    return "plus"
  if name == "-":
    return "minus"
  if name == "*":
    return "times"
  if name == "/":
    return "divide"
  return name

class ApplyExpression(Expression):
  """A class to represent the AST of an apply expression, applying a
     named function to a list of argument expressions"""
  def __init__(self, name, args):
    Expression.__init__(self)
    self.name = name
    self.args = args

  def convert_to_sbml(self):
    """return a string representing the sbml of an math apply expression"""
    result = "<apply>\n"
    result += "  <" + rename_operator(self.name) + "/>\n"
    for argument in self.args:
      result += "  " + argument.convert_to_sbml() + "\n"
    result += "</apply>"
    return result

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

def create_variable_dec(tokens):
  """The parse action for the variable_definition parser"""
  return VariableDeclaration (tokens[0], tokens[2])
variable_definition.setParseAction(create_variable_dec)

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

behaviour_op = pyparsing.oneOf ([ ">>", "<<", "(+)", "(-)", "(.)" ])
behaviour = variable_name + behaviour_op
def make_behaviour(tokens):
  """The parse action for the behaviour parser"""
  return Behaviour(tokens[0], tokens[1])
behaviour.setParseAction(make_behaviour)

behaviour_list_parser = pyparsing.delimitedList(behaviour,
                                                delim="+",
                                                combine=False)

component_definition = variable_name + "=" + behaviour_list_parser + ";"
def make_component_def(tokens):
  """The parse action for the component definition parser"""
  # We want the list from the 3rd token to everything but the last.
  # The last token is the semi-colon. The first is the component name
  # whilst the 2nd token (token[1]) is the equals sign.
  my_behaviours = tokens[2:len(tokens) - 1]
  return ComponentDefinition(tokens[0], my_behaviours)
component_definition.setParseAction(make_component_def)


any_definition = Or([variable_definition,
                     component_definition,
                    ])
definition_list = OneOrMore(any_definition)

# There should be a way for us to say that the parse result is
# simply that of 'definition_list'
model_parser = definition_list + pyparsing.StringEnd()


def parse_model(model_source):
  """A significant entry point to this module. Takes in the string
     of the source of a Bio-PEPA model and return the parse result
     of parsing the model
  """
  return model_parser.parseString(model_source)



def process_file(filename):
  """A simple method to process a Bio-PEPA model file and print
     out the parse result, mostly for debugging purposes"""
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
