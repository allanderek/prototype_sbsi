import parcon
from parcon import Forward, InfixExpr, Translate, Optional, ZeroOrMore

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
    operator_el = document.createElement(rename_operator(self.name))
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

