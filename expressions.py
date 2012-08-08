""" A module implementing expressions, these are generally used as part of
    an SBML model, for example the rate expressions.
"""
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
       We raise the exception 'KeyError', if the value cannot be derived,
       this will generally be because a name used in the expression is not
       defined in the given environment (or there is no given environment).
    """
    # pylint: disable=W0613
    # pylint: disable=R0201
    raise ValueError("Virtual method 'get_value' called")

  def get_value_none(self, environment=None):
    """The same as 'get_value' except that it returns None in the case
       that we cannot evaluate the expression (because the environment
       either is not given or does not define some name used in the
       expression)
    """
    try:
      return self.get_value(environment=environment)
    except KeyError:
      return None

  def reduce(self, environment=None):
    """Similar to 'get_value' except that we always return an expression,
       and in the case that the expression is not wholly reducible to a
       number, it may be reducible to a more simpler expression, for example
       R * (factor ^ 2)
       if we are given 'factor' as a constant, let's say mapped to '2', then
       we can return the reduced expression
       R * 4
       The idea is that if the expression is something like a rate expression
       which must be re-evaluated many times, then we can save time by 
       partially evaluating it.
       However if the expression cannot be reduced then we simply return
       the original expression.
    """
    # pylint: disable=W0613
    # pylint: disable=R0201
    return self


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

  def visit (self, visitor):
    """Implements the visit method allowing ExpressionVisitors to work"""
    visitor.visit_NumExpression(self)

  def convert_to_sbml(self):
    """Convert the number expression to SBML code for the expression"""
    return "<cn>" + str(self.number) + "</cn>"
  def show_expr(self):
    """Display the underlying number of the numerical expression"""
    return str(self.number)

  def get_value(self, environment=None):
    """Returns the underlying value of this expression"""
    return self.number

  def reduce(self, environment=None):
    """Returns a reduced expression, but this is as far as we can reduce it"""
    return self

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

  def visit (self, visitor):
    """Implements the visit method allowing ExpressionVisitors to work"""
    visitor.visit_NameExpression(self)


  def show_expr(self):
    """Format as a string the name expression"""
    return self.name

  def get_value(self, environment=None):
    """Evalutes this expression based on the given variable_dictionary,
       If the environment is given and defines the name which is this
       expression then we return whatever the environment defines this name
       to be. Otherwise, as per 'get_value' we raise the expception KeyError.
    """
    # This is simple, we just do whatever looking up in the environment
    # does, which is the advertised behaviour of get_value.
    if environment != None:
      return environment[self.name]
    else:
      raise KeyError(self.name)
     
  def reduce(self, environment=None):
    """Attempts to reduce this expression, since this is a name expression
       we can reduce it to a number if the value is in the constant
       value mapping provided, otherwise we just return ourselves.
    """
    try:
      value = self.get_value(environment=environment)
      return NumExpression(value)
    except KeyError:
      return self
      

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

def list_product(factors):
  """Simple utility the same as 'sum' but for the product of the arguments.
     Note: returns 1 for the empty list, which seems reasonable, given that
     sum([]) = 0.
  """
  result = 1
  for factor in factors:
    result *= factor
  return result

class ApplyExpression(Expression):
  """A class to represent the AST of an apply expression, applying a
     named function to a list of argument expressions"""
  def __init__(self, name, args):
    Expression.__init__(self)
    self.name = name
    self.args = args

  def visit (self, visitor):
    """Implements the visit method allowing ExpressionVisitors to work"""
    visitor.visit_ApplyExpression(self)


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
        raise ValueError("Square root function takes exactly one argument")
      argument = new_args[0]
      half_expr = NumExpression(0.5)
      new_expr = ApplyExpression("power", [argument, half_expr])
      return new_expr
    else:
      new_expr = ApplyExpression(self.name, new_args)
      return new_expr

  def get_value(self, environment=None):
    """Return the value to which the expression evaluates. If any
       environment is given it should be used to resolve any names in
       the sub-expressions of this expression. If the expression is
       irreducible, generally because it contains a name which is not
       in the given environment (or none is given) then None is returned.
    """
    arg_values = [ arg.get_value(environment=environment)
                     for arg in self.args ]
    if self.name == "plus":
      return sum(arg_values)
    elif self.name == "times":
      answer = 1
      for arg in arg_values:
        answer *= arg
      return answer
    elif self.name == "minus":
      # What should we do if there is only one argument, I'm not sure if
      # <apply><minus/><cn>1.0</cn></apply> is allowed but if so I'm guessing
      # it should evaluate to -1.0
      answer = arg_values[0]
      for arg in arg_values[1:]:
        answer -= arg
    elif self.name == "divide":
      answer = arg_values[0]
      for arg in arg_values[1:]:
        answer /= arg
    elif self.name == "power":
      # power is interesting because it associates to the right
      exponent = 1
      # counts downwards from the last index to the 0.
      # As an example, consider power(3,2,3), the answer should be
      # 3 ** (2 ** 3) = 3 ** 8 = 6561, not (3 ** 2) ** 3 = 9 ** 3 = 81
      # going through our loop here we have
      # exp = 1
      # exp = 3 ** exp = 3
      # exp = 2 ** exp = 2 ** 3 = 8
      # exp = 3 ** exp = 3 ** 8 = 6561
      for i in range(len(arg_values) - 1, -1, -1):
        exponent = arg_values[i] ** exponent
      return exponent 
    else:
      raise ValueError("Unknown function name: " + self.name)  

  def reduce(self, environment=None):
    """Attempts to reduce this expression, note that there are three
       possibilities, but the first is that the entire expression cannot be
       reduced any further, in which case we can just return the current
       expression, but we can ignore this possibility. Another possibility is
       that some of the argument expressions can be reduced but not all, in
       which case we return a new apply expression with the reduced (as far
       as they can be) argument expressions. Finally the case in which all
       expressions can be reduced, in which case we return the value.
    """
    # We can easily check the final case by simply calling 'get_value' on
    # this expression, if it doesn't reduce we can assume that at least
    # one argument expression doesn't reduce.
    # Note that this means we in some sense do a little of the work twice
    # in the worst case we may have many arguments which reduce to a value
    # and only one which does not, in which case we could have reduced
    # all arg expressions and then applied the function if they all reduced
    # to a NumExpression, otherwise build up the NameExpression with the
    # reduced expression. The point is that here we assume you are doing this
    # reduction once at the start of say a simulation and hence you don't
    # require this to be extremely fast, and this is a very nice definition
    # which means we need not write code to evaluate plus/minus etc twice.
    # Alternatively we could write 'get_value' in terms of 'reduce' but that
    # would mean building up NumExpressions more than we needed to.
    try:
      value = self.get_value(environment=environment)
      return NumExpression(value)
    except KeyError:
      # We have a special case for the commutative expression plus and times,
      # here we try to pull out all of the argument expressions which reduce
      # to a value and sum or product them together. This is simply so that
      # we can turn the expression:
      # R * 0.2 * 10 where R is a dynamic variable into the expression
      # R * 2
      # which may save a bit of time during the simulation.
      if self.name == "plus" or self.name == "times":
        factors = []
        arg_expressions = []
        for arg in self.args:
          try:
            factors.append(arg.get_value(environment=environment))
          except KeyError:
            arg_expressions.append(arg.reduce(environment=environment))
        
        # Based on whether it's plus or times we must (possibly) add the
        # correct factor argument into the list of argument_expressions.
        # At the end of this conditional arg_expressions will be the
        # correct set of argument expressions.
        if self.name == "plus":
          factor = sum(factors)
          # So if factor is not zero then we must add it as an arg expr.
          if factor != 0:
            factor_exp = NumExpression(factor)
            arg_expressions = [factor_exp] + arg_expressions
        else: # assume it equals "times", we above check this.
          factor = list_product(factors)
          if factor != 1:
            factor_exp = NumExpression(factor)
            arg_expressions = [factor_exp] + arg_expressions
        # Now that we have the correct set of argument expressions
        # we may return the reduced apply expression, but we first just
        # check that we have not reduced it to a single expression, in which
        # case we can simply return that expression, eg we may have started
        # with R * 0.1 * 10, which would reduce to R * 1, but since then the
        # factor would be 1 we would not have added it as an arg_expression.
        if len(arg_expressions) == 1:
          return arg_expressions[0]
        else:
          return ApplyExpression(self.name, arg_expressions)
          
      else:
        # The easy case for non-commutative, we could go deeper and
        # try to partially evaluate some of these, for example
        # R - 3  - 1 could be made into R - 4. But for now the above will
        # do, since I believe that multiplications by more than one constant
        # are fairly common.
        arg_expressions = [ arg.reduce(environment=environment)
                              for arg in self.args ]
        return ApplyExpression(self.name, arg_expressions)
     

class ExpressionVisitor(object):
  """A parent class for classes which descend through the abstract syntax
     of expressions, generally storing a result along the way.
  """
  def __init__(self):
    self.result = None

  def generic_visit(self, expression):
    """The main entry for visiting generic expression whose type we do
       not yet know, this is the most klutchy part of this, but there is
       no way around this.
    """
    expression.visit(self)

  def generic_visit_get_results(self, expression):
    """Performs the visit and also returns the result, sort of useful
       for doing this within a list comprehension.
    """
    self.generic_visit(expression)
    return self.result

  ###################################
  # These are the unimplemented methods that you would be likely
  # to override for your expression visitor.
  # pylint: disable=C0103
  def visit_NumExpression(self, _expression):
    """Visit a NumExpression element"""
    # pylint: disable=R0201
    message = "visit_NumExpression element for expression visitor"
    raise NotImplementedError(message)

  def visit_NameExpression(self, _expression):
    """Visit a NameExpression"""
    # pylint: disable=R0201
    message = "visit_NameExpression element for expression visitor"
    raise NotImplementedError(message)

  def visit_ApplyExpression(self, _expression):
    """Visit an ApplyExpression element"""
    # pylint: disable=R0201
    message = "visit_ApplyExpression element for expression visitor"
    raise NotImplementedError(message)



