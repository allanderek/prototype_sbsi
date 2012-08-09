"""A module implementing a formatting translation from Bio-PEPA
   to the LaTeX document processing system
"""

import sys

import biopepa.biopepa_parser as biopepa_parser
import utils
import expressions

class ExpressionLatexifer(expressions.ExpressionVisitor):
  """A parent class for classes which descend through the abstract syntax
     of expressions, generally storing a result along the way.
  """
  def __init__(self):
    super(ExpressionLatexifer, self).__init__()

  ###################################
  # These are the unimplemented methods that you would be likely
  # to override for your expression visitor.
  # pylint: disable=C0103
  def visit_NumExpression(self, expression):
    """Visit a NumExpression element"""
    self.result = str(expression.number)

  def visit_NameExpression(self, expression):
    """Visit a NameExpression"""
    self.result = expression.name

  def visit_ApplyExpression(self, expression):
    """Visit an ApplyExpression element"""
    result = expression.name + " ("
    arg_strings = [ self.generic_visit_get_results(arg)
                      for arg in expression.args ]
    result += ", ".join(arg_strings)
    result += ")"
    self.result = result

def latexify_expression(expression):
  """Produce a latex string representing the given expression"""
  latexifier = ExpressionLatexifer()
  latexifier.generic_visit(expression)
  return latexifier.result


def latexify_operator(operator):
  """Convert the ascii operator into a latex equivalent"""
  if operator == "<<":
    return "\\reactant"
  elif operator == ">>":
    return "\\product"
  elif operator == "(+)":
    return "\\activator"
  elif operator == "(-)":
    return "\\inhibitor"
  elif operator == "(.)":
    return "\\modifier"
  # And just to allow unknown operators, remember this is Bio-PEPA to
  # LaTeX so there is no guarantee we particularly wish to only allow
  # valid Bio-PEPA.
  else:
    return operator


def latexify_behaviour(behaviour, component_name):
  """Print out the behaviour as a string in Bio-PEPA format"""
  if behaviour.stoichiometry == ['1']:
    result = (behaviour.reaction_name + " " +
              latexify_operator(behaviour.operator[0]))
  else:
    stoichs = ", ".join([str(s) for s in behaviour.stoichiometry])
    opers = " ".join([ latexify_operator(op)
                         for op in behaviour.operator])
    result =  ("(" + behaviour.reaction_name + ", " + 
           stoichs + ") " + opers)

  if behaviour.location != None:
    result += " " + component_name + "@" + behaviour.location

  return result


def translate_biopepa_model(model, out_file):
  """Translate the parsed Bio-PEPA model into a LaTeX document"""

  # The variable declarations
  if model.var_decs != None and model.var_decs:
    out_file.write ("\\section{Variable Declarations}\n\n")
    out_file.write ("\\begin{gather*}\n")
    for var_dec in model.var_decs:
      out_file.write(var_dec.variable)
      out_file.write(" & = & ")
      out_file.write(latexify_expression(var_dec.expression))
      out_file.write("\\\\\n")
    out_file.write("\\end{gather*}\n\n\n")

  # The rate definitions
  if model.rate_defs != None and model.rate_defs:
    out_file.write("\\section{Rate Definitions}\n\n")
    out_file.write("\\begin{gather*}\n")
    for rate_def in model.rate_defs:
      out_file.write(rate_def.name)
      out_file.write (" & = & [ ")
      out_file.write(latexify_expression(rate_def.rate))
      out_file.write (" ] \\\\\n")
      
    out_file.write("\\end{gather*}\n\n\n")

  # Now on to component definitions 
 
  
  if model.component_defs != None and model.component_defs:
    out_file.write("\\section{Component Definitions}\n\n")
    out_file.write("\\begin{gather*}\n")
    for comp_def in model.component_defs:
      out_file.write(comp_def.name)
      out_file.write (" & = &  ")
      behaviour_strings = [ latexify_behaviour(b, comp_def.name)
                              for b in comp_def.behaviours ]
      out_file.write(" + ".join(behaviour_strings))
      out_file.write (" \\\\\n")
      
    out_file.write("\\end{gather*}\n\n\n")

  if model.system_equation != None and model.system_equation:
    comp_pops = [ comp_pop.format() for comp_pop in model.system_equation ]
    system_equation = " <*>\n".join(comp_pops)

    out_file.write ("\\section{The system equation}\n\n")
    out_file.write (system_equation)
    out_file.write ("\n")


def convert_source(source):
  """Converts the source of a Bio-PEPA model, represented as a string,
     into a string representing the LaTeX source
  """
  model = biopepa_parser.parse_model(source)
  outfile_imposter = utils.StringFile()
  translate_biopepa_model(model, outfile_imposter)
  return outfile_imposter.get_results()
  
def latexify_file(filename, arguments):
  """ Converts a file containing Bio-PEPA source into a file containing
      LaTeX source to format the same model.
  """
  model_file = open(filename, "r")
  parse_result = biopepa_parser.parse_model_file_exit_on_error(model_file)
  model_file.close()

  biopepa_model = parse_result

  out_filename = utils.get_output_filename(filename, arguments, ".tex")
  if out_filename == "stdout":
    translate_biopepa_model(biopepa_model, sys.stdout)
  else:
    out_file = open(out_filename, "w")
    translate_biopepa_model(biopepa_model, out_file)
    out_file.close()

def main():
  """The main work of processing the command-line arguments and then
     converting the Bio-PEPA file into LaTeX
  """
  description = "Translates Bio-PEPA model files to LaTeX"
  file_help = "A model file to translate to LaTeX"
  parser = utils.argparser_withfiles(description, file_help, True)
  utils.add_output_file_arg(parser)
  arguments = parser.parse_args()

  for biopepa_file in arguments.filenames:
    latexify_file(biopepa_file, arguments)


if __name__ == '__main__':
  main()


