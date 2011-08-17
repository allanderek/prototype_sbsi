"""A simple utility which uses the biopepa parser to parse in function
definitions and output them as sbml function definitions.
"""
import os.path
import argparse
import xml.dom.minidom
from parcon import Forward, InfixExpr, Translate, Optional, ZeroOrMore
import parcon
import biopepa.biopepa_parser as bparser

class FunctionDefinition:
  def __init__(self, name, parameters, body):
    self.name = name
    self.parameters = parameters
    self.body = body
    self.math_ns = "http://www.w3.org/1998/Math/MathML"

  def show_sbml(self):
    result = "<functionDefinition id=\"" + self.name + "\" >\n"
    result += "<math xmlns=\"http://www.w3.org/1998/Math/MathML\"\n"
    result += "      xmlns:sbml=\"http://www.sbml.org/sbml/level3/"
    result += "version1/core\">"
    result += "<lambda>\n"

    for parameter in self.parameters:
      result += "<bvar><ci>" + parameter + "</ci></bvar>"

    result += self.body.convert_to_sbml()

    result += "</lambda>\n"
    result += "</math>\n"
    result += "</functionDefinition>"
    return result

  def create_element(self, document):
    fun_def = document.createElement("functionDefinition")
    fun_def.setAttribute("id", self.name)

    e_math = document.createElementNS(self.math_ns, "math")
    e_math.setAttribute("xmlns", self.math_ns)
    fun_def.appendChild(e_math)

    lambda_el = document.createElement("lambda")
    e_math.appendChild(lambda_el)

    for parameter in self.parameters:
      bvar_el = document.createElement("bvar")
      lambda_el.appendChild(bvar_el)
      ci_el = document.createElement("ci")
      bvar_el.appendChild(ci_el)
      ci_text = document.createTextNode(parameter)
      ci_el.appendChild(ci_text)

    body_el = self.body.create_sbml_element(document)
    lambda_el.appendChild(body_el)
   
    return fun_def
 
  def show_fun_def (self):
    result = self.name + " "
    for param in self.params:
      result += param + " "
    result += " = "
    result += self.body.show_expr()
    return result
   
def make_fun_def(parse_result):
  return FunctionDefinition(parse_result[0],
                            parse_result[1],
                            parse_result[2])

parameter_list = ZeroOrMore(parcon.alpha_word)
fun_def = parcon.alpha_word + parameter_list + "=" + bparser.expr + ";"
fun_def_parser = Translate (fun_def, make_fun_def)

function_def_list_parser = ( parcon.OneOrMore(fun_def_parser) +
                             parcon.End()
                           )


def parse_function_list(source):
  """A significant entry point into this model, parse all of the
     string as though it were the contents of a file containing a list
     of function definitions
  """
  return function_def_list_parser.parse_string(source)

def add_fun_defs_to_sbml(fun_defs, sbml_file):
  dom = xml.dom.minidom.parse(sbml_file)
  model = dom.getElementsByTagName("model")[0]
  
  # Where loe = listOfEvents
  loe_elements = model.getElementsByTagName("listOfFunctionDefinitions")
  if not loe_elements:
    loe_element = dom.createElement("listOfFunctionDefinitions")
    # In general when placing nodes with the model we must take care
    # of where in the list we place them because in sbml the child elements
    # of the model element are in a specified order. So in general one
    # should get the list of child elements and then progress through
    # it until we recognise we are in the correct place and then insert
    # there. In this instance because 'listOfFunctionDefinitions' is the
    # first element we can just insert it at the front, which we can do
    # by obtaining whatever is the first element.
    first_child = model.firstChild
    model.insertBefore(loe_element, first_child)
  else:
    loe_element = loe_elements[0]

  for fun_def in fun_defs:
    loe_element.appendChild(fun_def.create_element(dom))

  # if arguments.pretty:
  document = dom.toprettyxml(indent="  ", encoding="UTF-8")
  # else:
  #  document = dom.toxml("UTF-8")
  print(document)



def process_file(filename, sbml_files):
  """A simple method to process function definition file and
     output the corresponding SBML
  """

  fun_file = open(filename, "r")
  contents = fun_file.read()
  parse_result = function_def_list_parser.parse_string(contents)
  fun_file.close()

  if not sbml_files:
    # print (parse_result.asXML())
    for fun_def in parse_result:
      print (fun_def.show_sbml())
  else:
    for sbml_file in sbml_files:
      add_fun_defs_to_sbml(parse_result, sbml_file)
  
def main():
  """A simple main function to parse in the arguments as files containing
     human readable function definitions and print those function
     definitions out as sbml function definitions.
  """
  description = "Convert human readable functions to sbml functions"
  parser = argparse.ArgumentParser(add_help=True,
                                   description=description)
   # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="A function definition file to translate")
  arguments = parser.parse_args()

  sbml_files = [ sbml for sbml in arguments.filenames 
                        if os.path.splitext(sbml)[1] in [".xml", ".sbml" ]
               ]

  # Doesn't work!
  # bparser.expr.graph().draw("syntax-expr.png")

  for filename in arguments.filenames:
    if filename not in sbml_files :
      process_file(filename, sbml_files)

if __name__ == '__main__':
  main()
