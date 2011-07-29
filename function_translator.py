"""A simple utility which uses the biopepa parser to parse in function
definitions and output them as sbml function definitions.
"""
import argparse
import pyparsing
from pyparsing import ( Literal, Word, OneOrMore, Optional, 
                        Combine, nums, Or
                      )
import biopepa.biopepa_parser as bparser


class FunctionDefinition:
  def __init__(self, name, parameters, body):
    self.name = name
    self.parameters = parameters
    self.body = body

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

parameter_list_parser = pyparsing.ZeroOrMore(bparser.variable_name)

function_def_parser = ( bparser.variable_name + 
                        parameter_list_parser +
                        "=" + 
                        bparser.expression_parser +
                        ";"
                      )
def make_function_definition(tokens):
  name = tokens[0]
  # The last token is the semi-colon, the second last is the body
  # expression, the third last then is the equals sign
  equals_index = len(tokens) - 3
  parameters = tokens[1:equals_index]
  body = tokens[equals_index + 1]
  return FunctionDefinition(name, parameters, body)
function_def_parser.setParseAction(make_function_definition)


function_def_list_parser = (OneOrMore(function_def_parser) + 
                            pyparsing.StringEnd())


def parse_function_list(source):
  """A significant entry point into this model, parse all of the
     string as though it were the contents of a file containing a list
     of function definitions
  """
  return function_def_list_parser.parseString(source)

def process_file(filename):
  """A simple method to process function definition file and
     output the corresponding SBML
  """
  parse_result = function_def_list_parser.parseFile(filename)

  # print (parse_result.asXML())
  for fun_def in parse_result:
    print (fun_def.show_sbml())
  
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
  for filename in arguments.filenames:
    process_file(filename)

if __name__ == '__main__':
  main()
