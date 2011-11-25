"""A module implementing a transformation from facile to sbml"""

import sys
import argparse
import xml.dom.minidom as minidom
import parcon

import facile_parser
import utils
import outline_sbml
import create_sbml

def translate_facile_model(parse_result):
  """Translate the parsed facile model into an SBML xml document"""
  sbml_model = create_sbml.SBML_Model()

  sbml_model.reactions = parse_result # [0]
  document = sbml_model.create_sbml_document()
  return document

def translate_file(filename, arguments):
  """Translate a facile file into sbml"""
  model_file = open(filename, "r")
  try:
    parse_result = facile_parser.parse_model_file(model_file)
  except parcon.ParseException as parse_except:
    print parse_except
    sys.exit(1)
  model_file.close()

  document = translate_facile_model(parse_result)
  create_sbml.output_to_sbml_file(filename, arguments, document)


def main():
  """A simple main function to parse in the arguments
     as facile model file
  """
  description = "Translate facile model file(s) into sbml"
  parser = argparse.ArgumentParser(add_help=True,
                                   description=description)
   # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="A facile model file to translate")
  parser.add_argument('--output-file', action='store',
    help="""Specify an output filename; "stdout" to print to the console""")
  arguments = parser.parse_args()

  for filename in arguments.filenames:
    translate_file(filename, arguments)


if __name__ == '__main__':
  main()


