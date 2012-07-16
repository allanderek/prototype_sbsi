"""A script intended to be able to check the rates for all reactions
   in an SBML model. We are looking for surprising references to
   non-reactant species, or surprisiring ommissions of references to
   reactant species
"""
import xml.dom.minidom
import argparse

import rate_checker_sbml

def analyse_model(model, arguments):
  """Perform static analysis over an SBML model"""
  num_warnings = rate_checker_sbml.check_rates_sbml_model(model, arguments)

  # num_warnings += other analyses
  return num_warnings

def analyse_document(dom, arguments):
  """Perform static analysis over an SBML document"""
  model = dom.getElementsByTagName("model")[0]
  return analyse_model(model, arguments)

def analyse_file(filename, arguments):
  """Perform analysis over the rate definitions of an SBML file"""
  dom = xml.dom.minidom.parse(filename)

  return analyse_document(dom, arguments)


def create_arguments_parser():
  """Create the argument parser, it's best to have this as a separate
     method such that other argument parsers can have this as a
     parent argument parser
  """
  description = "Statically analyse SBML files for modelling errors"
  parent_arg_parser = rate_checker_sbml.create_arguments_parser()
  parser = argparse.ArgumentParser(description=description,
                                   parents=[parent_arg_parser])
  return parser

def run():
  """perform the banalities of command-line argument processing and
     then go ahead and calculate the outline for each model file"""
  parser = create_arguments_parser()
  arguments = parser.parse_args()

  num_warnings = 0
  for filename in arguments.filenames:
    num_warnings += analyse_file(filename, arguments)
 
  if num_warnings == 0:
    print ("No warnings")
  elif num_warnings == 1:
    print ("There was a single warning")
  else:
    print ("There were " + str(num_warnings))

if __name__ == "__main__":
  run()
