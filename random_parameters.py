"""A very simple module to generate some random parameters for an
   sbml file. Generally this is useful for setting up an
   optimisation challenge which you might test your model on to see
   how identifiable it is
"""
import argparse
import sys

import parameters
import guess_initparams
import utils

def run():
  """Do the simple task of generating a bunch of random parameters
     for an sbml model. The input can either be the sbml model or
     an initparams file
  """
  description = "Create random parameters from an init params file"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="An initial parameters file")
  arguments = parser.parse_args()

  if len(arguments.filenames) < 1:
    print ("Must provide an initial parameters file")
    sys.exit(1)

  filename = arguments.filenames[0]
  if utils.has_xml_or_sbml_ext(filename):
    params = guess_initparams.init_params_model_file(filename)
  else:
    params = parameters.get_init_param_parameters(filename)

  for param in params:
    print(param.name + "\t" + str(param.get_random_value_full_range())) 

if __name__ == "__main__":
  run()
