""" A simple script to replace the values of parameters within a
    biopepa file. Currently only works for parameters defined on
    a single line.
"""
     
import sys
import argparse
import logging

import utils
import parameters

def parameterise_model_file(dictionary, biopepa_filename, output_filename):
  """Parameterise a model file with the given dictionary"""
  biopepa_file = open(biopepa_filename, "r")
  if output_filename == "stdout":
    output_file = sys.stdout
  else:
    output_file = open(output_filename, "w")
  parameterised_names = []

  for line in biopepa_file:
    if '=' in line:
      name = line.partition('=')[0]
      name = name.lstrip().rstrip()
      if dictionary.has_key (name):
        new_value = dictionary[name]
        output_file.write(name + " = " + str(new_value) + " ;\n")
        parameterised_names.append(name)
      else:
        output_file.write(line)
    else:
      output_file.write(line)

  unparameterised = [ x for x in dictionary.keys()
                          if x not in parameterised_names ]
  if unparameterised:
    logging.error ("Failed to fully parameterise Bio-PEPA file")
    sys.exit(1)

  biopepa_file.close()
  output_file.close()


def run():
  """Perform the banalities of command-line argument processing and
     and then get under way in parameterising the model"""
  description = "Parameterise an SBML model based on a given param file"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="Bio-PEPA and parameter files")

  arguments = parser.parse_args()

  biopepa_extentions = [ ".biopepa" ]
  param_files = [ x for x in arguments.filenames
                        if not utils.has_extension(x, biopepa_extentions) ]

  biopepa_files = [ x for x in arguments.filenames
                      if utils.has_extension(x, biopepa_extentions) ]

  dictionary = dict()
  for param_file in param_files:
    parameters.parse_param_file(param_file, dictionary)
  for biopepa_file in biopepa_files:
    parameterise_model_file(dictionary, biopepa_file, "stdout")

if __name__ == "__main__":
  run()
