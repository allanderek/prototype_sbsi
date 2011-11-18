""" A simple script to replace the values of parameters within a
    biopepa file. Currently only works for parameters defined on
    a single line.
"""
     
import sys
import os
import argparse
import logging

def parameterise_model_file(filename, dictionary):
  """Parameterise a model file with the given dictionary"""
  biopepa_file = open(filename, "r")
  output_file = sys.stdout
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


def parse_param_file(param_filename, dictionary):
  """Parse a parameter file into the given dictionary"""
  param_file = open(param_filename, "r")

  for line in param_file:
    separator = "\t"
    if not separator in line and ":" in line:
      separator = ":"
    (name, value_string) = line.split(separator, 1)
    value = float(value_string.lstrip().rstrip())
    # We should also check if the name already exists in the dictionary
    dictionary[name] = value

  param_file.close()

def has_extension(filename, extensions):
  """Returns true if the given filename has one of the given
     list of extensions"""
  file_extension = os.path.splitext(filename)[1]
  return file_extension in extensions

def run():
  """Perform the banalities of command-line argument processing and
     and then get under way in parameterising the model"""
  description = "Parameterise an SBML model based on a given param file"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="input files, parameters and sbml model files")
  # parser.add_argument('--pretty', action='store_true',
  #                     help="Pretty print the xml")

  arguments = parser.parse_args()

  biopepa_extentions = [ ".biopepa" ]
  param_files = [ x for x in arguments.filenames
                        if not has_extension(x, biopepa_extentions) ]

  biopepa_files = [ x for x in arguments.filenames
                      if has_extension(x, biopepa_extentions) ]

  dictionary = dict()
  for param_file in param_files:
    parse_param_file(param_file, dictionary)
  for biopepa_file in biopepa_files:
    parameterise_model_file(biopepa_file, dictionary)

if __name__ == "__main__":
  run()
