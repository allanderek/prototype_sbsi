"""A module implementing a transformation from biopepa to xml"""

import sys
import argparse
import parcon

import utils
import biopepa.biopepa_parser as biopepa_parser
import sbml_ast

def flatten_model(loc_model):
  """Removes the locations from the given model"""
  def_map = dict()


  def munge_name(name):
    return name.replace("@", "__")

  # Rate definitions
  for rate_def in loc_model.rate_defs:
    rate_def.rate.munge_names(munge_name)
     

  # Component definitions
  for component_def in loc_model.component_defs:
    def add_behaviour(name, behaviour):
      if name in def_map:
        comp_def = def_map[name]
      else:
        comp_def = biopepa_parser.ComponentDefinition(name, [])
        def_map[name] = comp_def
    
      comp_def.behaviours.append(behaviour) 

    for behaviour in component_def.behaviours:
      participant_name = component_def.name 
      if behaviour.location != None:
        participant_name += "__" + behaviour.location
        behaviour.location = None
      add_behaviour(participant_name, behaviour)

  loc_model.component_defs = def_map.values()

  # Component populations are easy we just munge the name
  # in the same way as above.
  # TODO: population expressions should be name-munged as well
  for comp_pop in loc_model.system_equation:
    # This is horribly bad helicoptor code.
    # The default should None and then we should detect that.
    location = comp_pop.location
    if location != None and location != "default location":
      comp_pop.location = "default location"
      comp_pop.name = comp_pop.name + "__" + location


def process_file(filename, arguments):
  """Parse in a Bio-PEPA file, translate to SBML and create an
     SBML file which should be the translated Bio-PEPA model.
  """
  model_file = open(filename, "r")
  try:
    parse_result = biopepa_parser.parse_model_file(model_file)
  except parcon.ParseException as parse_except:
    print parse_except
    sys.exit(1)
  model_file.close()

  output_filename = utils.get_output_filename(filename, arguments,
                                              "_flattened.biopepa") 
  if output_filename != "stdout":
    output_file = open(output_filename, "w")
  else:
    output_file = sys.stdout


  flatten_model(parse_result)
  parse_result.output_model(output_file)

  if output_filename != "stdout":
    output_file.close()

  


def main():
  """A simple main function to parse in the arguments as
     Bio-PEPA model files"""
  description = "Flattens a Bio-PEPA model to remove compartments"
  parser = argparse.ArgumentParser(add_help=True,
                                   description=description)
   # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="A Bio-PEPA model file to flatten")
  utils.add_output_file_arg(parser)
  arguments = parser.parse_args()
  for filename in arguments.filenames:
    process_file(filename, arguments)


if __name__ == '__main__':
  main()


