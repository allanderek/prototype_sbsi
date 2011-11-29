"""A script to create a C model suitable for use with SBSI numerics"""
import xml.dom.minidom
import argparse

import outline_sbml


def convert_file(filename):
  """Given an sbml file convert it into the corresponding UserModel.{ch}
     files for use with sbsi numerics
  """
  dom = xml.dom.minidom.parse(filename)
  model = dom.getElementsByTagName("model")[0]
  species = outline_sbml.get_list_of_species(model)

  model_name = model.getAttribute("name")
  if not model_name:
    model_name = model.getAttribute("id")
  if not model_name:
    model_name = "nameless_model"

  c_file = open ("UserModel.C", "w")
  c_file.write(
  """#include <cmath>
using namespace std;

#include <cModel.h>
#include "UserModel.h"
#define pi M_PI 
  """)

  c_file.write("/**** The model has " + str(len(species)) + 
               " species ****/\n")
  c_file.write("void " + model_name + "::inputModel(){\n")

  # For each species 
  for component in species:
    if component.initial_amount:
      c_file.write ("  setStateByName(\"" + component.ident +
                    "\", " + str(component.initial_amount) + ");\n")
    else:
      c_file.write ("  // Setting the state to zero as it should be\n")
      c_file.write ("  // updated by an initial assignment\n")
      c_file.write ("  setStateByName(\"" + component.ident +
                    "\", 0);\n")

  # We also need to do the same for variables and parameters
  parameters = outline_sbml.get_list_of_parameters(model)
  for param in parameters:
    if param.value:
      c_file.write("  setParamByName(\"" + param.ident +
                   "\", " + param.value + ");\n")
    else:
      c_file.write("  // Parameter's value not set in model file\n")
      c_file.write("  // Presumably set elsewhere, eg initialAssignment\n")
      c_file.write("  setParamByName(\"" + param.ident +
                   "\", 0);\n")
      
   

  c_file.write(
"""
	/*  0th Species with Id A  in main */   
	// Setting the state to zero as it should be 
	// updated by an initial assinment
	setStateByName("A",0);

	/*  1th Species with Id B  in main */   
	// Setting the state to zero as it should be 
	// updated by an initial assinment
	setStateByName("B",0);

	/*  2th Species with Id AB  in main */   
	// Setting the state to zero as it should be 
	// updated by an initial assinment
	setStateByName("AB",0);


	setVarByName("AB",0);

	setParamByName("k1",1);
	setParamByName("k2",1);
	setParamByName("k3",1);
	setParamByName("main",1);
	};
""")


  c_file.write("\n")
  c_file.close()
  
 
def run():
  """Perform the banalities of command-line argument processing
     and then do the actual work on each argument file
  """
  description = "Analyse SBML files for invariants"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="an sbml file to check invariants for")

  arguments = parser.parse_args()

  for filename in arguments.filenames:
    # Since it always converts to UserModel/UserModel.{ch}
    # more than one file will just overwrite the others, so
    # perhaps we need to think about this somewhat more.
    convert_file(filename)

if __name__ == "__main__":
  run()
