"""
A module to ease creation of sbml documents
"""

import sys
import xml.dom.minidom as minidom

import utils

def output_to_sbml_file(filename, arguments, document):
  """A utility function for consumers of this module to calculate
     the output file for given a set of parsed arguments which includes
     the possible argument "--output-file". We then output the given
     document to that file. This is basically a suitable place to put
     this functionality that is common to nearly all x_to_sbml translators
     that we have.
     Note that the first argument is the filename of the file from which
     we have translated into sbml and not the sbml filename itself.
  """

  if arguments.output_file:
    sbml_filename = arguments.output_file
  else:
    sbml_filename = utils.change_filename_ext(filename, ".sbml")
  if sbml_filename == "stdout":
    sbml_file = sys.stdout
  else:
    sbml_file =  open(sbml_filename, "w")

  document.writexml(sbml_file, indent="", addindent="  ", newl="\n")
  if sbml_filename != "stdout":
    sbml_file.close()

 
def create_compartment_elements(document, model_element):
  """This is likely to change signature when we allow for user defined
     compartments but for now we have only a single default compartment
     to add to the model
  """
  list_of_compartments = document.createElement("listOfCompartments")
  model_element.appendChild(list_of_compartments)

  default_compartment = document.createElement("compartment")
  default_compartment.setAttribute("id", "default_compartment")
  default_compartment.setAttribute("constant", "true")
  list_of_compartments.appendChild(default_compartment)

def create_species_elements(document, model_element, component_defs):
  """Creates the listOfSpecies element and all of the associated
     species elements under it for a set of component definitions"""
  list_of_species = document.createElement("listOfSpecies")
  model_element.appendChild(list_of_species)
  if component_defs:
    for comp_def in component_defs:
      species_element = document.createElement("species")
      name = comp_def.name
      species_element.setAttribute("id", name)
      species_element.setAttribute("name", name)
      species_element.setAttribute("compartment", comp_def.location)
      species_element.setAttribute("hasOnlySubstanceUnits", "true")
      species_element.setAttribute("constant", "false")
      species_element.setAttribute("boundaryCondition", "false")
      list_of_species.appendChild(species_element) 

def convert_variable_declarations(document, top_element, var_decs):
  """Convert a set of variable declarations into SBML elements and
     add them to an sbml document
  """
  if var_decs:
    list_of_params = document.createElement("listOfParameters")
    top_element.appendChild(list_of_params)
 
    init_assigns = document.createElement("listOfInitialAssignments")
    top_element.appendChild(init_assigns)
    for var_dec in var_decs:
      param_element = var_dec.create_parameter_element(document)
      list_of_params.appendChild(param_element)

      init_assign = var_dec.create_initial_assignment(document)
      init_assigns.appendChild(init_assign)


class SBML_Model(object):
  """A class representing an SBML model which could be written
     out as an sbml document
  """
  def __init__(self):
    self.reactions = None
    self.component_defs = None
    self.var_decs = None
    

  def create_sbml_document(self):
    """Create an sbml model document for this SBML model"""
                          
    xml_implementation = minidom.getDOMImplementation()
    document = xml_implementation.createDocument(None, "sbml", None)
    top_element = document.documentElement 
    xmlns = "http://www.sbml.org/sbml/level3/version1/core" 
    top_element.setAttribute("xmlns", xmlns)
    top_element.setAttribute("level", "3")
    top_element.setAttribute("version", "1")

    model_element = document.createElement("model")
    top_element.appendChild(model_element)

    create_compartment_elements(document, model_element)
    create_species_elements(document, model_element, self.component_defs)
    convert_variable_declarations(document, model_element, self.var_decs)

    list_of_reactions = document.createElement("listOfReactions")
    model_element.appendChild(list_of_reactions)

    if self.reactions:
      for reaction in self.reactions:
        reaction_element = reaction.create_element(document)
        list_of_reactions.appendChild(reaction_element)

    return document

