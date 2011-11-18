"""A module implementing a transformation from biopepa to xml"""

import sys
import argparse
import xml.dom.minidom as minidom
import parcon

from biopepa_parser import VariableDeclaration, RateDefinition, \
                           ComponentDefinition
import biopepa_parser
import utils
import outline_sbml

def create_compartment_elements(document, model_element):
  """This is likely to change signature when we allow for user defined
     compartments but for now we have only a single default compartment
     to add to the model"""
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
  for comp_def in component_defs:
    species_element = document.createElement("species")
    name = comp_def.get_name()
    species_element.setAttribute("id", name)
    species_element.setAttribute("name", name)
    species_element.setAttribute("compartment", comp_def.get_location())
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

class Reaction(outline_sbml.Reaction):
  """A class representing reactions within a Bio-PEPA model. These
     reactions have been derived from the rate laws and the component
     definitions (and their associated reaction behaviours)
  """
  def __init__(self, name):
    super(Reaction, self).__init__(name)

  @property
  def rate(self):
    """The Bio-PEPA code essentially called the kinetic_law part
       the rate part so we set up a property to mimic that.
    """
    return self.kinetic_law
  @rate.setter
  def rate(self, rate):
    """The setter for the rate property"""
    self.kinetic_law = rate

  def get_mass_action_participants(self):
    """Return the left hand side participants which contribute to the
       rate in a mass action rate method. For example A + B --> C 
       would return A and B. Essentially this is so that fMA(r) could
       be translated to r * A * B.
    """
    # I think we return all the modifiers but perhaps not the inhibitors?
    return self.reactants + self.modifiers


  def add_behaviour(self, comp_name, behaviour):
    """Add the given behaviour to the reaction, essentially this means
       we are adding a component as a participant (reactant, product
       or modifier) to the reaction.
    """
    stoichimetries = behaviour.stoichiometry
    operators = behaviour.operator
    if len(stoichimetries) != len(operators):
      print ("Must have the same number of operators as stoichimetries")
      sys.exit(1)

    # Not sure if this is really the best way to do this, this means
    # we will be adding multiple behaviours for the same component.
    # eg, we'll have A + B -> B + B, where as arguably we should retain
    # the knowledge, perhaps we should just have the operator <> or ><
    for stoich, operator in zip(stoichimetries, operators):
      participant = outline_sbml.ReactionParticipant(comp_name)
      participant.stoich = stoich
      if operator == ">>" :
        self.products.append(participant)
      elif operator == "<<": 
        self.reactants.append(participant)
      elif operator == "(+)":
        self.modifiers.append(participant)
      elif operator == "(-)":
        self.modifiers.append(participant)
      elif operator == "(.)":
        self.modifiers.append(participant)
      else:
        print ("Unrecognised behaviour operator: " + operator)
        sys.exit(1)

  def create_element(self, document):
    """Create an xml element representing this reaction"""
    reaction_element = document.createElement("reaction")
    reaction_element.setAttribute("id", self.name)
    reaction_element.setAttribute("reversible", "false")
    reaction_element.setAttribute("fast", "false")
    if self.location:
      reaction_element.setAttribute("compartment", self.location)
    if self.reactants:
      list_of_reactants = document.createElement("listOfReactants")
      reaction_element.appendChild(list_of_reactants)
      for reactant in self.reactants:
        spec_ref = document.createElement("speciesReference")
        # spec_ref.setAttribute("id", reactant.name)
        spec_ref.setAttribute("species", reactant.name)
        # for biopepa models the stoichiometry values cannot
        # change during the simulation so this 'constant' attribute is
        # always true.
        spec_ref.setAttribute("constant", "true")
        spec_ref.setAttribute("stoichiometry", str(reactant.stoich))
        list_of_reactants.appendChild(spec_ref)
    if self.products:
      list_of_products = document.createElement("listOfProducts")
      reaction_element.appendChild(list_of_products)
      for product in self.products:
        spec_ref = document.createElement("speciesReference")
        # spec_ref.setAttribute("id", product.name)
        spec_ref.setAttribute("species", product.name)
        # for biopepa models the stoichiometry values cannot
        # change during the simulation so this 'constant' attribute is
        # always true.
        spec_ref.setAttribute("constant", "true")
        spec_ref.setAttribute("stoichiometry", str(product.stoichiometry))
        list_of_products.appendChild(spec_ref)
    if self.modifiers:
      list_of_modifiers = document.createElement("listOfModifiers")
      reaction_element.appendChild(list_of_modifiers)
      for modifier in self.modifiers:
        spec_ref = document.createElement("modifierSpeciesReference")
        # spec_ref.setAttribute("id", modifier.name)
        # Note there is no stoichiometry attribute on
        # modifierSpeciesReference elements in SBML.
        spec_ref.setAttribute("species", modifier.name)
        list_of_modifiers.appendChild(spec_ref)
    if self.rate:
      kinetic_law = document.createElement("kineticLaw")
      reaction_element.appendChild(kinetic_law)
      math_element = document.createElement("math")
      kinetic_law.appendChild(math_element)
      mathxmlns = "http://www.w3.org/1998/Math/MathML"
      math_element.setAttribute("xmlns", mathxmlns)
      simplified_rate = self.rate.remove_rate_law_sugar(self)
      expr_element = simplified_rate.create_sbml_element(document)
      math_element.appendChild(expr_element)
    return reaction_element



def build_reaction_dictionary(components, rate_definitions):
  """From a list of components build a reaction dictionary which maps
     reaction names to reaction representations"""
  reaction_dictionary = dict()
  for rate_def in rate_definitions:
    name = rate_def.get_name()
    reaction = Reaction(name)
    reaction.rate = rate_def.get_rate()
    reaction_dictionary[name] = reaction

  def have_same_locations(reaction, reaction_behaviour):
    """returns true if the reaction has the same location as
       the reaction behaviour
    """
    return reaction.location == reaction_behaviour.location

  def new_reaction_for_behaviour(behaviour):
    """If there is no existing reaction to which to add a given
       behaviour then call this method to create a new such reaction
       and add it to the reaction dictionary
    """
    reaction = Reaction(b_name)
    reaction_dictionary[b_name] = reaction
    if behaviour.location:
      reaction.set_location(behaviour.location)
    return reaction

  for component_def in components:
    for behaviour in component_def.get_behaviours():
      b_name = behaviour.get_name() 
      # So basically we add the current behaviour to a reaction
      # if we have already created a reaction with the given name
      # at the given location (if any) of the behaviour. Otherwise
      # we create a new reaction and add it to the dictionary, but either
      # way we add the current behaviour to some reaction.

      # Ach this totally won't work if we have the same reaction name
      # with different locations because we will overwrite the mapping
      # with the second (and subsequent) ones we encounter.

      if b_name in reaction_dictionary:
        reaction = reaction_dictionary[b_name]
        if not have_same_locations(reaction, behaviour):
          # So there exists a reaction with the correct name, but not
          # at the correct location
          reaction = new_reaction_for_behaviour(behaviour)
      else:
        # There is no reaction with the correct name, regardless of
        # the location of such a reaction.
        reaction = new_reaction_for_behaviour(behaviour)
      # Whether we have created a new reaction or we are adding
      # behaviour to an existing one we must add the current
      # behaviour to the reaction.
      reaction.add_behaviour(component_def.get_name(), behaviour)
  return reaction_dictionary


def translate_biopepa_model(parse_result):
  """Translate the parsed Bio-PEPA model into an SBML xml document"""
  var_decs = [ x for x in parse_result.definitions
                     if isinstance(x, VariableDeclaration) ]
  rate_defs = [ x for x in parse_result.definitions
                      if isinstance(x, RateDefinition) ]
  components = [ x for x in parse_result.definitions
                       if isinstance(x, ComponentDefinition) ]

  reaction_dictionary = build_reaction_dictionary(components, rate_defs)

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

  create_species_elements(document, model_element,
                          parse_result.system_equation)

  convert_variable_declarations(document, model_element, var_decs)

  list_of_reactions = document.createElement("listOfReactions")
  model_element.appendChild(list_of_reactions)

  for reaction in reaction_dictionary.values():
    reaction_element = reaction.create_element(document)
    list_of_reactions.appendChild(reaction_element)

  return document


def process_file(filename, arguments):
  """Parse in a Bio-PEPA file, translate to SBML and create an
     SBML file which should be the translated Bio-PEPA model.
  """
  # I'm not exactly sure how to do this, there is no parseFile
  # in parcon, I think we just need to open the file and pass in
  # the file handle but I haven't tried that yet. 
  model_file = open(filename, "r")
  try:
    parse_result = biopepa_parser.parse_model_file(model_file)
  except parcon.ParseException as parse_except:
    print parse_except
    sys.exit(1)
  model_file.close()

  document = translate_biopepa_model(parse_result)

  if arguments.output_file:
    sbml_filename = arguments.output_file
  else:
    sbml_filename = utils.change_filename_ext(filename, ".sbml")
  if sbml_filename == "stdout":
    sbml_file = sys.stdout
  else:
    sbml_file = open(sbml_filename, "w")

  document.writexml(sbml_file, indent="", addindent="  ", newl="\n")
  if sbml_filename != "stdout":
    sbml_file.close()

def main():
  """A simple main function to parse in the arguments as
     Bio-PEPA model files"""
  description = "Parse a Bio-PEPA model file(s)"
  parser = argparse.ArgumentParser(add_help=True,
                                   description=description)
   # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="A Bio-PEPA model file to parse")
  parser.add_argument('--output-file', action='store',
                      help="The file to output to, stdout to print")
  arguments = parser.parse_args()
  for filename in arguments.filenames:
    process_file(filename, arguments)


if __name__ == '__main__':
  main()


