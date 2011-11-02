"""A script intended to be able to check the rates for all reactions
   in an SBML model. We are looking for surprising references to
   non-reactant species, or surprisiring ommissions of references to
   reactant species
"""
import xml.dom.minidom
import argparse

import outline_sbml


def get_rate_affectors(rate_element):
  """Returns the identifiers used the rate name"""
  ci_elements = rate_element.getElementsByTagName("ci")
  identifiers = [] 
  for ident_element in ci_elements:
    name = ident_element.firstChild.data
    identifiers.append(name.lstrip().rstrip())

  return identifiers


def check_reaction(reaction, rate_analyser):
  """Check the kinetic law of a reaction"""
  num_warnings = 0
  kin_law = reaction.get_kinetic_law()
  # All the identifiers which affect the rate of the reaction.
  # Of course some of these may just be rate constants, some may not
  # be species. Additionally some may be dynamic variables which are
  # the sum of some species. So 'get_rate_affectors' should be a 'deep'
  # method which chases variable definitions. And here we should filter
  # out any that are not species names.
  species_names = [ unicode(s.get_name) for s in rate_analyser.species ]
  if not kin_law:
    return 0
  identifiers = rate_analyser.get_rate_affectors(kin_law)
  referenced_species = [ s for s in identifiers if s in species_names ]

  reactant_names = [ unicode(r.get_name()) 
                     for r in reaction.get_reactants() ]
  modifier_names = [ unicode(m.get_name())
                     for m in reaction.get_modifiers() ]
  # by this I mean all names on the left-hand side of the reaction
  all_lhs_names = reactant_names + modifier_names

  for reactant_name in all_lhs_names:
    if reactant_name not in identifiers: # referenced_species?
      num_warnings += 1
      print (reactant_name + " is a reactant of reaction " +
             reaction.get_name() + " but does not affect the rate")

  for identifier in referenced_species:
    if identifier not in all_lhs_names:
      num_warnings += 1
      print (identifier + " modifies the rate of reaction " +
             reaction.get_name() + " but is not a reactant or modifier")
    

  return num_warnings


class RateAnalyser:
  """A class implementing rate analysis for an sbml model"""
  def __init__(self, model):
    self.model = model
    self.species = outline_sbml.get_list_of_species(model)
    self.assignment_rules = outline_sbml.get_list_of_assignment_rules(model)
    self.assignments = dict()
    for assignment_rule in self.assignment_rules:
      name = assignment_rule.get_variable_name()
      self.assignments[name] = assignment_rule.get_assigned_expr()

  def get_rate_affectors(self, rate_element):
    """Returns the identifiers used within a rate (or math expression)
       element. This is deep in that it checks the names used for
       assignments which may themselve use other names.
    """
    ci_elements = rate_element.getElementsByTagName("ci")
    identifiers = []
    for ident_element in ci_elements:
      name = ident_element.firstChild.data
      identifiers.append(name.lstrip().rstrip())
  
    visit_stack = [ i for i in identifiers ]
    while visit_stack:
      name = visit_stack.pop()
      if name in self.assignments:
        expr = self.assignments[name]
        ci_elements = expr.getElementsByTagName("ci")
        for ident_element in ci_elements:
          new_name = ident_element.firstChild.data.lstrip().rstrip()
          if new_name not in identifiers:
            identifiers.append(new_name)
            visit_stack.append(new_name)

    return identifiers

   

def check_rates_sbml_file(filename):
  """Perform analysis over the rate definitions of an SBML model"""
  dom = xml.dom.minidom.parse(filename)
  model = dom.getElementsByTagName("model")[0]
  reactions = outline_sbml.get_list_of_reactions(model)
  rate_analyser = RateAnalyser(model)

  num_warnings = 0
 
  for reaction in reactions:
    num_warnings += check_reaction(reaction, rate_analyser)

  return num_warnings

 
  
def run():
  """perform the banalities of command-line argument processing and
     then go ahead and calculate the outline for each model file"""
  description = "Print out an outline of an SBML file"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  # This should instead just be a child of the outline parser.
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="an sbml file to outline")
  parser.add_argument("--ignore-sources",
                      action="store_true", default=False,
    help="Ingore source reactions")
  parser.add_argument("--ignore-sinks",
                      action="store_true", default=False,
    help="Ignore sink reactions")
  arguments = parser.parse_args()

  for filename in arguments.filenames:
    check_rates_sbml_file(filename)

if __name__ == "__main__":
  run()