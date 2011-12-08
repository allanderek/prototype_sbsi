"""A script to calculate the invariants of an SBML model"""
import xml.dom.minidom
import argparse

import outline_sbml
import timeseries
import utils

class KinecticIndependenceGraph:
  """A class representing the kinectic independence graph allowing
     the addition of reactions"""
  def __init__(self, species, reactions):
    self.reactions = reactions
    self.rows = dict()
    num_reactions = len(reactions)
    for spec in species:
      self.rows[spec] = [0] * num_reactions

  def add_reaction_info(self, reaction):
    """Add a reaction to the kinectic independence graph"""
    reactants = reaction.reactants
    column_index = self.reactions.index(reaction.name)
    for reactant in reactants:
      value = - reactant.get_stoichiometry()
      self.rows[reactant.get_name()][column_index] = value
    for product in reaction.products:
      value = product.get_stoichiometry()
      self.rows[product.get_name()][column_index] = value

  def get_rows_dictionary(self):
    """Return the rows of the kig as a dictionary mapping reaction
       names to a list of values, which correspond to how each
       species is affected by the given reaction"""
    return self.rows

  def get_species(self):
    """Return the set of species within the kig"""
    return self.rows.keys()

  def get_num_reactions(self):
    """Return the number of reactions in the kig"""
    return len(self.reactions)

  def print_kig(self):
    """Print to the console the kig"""
    prefix = "    "
    head_line = ""
    for reaction_name in self.reactions:
      head_line += prefix
      head_line += reaction_name
      prefix = ", "
    print (head_line)
    for row_name, row in self.rows.items():
      line = row_name
      for value in row:
        line += ", "
        line += str(value)
      print (line)

class InvariantInferer:
  """A class implementing the invariant inferer over the given 
     kinectic independence graph"""
  def __init__(self, kig):
    self.kig = kig
    self.species = kig.get_species()
    # These will be set when we run 'calculate_invariants'
    self.invariants = []
    self.have_calculated_invariants = False

  class InvRow:
    """A simple class representing an invariant row"""
    def __init__(self, species, row):
      self.species = species 
      self.row = row

    def print_row(self):
      """Print out the invariant row"""
      print (str(self.species) + str(self.row))

     
  def combine_invariant_rows(self, row1, row2):
    """A (private) helper function to combine to rows into a
       single row by summing the columns"""
    species = utils.add_lists(row1.species, row2.species)
    row = utils.add_lists(row1.row, row2.row)
    return self.InvRow(species, row)

  def calculate_initial_rows(self):
    """A (private) function to calculate the initial rows of the
       invariant matrix prior to any computation."""
    inv_rows = []
    species = self.species
    num_species = len(species)
    kig_rows_dict = self.kig.get_rows_dictionary()
    for index in range(num_species):
      inv_species = [0] * num_species
      inv_species[index] = 1
      inv_row = self.InvRow(inv_species, kig_rows_dict[species[index]])
      inv_rows.append(inv_row)

    return inv_rows

  def calculate_invariant_rows(self):
    """Calculate the rows of the invariant matrix"""
    inv_rows = self.calculate_initial_rows()
    # Lets begin by just trying to remove 
    for index in range(self.kig.get_num_reactions()):
      num_rows = len(inv_rows)
      new_inv_rows = []
      for i in range(num_rows):
        i_row   = inv_rows[i]
        i_value = i_row.row[index]
        if i_value == 0:
          new_inv_rows.append(i_row)
        else: 
          for j in range(i+1, num_rows):
            j_row   = inv_rows[j]
            j_value = j_row.row[index]
            sum_value = i_value + j_value 
            if i_value != 0 and sum_value == 0:
              new_row = self.combine_invariant_rows(i_row, j_row)
              new_inv_rows.append(new_row)
      # new_inv_rows = [ r for r in inv_rows if r.row[index] == 0 ]
      inv_rows = new_inv_rows
    return inv_rows

  class Invariant:
    """A simple class representing an invariant"""
    def __init__(self, species, invariant_row):
      self.species = species
      self.coeffs  = invariant_row.species
      assert(len(self.species) == len(self.coeffs))
      self.involved_species = []
      for i in range(len(species)):
        if self.coeffs[i] != 0:
          self.involved_species.append(species[i])

    def format_as_expression (self):
      """Return a string representing the stored invariant as an
         expression summing to the value of the invariant"""
      items = []
      for i in range(len(self.species)):
        coeff = self.coeffs[i]
        if coeff != 0:
          name = self.species[i]
          if coeff == 1:
            items.append(name)
          else:
            items.append("(" + str(coeff) + name + ")")
        
      return " + ".join(items)

  def calculate_invariants(self):
    """Calculate the set of species invariants for this instance
       which should already contain the kig (from the initialiser)"""
    inv_rows = self.calculate_invariant_rows()
    self.invariants = [ self.Invariant(self.species, inv_row) 
                          for inv_row in inv_rows ]
    self.have_calculated_invariants = True
    return self.invariants

  def get_uncovered(self):
    """Return the list of species which are not covered by any
       invariant.
    """
    if not self.have_calculated_invariants :
      self.calculate_invariants()
    unused = self.species[:]
    for invariant in self.invariants:
      for name in invariant.involved_species:
        try:
          unused.remove(name)
        except ValueError:
          pass
    return unused
      

def kig_of_model(model, ignore_sources, ignore_sinks):
  """Compute and return the kinetic independence graph of 
     a model. The kig is simply a mapping from reactions to
     their effects on the population of each species"""
  species = outline_sbml.get_list_of_species(model)
  species_names = [ spec.get_name() for spec in species ]
  all_reactions = outline_sbml.get_list_of_reactions(model)
  reactions = [ r for r in all_reactions 
                    if (not r.is_source() or not ignore_sources) 
                        and 
                       (not r.is_sink() or not ignore_sinks) ]

  reaction_names = [ r.name for r in reactions ]
  kig = KinecticIndependenceGraph(species_names, reaction_names)
  for reaction in reactions:
    kig.add_reaction_info(reaction)
  return kig

def remove_inversed_reactions(all_reactions):
  """Returns the given list of reactions such that for any pair of
     reactions which are the inverse of each other, only one is
     in the returned list of reactions. In other words it removes
     reactions which are the inverse of something else in the list, but
     keeps one of the pair.
  """
  reactions = []
  # So basically for each reaction check if its reverse has already
  # been added to 'reactions' and if not then add it to reactions, but
  # if so, do not. This means we only get one of two reactions which are
  # the inverse of each other.
  for reaction in all_reactions:
    if any(already_in.is_reverse(reaction) for already_in in reactions):
      break
    reactions.append(reaction)
  return reactions

def create_knocked_out_kig(reaction, reactions, species_names):
  """Create a kig from a set of reaction, but ignoring the given
     reaction. 
  """      
  reaction_names = [ r.get_name() for r in reactions if r != reaction ]  
  kig = KinecticIndependenceGraph(species_names, reaction_names)
  for other_reaction in reactions:
    if other_reaction != reaction:
      kig.add_reaction_info(other_reaction)
  return kig


def reaction_knockout_kigs(filename, ignore_sources, ignore_sinks):
  """From an sbml file, return a list of kinetic information graphs,
     where for each reaction in the model there is a single kig which
     is derived from the model whilst ignoring that one reaction
     (and its inverse if present) and the source and/or sink reactions if
     those parameters are set to True
  """
  dom = xml.dom.minidom.parse(filename)
  model = dom.getElementsByTagName("model")[0]
  species = outline_sbml.get_list_of_species(model)
  species_names = [ spec.get_name() for spec in species ]
  all_reactions = outline_sbml.get_list_of_reactions(model)
  non_ignored_reactions = [ r for r in all_reactions 
                              if (not r.is_source() or not ignore_sources) 
                              and 
                              (not r.is_sink() or not ignore_sinks) 
                          ]
  reactions = remove_inversed_reactions(non_ignored_reactions)

  kig_dictionary = dict()
  for reaction in reactions:
    kig = create_knocked_out_kig(reaction, reactions, species_names)
    kig_dictionary[reaction.get_name()] = kig
  return kig_dictionary


# We should also check if there is anything in the model, for example
# rateRules which may violate what is otherwise an invariant within
# the model.
def invariant_inferer_model_file(filename, 
                                 ignore_sources, 
                                 ignore_sinks):
  """Parse in a file as an SBML model, extract the reaction information
     and from this create an invariant inferer.
  """
  dom = xml.dom.minidom.parse(filename)
  model = dom.getElementsByTagName("model")[0]
  kig = kig_of_model(model, ignore_sources, ignore_sinks)
  # kig.print_kig()

  return InvariantInferer(kig)


def get_uncovered_model_file(filename, ignore_sources, ignore_sinks):
  """Return the list of species in a model which are uncovered by
     any invariant
  """
  inv_inferer = invariant_inferer_model_file(filename,
                                             ignore_sources, 
                                             ignore_sinks)
  uncovered  = inv_inferer.get_uncovered()
  return uncovered 
 

def display_invariants_model_file(filename, ignore_sources, ignore_sinks):
  """Given an sbml model file calculate the invariants (whilst ignoring
     source and/or sink reactions according to the arguments) and print
     out those invariants and also some information such as any species
     which are not covered by invariants.
  """
  inv_inferer = invariant_inferer_model_file(filename,
                                             ignore_sources, 
                                             ignore_sinks)
  invariants  = inv_inferer.calculate_invariants()
  if not invariants:
    print ("There are no invariants in this model")
  for invariant in invariants:
    print (invariant.format_as_expression()) 

  uncovered = inv_inferer.get_uncovered()
  if uncovered:
    print ("The following components are not covered by any invariant")
    for name in uncovered:
      print (name)
    if not ignore_sources or not ignore_sinks:
      print ("Try again whilst ignore source and sink reactions using the")
      print ("command-line flags: --ignore-sources --ignore-sinks")
  else:
    print ("There are no components uncovered by invariants")

def reaction_knockout_table(filename, ignore_sources, ignore_sinks):
  """Given an sbml file, perform composite reduction analysis and
     print the results. This basically means perform invariant analysis
     once for each reaction in the model. For a reaction's analysis we
     ignore that reaction and calculate the set of invariants and the
     set of species not covered by any invariant. The idea is to compare
     this with the set of invariants for all reactions thus giving us
     some idea which reactions may be violating an expected invariant
  """
  kig_dictionary = reaction_knockout_kigs(filename,
                                          ignore_sources,
                                          ignore_sinks)

  def get_length_uncovered(my_pair):
    """Return the number of uncovered species"""
    return len(my_pair[1])

  def get_uncovered(kig):
    """Return the uncovered species from a kig"""
    inv_inferer = InvariantInferer(kig)
    return inv_inferer.get_uncovered()

  uncovered_map = [ (rname, get_uncovered(kig) )
                    for (rname, kig) in kig_dictionary.items() ]
  sorted_uncovereds = sorted(uncovered_map, key=get_length_uncovered)
  for (reaction_name, uncovered) in sorted_uncovereds:
    print ("----- " + reaction_name + " -------")
    print (", ".join(uncovered))


def calculate_invariant_timecourse(invariant, timecourse):
  """For a given invariant (that is a set of species and their
     coefficients which the model suggests should sum to a constant
     value of a simulation), and time course data, return the time course
     of the invariant's value according to the time course data.
  """
  columns = [ timecourse.get_column_data(species_name)
                for species_name in invariant.species ]
 
  invariant_timecourse_rows = [] 
  for i in range(0, len(columns[0])):
    invariant_value = 0
    for j in range(0, len(columns)):
      coeff = invariant.coeffs[j]
      value = columns[j][i] 
      invariant_value += value * coeff
    invariant_timecourse_rows.append(invariant_value)
  
  return invariant_timecourse_rows
    
    

def check_experimental_data(sbml_file, timecourse_file):
  """Parse in both an sbml file and a time course data file. Calculate
     the invariants of the sbml file and attempt to determine whether
     or not the experimental data conforms to the invariant. We could also
     get the initial populations of the species and determine whether or
     not the model is capable of matching the experimental data.
  """
  inv_inferer = invariant_inferer_model_file(sbml_file, False, False)
  timecourse = timeseries.get_timecourse_from_file(timecourse_file)

  invariants = inv_inferer.calculate_invariants()
 
  for invariant in invariants:
    invariant_timecourse_rows = calculate_invariant_timecourse(invariant,
                                                               timecourse)
    print(invariant_timecourse_rows) 
  
  
  

def run():
  """Perform the banalities of command-line argument processing
     and then go ahead an calculate the invariants for each given
     SBML model"""
  description = "Analyse SBML files for invariants"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="an sbml file to check invariants for")
  parser.add_argument("--ignore-sources",
                      action="store_true", default=False,
    help="Ignore source reactions when computing invariants")
  parser.add_argument("--ignore-sinks",
                      action="store_true", default=False,
    help="Ignore sink reactions when computing invariants")
  parser.add_argument("--reaction-knockout-table",
                      action="store_true", default=False,
    help="Selectively knock-out each reaction and report the uncovereds")
  parser.add_argument('--check-exp-data', action='store',
                      help="Time course file to check invariants against")

  arguments = parser.parse_args()

  for filename in arguments.filenames:
    if arguments.reaction_knockout_table:
      reaction_knockout_table(filename, 
                              arguments.ignore_sources,
                              arguments.ignore_sinks)
    elif arguments.check_exp_data:
      check_experimental_data(filename, arguments.check_exp_data)
    else :
      display_invariants_model_file(filename,
                                    arguments.ignore_sources,
                                    arguments.ignore_sinks)

if __name__ == "__main__":
  run()
