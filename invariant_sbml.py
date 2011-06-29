"""A script to calculate the invariants of an SBML model"""
import xml.dom.minidom
import outline_sbml
import argparse

def add_lists(left, right):
  """A utility function to add the elements of two lists together
     to form a new list. The result list will be the same length as
     the two given lists which must be of equal length"""
  assert(len(left) == len(right))
  result = []
  for i in range(len(left)):
    result.append(left[i] + right[i])
  return result

def format_list(separator, items):
  """Formats a list of strings as a single string containing all
     items separated by 'separator'"""
  result = ""
  prefix = ""
  for item in items:
    result += prefix
    result += item
    prefix = separator
  return result


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
    reactants = reaction.get_reactants()
    column_index = self.reactions.index(reaction.get_name())
    for reactant in reactants:
      value = - reactant.get_stoichiometry()
      self.rows[reactant.get_name()][column_index] = value
    for product in reaction.get_products():
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
    species = add_lists(row1.species, row2.species)
    row = add_lists(row1.row, row2.row)
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
    """A simple class representing on invariant"""
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
        
      return format_list(" + ", items)
       

  def calculate_invariants(self):
    """Calculate the set of species invariants for this instance
       which should already contain the kig (from the initialiser)"""
    inv_rows = self.calculate_invariant_rows()
    self.invariants = [ self.Invariant(self.species, inv_row) 
                          for inv_row in inv_rows ]
    return self.invariants

  def get_uncovered(self):
    """Return the list of species which are not covered by any
       invariant, relies upon the 'calculate_invariants' method
       having been invoked first"""
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
  reaction_names = [ r.get_name() for r in reactions ]
  kig = KinecticIndependenceGraph(species_names, reaction_names)
  for reaction in reactions:
    kig.add_reaction_info(reaction)
  return kig

# We should also check if there is anything in the model, for example
# rateRules which may violate what is otherwise an invariant within
# the model.
def invariants_model_file(filename, ignore_sources, ignore_sinks):
  """Parse in a file as an SBML model, extract the reaction information
     and then attempt to calculate the invariants for the given model"""
  dom = xml.dom.minidom.parse(filename)
  model = dom.getElementsByTagName("model")[0]
  kig = kig_of_model(model, ignore_sources, ignore_sinks)
  # kig.print_kig()

  inv_inferer = InvariantInferer(kig)
  invariants  = inv_inferer.calculate_invariants()
  for invariant in invariants:
    print (invariant.format_as_expression()) 

  uncovered = inv_inferer.get_uncovered()
  if uncovered:
    print ("The following components are not covered by any invariant")
    for name in uncovered:
      print (name)
    if not ignore_sources or not ignore_sinks:
      print ("Try again whilst ignore source and sink reactions using the")
      print ("command-line flags: --ignore_sources --ignore_sinks")
  else:
    print ("There are no components uncovered by invariants")

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
                      action="store_true", default=False)
  parser.add_argument("--ignore-sinks",
                      action="store_true", default=False)
 
  arguments = parser.parse_args()

  for filename in arguments.filenames:
    invariants_model_file(filename,
                          arguments.ignore_sources,
                          arguments.ignore_sinks)

if __name__ == "__main__":
  run()
