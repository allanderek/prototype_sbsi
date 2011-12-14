"""A simple module to put some simple utility functions in"""
import os.path
import argparse

def change_filename_ext(filename, new_ext):
  """Returns a new file name based on the first but with the extension
     replaced by the given extension. The new_ext should include the
     '.' separator if that is desired"""
  basename = os.path.splitext(filename)[0]
  return basename + new_ext
 
def get_new_directory(desired_name):
  """Returns a newly created directory with the desired name. If that
     directory already exists it appends a number on to the name"""
  dir_name = desired_name
  number = 0
  while os.path.exists(dir_name):
    number += 1
    dir_name = desired_name + "_" + str(number)  
    
  # create directory and return path
  os.makedirs(dir_name)
  return dir_name

def has_extension(filename, extensions):
  """Returns true if the given filename has one of the given
     list of extensions. The extensions should include the .
     to begin them, eg [ ".xml" ]
  """
  file_extension = os.path.splitext(filename)[1]
  return file_extension in extensions


def has_xml_or_sbml_ext(filename):
  """Returns true if we believe the file to be an SBML file
     based on the file's extension"""
  return has_extension(filename, [ ".xml", ".sbml" ])

def has_copasi_ext(filename):
  """Returns true if we believe the file to be a copasi model file
     based on the file's extension"""
  return has_extension(filename, [ ".cps" ])


def get_non_ignored (all_names, in_names, ignored_names):
  """Many of these utilities have a --column and a --mcolumn flag.
     The semantics of this is supposed to be that we plot/output/use
     all names specified under --column and ignore those specified
     under --mcolumn. But the default for --column is ALL names and
     the default for --mcolumn is the empty set. This function given
     the list of all names and the list of --column and --mcolumn names
     does the right thing and returns the correct set."""
  if not in_names and not ignored_names:
    return all_names
  elif not in_names and ignored_names:
    return [ n for n in all_names if n not in ignored_names ]
  elif in_names and not ignored_names:
    return [ n for n in all_names if n in in_names ]
  elif in_names and ignored_names:
    return [ n for n in all_names 
                 if n in in_names and n not in ignored_names
           ]


def add_lists(left, right, left_scale=1, right_scale=1):
  """A utility function to add the elements of two lists together
     to form a new list. The result list will be the same length as
     the two given lists which must be of equal length"""
  assert(len(left) == len(right))
  result = []
  for i in range(len(left)):
    left_value = left_scale * left[i]
    right_value = right_scale * right[i]
    result.append(left_value + right_value)
  return result

def equal_lists(left, right):
  """Checks if two lists are 'equal', equal if they are considered
     to be sets
  """
  # Could check the lengths here, but I think we want
  # 'l,l,r' to be equal to 'l,r', so we're ignoring duplicates.
  # Not sure if that's correct to do so though.
  for l_item in left:
    if l_item not in right:
      return False
  for r_item in right :
    if r_item not in left :
      return False
  return True


# the function to calculate the GCD
def gcd(num1, num2):
    if num1 > num2:
        for i in range(1,num2+1):
            if num2 % i == 0:
                if num1 % i == 0:
                    result = i
        return result

    elif num2 > num1:
        for i in range(1,num1+1):
            if num1 % i == 0:
                if num2 % i == 0:
                    result = i
        return result

    else:
        result = num1*num2/num1
        return result

# the function to calculate the LCM
def lcm(num1, num2):
    result = num1*num2/gcd(num1,num2)
    return result



class ListArgumentAction(argparse.Action):
  """ A simple argument action to allow a list of arguments to be
      separated by a comma. This allows slightly more concise command
      lines. For example --column p1 --column p2 --column p3 can be
      reduced to --column p1,p2,p3
      To be used as: 
  parser.add_argument('--column', action=utils.ListArgumentAction,
                      help="Specify a column to be plotted")
  """
  def __call__(self, parser, namespace, values, option_string=None):
    """The action called when the parser encounters the associated
       command-line argument.
    """
    # First of all get the list of values specified as this argument
    list_values = values.split(",")
    # Then get the current list already specified, this allows for
    # example: --column p1,p2,p3 --column p4,p5,p6 or more likely just
    # that some user does not know about the comma separator and does
    # --column p1 --column p2 --column p3 ... etc
    current = getattr(namespace, self.dest)
    if current:
      # If there are values there already extend the list
      current.extend(list_values)
    else: 
      # otherwise set the list to the current set of values
      setattr(namespace, self.dest, list_values)


