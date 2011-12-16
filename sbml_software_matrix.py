"""A script to analyse the sbml software matrix"""
import argparse
from BeautifulSoup import BeautifulSoup

class SoftwarePackage(object):
  """A class for representing a single software package"""
  def __init__(self):
    self.name = None
    self.recent_contact = False
    self.creation = False
    self.simulation = False
    self.analysis = False
    self.database = False
    self.utility = False
 
    self.ode = False
    self.dae = False
    self.pde = False
    self.stochastic = False
    self.events = False
    self.logical = False
    self.other = False

    self.api = None
    self.dep = None
  
    self.linux = False
    self.mac = False
    self.windows = False
    self.web_browser = False

    self.sbml_import = False
    self.sbml_export = False

    self.open_source = False
    self.academic_use = False
    self.commericial = False

def get_percentage(value_1, value_2):
  """Return x compared to y as a percentage, 3,5 = 60%"""
  return int(100.0 * (float(value_1) / float(value_2)))

def boolean_column(column):
  """Returns true if the given column has some content indicating that
     the given cell of the matrix is 'true'
  """
  return bool(column.findAll("span"))
   
def compare_attributes(packages, bool_f1, bool_f2):
  """Returns two numbers the first is the number of the given packages
     for which the first function returns true, the second is the
     number of given packages for which both are true. Useful for working
     out the percentage of packages which satisfy one condition also
     satisfy the second condition
  """
  value_1 = 0
  value_2 = 0
  for package in packages:
    if bool_f1(package):
      value_1 += 1
      if bool_f2(package):
        value_2 += 1

  return (value_1, value_2)

def return_number(packages, criterion_f):
  """Return the number of packages which satisfy a given criteria
     implemented by the given function
  """
  number = 0
  for package in packages:
    if criterion_f(package):
      number += 1
  return number


def parse_software_matrix(filename):
  """Parse in the software matrix file and return a list of
     SoftwarePackage's
  """
  page = open (filename, "r")

  soup = BeautifulSoup(page)

  tables = soup("table")

  software_table = tables[1]

  rows = software_table("tr")
  software_packages = []
  for row in rows:
    columns = row("td")
    software_package = SoftwarePackage()
    software_package.name = columns[0].a.string

    software_package.recent_contact = boolean_column(columns[1])
    software_package.creation = boolean_column(columns[2])
    
    software_package.creation = boolean_column(columns[2])
    software_package.simulation = boolean_column(columns[3])
    software_package.analysis = boolean_column(columns[4])
    software_package.database = boolean_column(columns[5])
    software_package.utility = boolean_column(columns[6])

    software_package.ode = boolean_column(columns[7])
    software_package.dae = boolean_column(columns[8])
    software_package.pde = boolean_column(columns[9])
    software_package.stochastic = boolean_column(columns[10])
    software_package.events = boolean_column(columns[11])
    software_package.logical = boolean_column(columns[12])
    software_package.other = boolean_column(columns[13])

    software_package.api = boolean_column(columns[14])
    software_package.dep = boolean_column(columns[15])

    software_package.linux = boolean_column(columns[16])
    software_package.mac = boolean_column(columns[17])
    software_package.windows = boolean_column(columns[18])
    software_package.web_browser = boolean_column(columns[19])

    software_package.sbml_import = boolean_column(columns[20])
    software_package.sbml_export = boolean_column(columns[21])

    software_package.open_source = boolean_column(columns[22])
    software_package.academic_use = boolean_column(columns[23])
    software_package.commericial = boolean_column(columns[24])
 
    software_packages.append(software_package)
  page.close()
  return software_packages


def analyse_file(filename):
  """Analyse the given software matrix file"""
  software_packages = parse_software_matrix(filename)
  contact_packages = [ s for s in software_packages if s.recent_contact ]

  def report_recent_contact(boolean_f, desc):
    """Report on the number/percentage of packages which have a given
       feature tested for by 'boolean_f', which also have recent contact
    """
    in_general = return_number(software_packages, boolean_f)
    of_have_contact = return_number(contact_packages, boolean_f)
    percentage = get_percentage(of_have_contact, in_general)
    print(str(percentage) + "% of " + desc + " have recent contact")


  print ("----Cabilities-----")
  report_recent_contact(lambda x: x.creation, "creation")
  report_recent_contact(lambda x: x.simulation, "simulation")
  report_recent_contact(lambda x: x.analysis, "analysis")
  report_recent_contact(lambda x: x.database, "database")
  report_recent_contact(lambda x: x.utility, "utility")
  def has_no_capabilities(pack):
    """Returns true if the package has no SBML capabilities"""
    return not (pack.creation or pack.simulation or pack.analysis
                or pack.database or pack.utility)
  report_recent_contact(has_no_capabilities, "no capabilities")

  print ("----Frameworks-----")
  report_recent_contact(lambda x: x.ode, "ode")
  report_recent_contact(lambda x: x.dae, "dae")
  report_recent_contact(lambda x: x.pde, "pde")
  report_recent_contact(lambda x: x.stochastic, "stochastic")
  report_recent_contact(lambda x: x.events, "events")
  report_recent_contact(lambda x: x.logical, "logical")
  report_recent_contact(lambda x: x.other, "other")
  def has_no_framework(pack):
    """Returns true if there is no specified frame work for the package"""
    return not (pack.ode or pack.dae or pack.pde or pack.stochastic or
                pack.events or pack.logical or pack.other)
  report_recent_contact(has_no_framework, "no framework")

  print (" ---- Platform ---- ")
  report_recent_contact(lambda x: x.mac, "mac available")
  report_recent_contact(lambda x: x.linux, "linux available")
  report_recent_contact(lambda x: x.windows, "windows available")
  report_recent_contact(lambda x: x.web_browser, "web available")
  def has_no_platform(package):
    """Returns true if the package has no specified platform"""
    return not(package.linux or package.mac or 
               package.windows or package.web_browser)
  report_recent_contact(has_no_platform, "no platform available")

  def has_multiple_platforms(package):
    """returns true if the package is available for multiple platforms"""
    number = 0
    if package.linux:
      number += 1
    if package.mac:
      number += 1
    if package.windows:
      number += 1
    if package.web_browser:
      number += 1
    return number > 1
  report_recent_contact(has_multiple_platforms, "multi-platform")


 
def run():
  """ The main program for the sbml software matrix analyser
  """
  description = "Analyse the sbml software matrix"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="The software matrix file")

  arguments = parser.parse_args()

  filenames = arguments.filenames
  if not filenames:
    filenames = [ "SBMLMatrix-static-matrix-pure.html" ]

  for filename in arguments.filenames:
    # Since it always converts to UserModel/UserModel.{ch}
    # more than one file will just overwrite the others, so
    # perhaps we need to think about this somewhat more.
    analyse_file(filename)

if __name__ == "__main__":
  run()
