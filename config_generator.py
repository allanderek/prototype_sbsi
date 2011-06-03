"""A script which will generate both an xml schema file for the
   configuration, and a python module to parse and store such a valid
   configuration"""

def print_imports(module_file):
  module_file.write("import sys\n")
  module_file.write("import xml.dom.minidom\n\n\n")


def print_numerical_types(schema_file):
  numerical_schema_types = """
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:sbsi_ofwconfig="http://www.uk.ac.ed.csbe.sbsivisual/ofwconfig/ns/level1" 
	 targetNamespace="http://www.uk.ac.ed.csbe.sbsivisual/ofwconfig/ns/level1" elementFormDefault="qualified">

<xs:simpleType name="nonZeroMultipleOfTwo"  xmlns:sbsi_ofwconfig="http://www.uk.ac.ed.csbe.sbsivisual/ofwconfig/ns/level1">
		<xs:restriction base="xs:positiveInteger">
			<xs:pattern value=".*2" />
			<xs:pattern value=".*4" />
			<xs:pattern value=".*6" />
			<xs:pattern value=".*8" />
			<xs:pattern value=".*0" />
		</xs:restriction>
	</xs:simpleType>
	
	<xs:simpleType name="nonZeroMultipleOfTwoAndGtThree" xmlns:sbsi_ofwconfig="http://www.uk.ac.ed.csbe.sbsivisual/ofwconfig/ns/level1">
		<xs:restriction base="xs:positiveInteger">
			<xs:pattern value=".*2" />
			<xs:pattern value=".*4" />
			<xs:pattern value=".*6" />
			<xs:pattern value=".*8" />
			<xs:pattern value=".*0" />
			<xs:minInclusive value="4" />
		</xs:restriction>
	</xs:simpleType>
	
	<xs:simpleType name="AzusaPercentage" xmlns:sbsi_ofwconfig="http://www.uk.ac.ed.csbe.sbsivisual/ofwconfig/ns/level1">
		<xs:restriction base="xs:int">
			<xs:minInclusive value="0" />
			<xs:maxInclusive value="100" />
		</xs:restriction>
	</xs:simpleType>
	
	<xs:simpleType name="probability" xmlns:sbsi_ofwconfig="http://www.uk.ac.ed.csbe.sbsivisual/ofwconfig/ns/level1">
		<xs:restriction base="xs:float">
			<xs:maxInclusive value="1" />
		</xs:restriction>
	</xs:simpleType>

	<xs:simpleType name="positiveDecimalIncludingZero" xmlns:sbsi_ofwconfig="http://www.uk.ac.ed.csbe.sbsivisual/ofwconfig/ns/level1">
		<xs:restriction base="xs:double">
			<xs:minInclusive value="0" />
		</xs:restriction>
	</xs:simpleType>

	<xs:simpleType name="positiveDecimalNotZero" xmlns:sbsi_ofwconfig="http://www.uk.ac.ed.csbe.sbsivisual/ofwconfig/ns/level1">
		<xs:restriction base="xs:double">
			<xs:minExclusive value="0" />
		</xs:restriction>
	</xs:simpleType>

<xs:simpleType name="minusOne" xmlns:sbsi_ofwconfig="http://www.uk.ac.ed.csbe.sbsivisual/ofwconfig/ns/level1">
 		<xs:restriction base="xs:double">
 			<xs:enumeration value="-1"></xs:enumeration>
 		</xs:restriction>
 	</xs:simpleType>
 	
 	<xs:simpleType name="positiveDecimalNotZeroButCanBeMinusOne" xmlns:sbsi_ofwconfig="http://www.uk.ac.ed.csbe.sbsivisual/ofwconfig/ns/level1">
 		<xs:union memberTypes="sbsi_ofwconfig:positiveDecimalNotZero sbsi_ofwconfig:minusOne"/>
 	</xs:simpleType>

	
	<xs:simpleType name="positiveFloatNotZero" xmlns:sbsi_ofwconfig="http://www.uk.ac.ed.csbe.sbsivisual/ofwconfig/ns/level1">
			<xs:restriction base="xs:float">
				<xs:minExclusive value="0" />
			</xs:restriction>
		</xs:simpleType>

	<xs:simpleType name="positiveFloatIncludingZero" xmlns:sbsi_ofwconfig="http://www.uk.ac.ed.csbe.sbsivisual/ofwconfig/ns/level1">
			<xs:restriction base="xs:float">
				<xs:minInclusive value="0" />
			</xs:restriction>
		</xs:simpleType>
	
	<xs:simpleType  name="positiveIntGTOne" xmlns:sbsi_ofwconfig="http://www.uk.ac.ed.csbe.sbsivisual/ofwconfig/ns/level1">
		<xs:restriction base="xs:int">
		<xs:minExclusive value="1" />
		</xs:restriction>
	</xs:simpleType>
	
	<xs:simpleType  name="positiveIntNotZero" xmlns:sbsi_ofwconfig="http://www.uk.ac.ed.csbe.sbsivisual/ofwconfig/ns/level1">
		<xs:restriction base="xs:int">
		<xs:minExclusive value="0" />
		</xs:restriction>
	</xs:simpleType>
	
	<xs:simpleType  name="positiveIntIncludingZero" xmlns:sbsi_ofwconfig="http://www.uk.ac.ed.csbe.sbsivisual/ofwconfig/ns/level1">
		<xs:restriction base="xs:int">
		<xs:minInclusive value="0" />
		</xs:restriction>
	</xs:simpleType>
	
	<xs:simpleType name="nonEmptyString" xmlns:sbsi_ofwconfig="http://www.uk.ac.ed.csbe.sbsivisual/ofwconfig/ns/level1">
      	<xs:restriction base="xs:string">
      	<xs:minLength value="1"/>
    	</xs:restriction>
    </xs:simpleType>
	
  """
  schema_file.write(numerical_schema_types)


class Element:
  def __init__(self, name):
    self.name = name

  def has_xml_children(self):
    return False
 
  def parse_method_name(self):
    return "parse_" + self.name
 
  def print_class_definition(self, module_file):
    module_file.write("class " + self.name + ":\n")
    self.print_parse_function(module_file)
    
class ContainerElement(Element):
  
  def __init__(self, name, children):
    Element.__init__(self, name)
    self.children = children


  def has_xml_children(self):
    return self.children != []

  def add_child(self, child):
    self.children.append(child)

  def print_parse_function(self, module_file):
    module_file.write("  def " + self.parse_method_name() +
                      "(self, " + self.name + "):\n")
    for child in self.children:
      module_file.write("    " + child.name + 
                        "s = " + self.name +
                        ".getElementsByTagName(\"" +
                        child.name + "\")\n")
      module_file.write("    for " + child.name + "_element in " + 
                        child.name + "s :\n")
      module_file.write("      self." + child.name + " = " +
                        child.name + "()\n")
      module_file.write("      " + child.name + "." +
                        child.parse_method_name() +
                        "(" + child.name + "_element)\n")
##
## Note where you have more than one sub-element of the same name
## eg. <FFT> <State> ... </State><State> ... </State> </FFT>
## We'll build a list element, where actually there is no
## actual xml element of that kind, but we fake it here.
class ListOfElements(Element):
  def __init__(self, name, item):
    Element.__init__(self, name)
    self.item = item

  # So the parse function here, unlike the others will be passed
  # in the parent element rather than the actual element.


class EnumerationElement(Element):
  def __init__(self, name, choices):
    Element.__init__(self, name)
    self.choices = choices

class LiteralElement(Element):
  def __init__(self, name, simple_type):
    Element.__init__(self, name)
    self.simple_type = simple_type
 
  def print_parse_function(self, module_file):
    arg_name = self.name
    module_file.write("  def " + self.parse_method_name() +
                      "(self, " + arg_name + "):\n")
    module_file.write("    text = " + arg_name + ".firstChild\n")
    # For most of them it will be float, but actually we will
    # have to check the simple type for some.
    module_file.write("    self.value = float(text)\n")
    # We should also have some code to check for the actual constraints
    # imposed by the particular simpleType.
    

class OneOfElement(Element):
  def __init__(self, name, choices):
    Element.__init__(self, name)
    self.choices = choices

  # Hmm, it is something similar to this, but not quite the
  # same as this.
  def print_parse_function(self, module_file):
    module_file.write("  def " + self.parse_method_name() +
                      "(" + self.name + "):\n")
    # Not quite right, we need to find out how many children we have.
    module_file.write("    num_children = 0\n")
    module_file.write("    if num_children != 1:\n")
    module_file.write("      print(\"Configuration parse error\")\n")
    module_file.write("      sys.exit(1)\n")
    module_file.write("    child_name = \"\"\n")
    for choice in self.choices:
      module_file.write("    " + choice.name + 
                        "s = " + self.name +
                        ".getElementsByTagName(\"" +
                        choice.name + "\")\n")
      module_file.write("    for " + choice.name + " in " + 
                        choice.name + "s :\n")
      module_file.write("      self." +
               
                        choice.parse_method_name() +
                        "(" + choice.name + ")\n")
   


      
def get_element_dictionary():
  """we build up a list of elements for which class definitions
     must be output. this avoids us outputting two class definitions
     for the same element because it is a child of more than one
     other element. Note that there is currently no checking that
     we have not already added an element with the same name, but
     there should be"""
  # we keep a mapping from names to 
  element_dictionary = dict()
  # This function appears to be being written from the bottom up
  # I didn't really mean this to be the case, but ah well dem's
  # the breaks.
  t_final_simple_type = "positiveDecimalNotZeroButCanBeMinusOne"
  t_final_element = LiteralElement("TFinal", t_final_simple_type) 
  element_dictionary[t_final_element.name] = t_final_element

  cvode_children = [t_final_element]
  cvode_element = ContainerElement("CVODESolver", cvode_children)
  element_dictionary[cvode_element.name] = cvode_element

  solver_choices = [ cvode_element ]
  solver_element = OneOfElement("Solver", solver_choices)
  element_dictionary[solver_element.name] = solver_element


  config_children = [ solver_element ]
  config_element = ContainerElement("OFWConfig", config_children)
  element_dictionary[config_element.name] = config_element
  return element_dictionary
  

def print_parse_config(module_file):
  parse_config = """def parse_config_file(filename):
  dom = xml.dom.minidom.parse(filename)
  config_element = dom.getElementsByTagName("OFWConfig")[0]
  ofwconfig = OFWConfig()
  ofwconfig.parse_OFWConfig(config_element)
  """ 
  module_file.write(parse_config)

def print_run_def(module_file):
  run_def = """
def run():
  arguments = [ x for x in sys.argv if not x.endswith(".py") ]
  for filename in arguments:
    parse_config_file(filename)
"""
  module_file.write(run_def)

def run():
  schema_file = open ("ofwconfig.xsd", "w")
  module_file = open ("configuration.py", "w")

  print_imports(module_file)

  print_numerical_types(schema_file)

  element_dictionary = get_element_dictionary()
  for element in element_dictionary.values():
    element.print_class_definition(module_file)

  print_parse_config(module_file) 

  print_run_def(module_file)
  main_file_run = "if __name__ == \"__main__\":\n  run()"
  module_file.write(main_file_run)

  schema_file.write("\n")
  schema_file.close()
  module_file.write("\n")
  module_file.close()
  

if __name__ == "__main__":
  run()
