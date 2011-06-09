"""A script which will generate both an xml schema file for the
   configuration, and a python module to parse and store such a valid
   configuration"""

def print_imports(module_file):
  module_file.write("import sys\n")
  module_file.write("import xml.dom.minidom\n\n\n")

def print_exception_definition(module_file):
  except_def = """
class ConfigurationError(Exception):
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)
"""
  module_file.write(except_def)
  module_file.write("\n\n\n")

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
    self.varname = name.lower()

  def has_xml_children(self):
    return False
 
  def parse_method_name(self):
    return "parse_" + self.name
 
  def print_class_definition(self, module_file):
    module_file.write("class " + self.name + ":\n")
    self.print_parse_function(module_file)



## Note where you have more than one sub-element of the same name
## eg. <FFT> <State> ... </State><State> ... </State> </FFT>
## We'll build a list element, where actually there is no
## actual xml element of that kind, but we fake it here.
class ListOfElements(Element):
  # This should also have a maximum and minimum number
  def __init__(self, item):
    Element.__init__(self, item.name)
    self.item = item

  # print parse function isn't defined and it should never
  # be called on this kind of element. ListOfElements exists
  # so that the parent parsing function can decide whether to
  # parse in a list of such elements or a single one.

   
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
                      "(self, " + self.varname + "):\n")
    for child in self.children:
      module_file.write("    " + child.varname + 
                        "_elements = " + self.varname +
                        ".getElementsByTagName(\"" +
                        child.name + "\")\n")
      if isinstance(child, ListOfElements):
        # If you need to get other data from the ListOfElements,
        # for example the max and min number of occurrences then
        # do that before you make child = child.item. 
        child = child.item
        module_file.write("    " + child.varname + "s = []\n")
        module_file.write("    for " + child.varname + "_element in " + 
                          child.varname + "_elements :\n")
        module_file.write("      " + child.varname + " = " +
                          child.name + "()\n")
        module_file.write("      " + child.varname + "." +
                          child.parse_method_name() +
                          "(" + child.varname + "_element)\n")
        module_file.write("      " + child.varname + "s.append(" +
                          child.varname + ")\n")
      else:
        module_file.write("    if len(" + child.varname + 
                                      "_elements) > 1:\n")
        module_file.write("      raise ConfigurationError (\"too many\")\n")
 
        # We should do something in the case that there are no
        # child elements, but we must first find a way to say
        # that the child element may be not there at all
        module_file.write("    elif len(" + child.varname +
                                        "_elements) == 1:\n")
        module_file.write("      self." + child.varname + " = " +
                          child.name + "()\n")
        module_file.write("      self." + child.varname + "." +
                          child.parse_method_name() +
                          "(" + child.varname + "_elements[0])\n")
        module_file.write("    else:\n")
        module_file.write("      self." + child.varname + " = None\n")
     

class EnumerationElement(Element):
  def __init__(self, name, choices):
    Element.__init__(self, name)
    self.choices = choices

  def print_parse_function(self, module_file):
    arg_name = self.varname
    module_file.write("  def " + self.parse_method_name() +
                      "(self, " + arg_name + "):\n")
    module_file.write("    text = " + arg_name + ".firstChild\n")
    module_file.write("    if text not in [ ")
    module_file.write("\"" + self.choices[0] + "\"")
    for i in range(1, len(self.choices)):
      module_file.write(", \"" + self.choices[i] + "\"")
    module_file.write("]:\n")
    # Should give more information than just the text
    module_file.write("      raise ConfigurationError(text)\n")
    module_file.write("    self.value = text\n\n")
   

class LiteralElement(Element):
  def __init__(self, name, simple_type):
    Element.__init__(self, name)
    self.simple_type = simple_type
 
  def get_data_parser(self, argument):
    # Obviously for each kind it could do some more checking
    if self.simple_type in [ "xs:int",
                             "nonZeroMultipleOfTwo",
                             "nonZeroMultipleOfTwoAndGtThree",
                             "AzusaPercentage",
                             "positiveIntGTOne",
                             "positiveIntNotZero",
                             "positiveIntIncludingZero",
                           ]:
      return "int(" + argument + ")"
    elif self.simple_type in [ "positiveDecimalIncludingZero",
                               "positiveDecimalNotZeroButCanBeMinusOne",
                               "positiveDecimalIncludingZero",
                               "positiveDecimalNotZero",
                               "probability",
                               "minusOne",
                               "positiveFloatNotZero",
                               "positiveFloatIncludingZero",
                             ]:
      return "float(" + argument + ")"
    elif self.simple_type in [ "nonEmptyString" ]:
      # The argument should already be a string so no parsing
      # function required
      return argument
    else:
      print ("Undefined simple type: " + self.simple_type)
      sys.exit(1)

  def print_parse_function(self, module_file):
    arg_name = self.varname
    module_file.write("  def " + self.parse_method_name() +
                      "(self, " + arg_name + "):\n")
    module_file.write("    text = " + arg_name + ".firstChild\n")
    # For most of them it will be float, but actually we will
    # have to check the simple type for some.
    module_file.write("    self.value = " + 
                      self.get_data_parser("text.data") + "\n")
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
                      "(self, " + self.varname + "):\n")
    # Not quite right, we need to find out how many children we have.
    module_file.write("    # num_children = len(" + self.varname +
                      ".childNodes)\n")
    module_file.write("    # if num_children != 1:\n")
    module_file.write("    #   print(\"Configuration parse error\")\n")
    module_file.write("    #   print(\"Number of children = \" +\n")
    module_file.write("    #           str(num_children))\n")
    module_file.write("    #   sys.exit(1)\n")
    module_file.write("    child_name = \"\"\n")
    for choice in self.choices:
      module_file.write("    " + choice.varname + 
                        "s = " + self.varname +
                        ".getElementsByTagName(\"" +
                        choice.name + "\")\n")
      module_file.write("    for " + choice.varname + "_element in " + 
                        choice.varname + "s :\n")
      module_file.write("      " + choice.varname + " = " +
                        choice.name  + "()\n")
      module_file.write("      " + choice.varname + "." +

                        choice.parse_method_name() +
                        "(" + choice.varname + "_element)\n")


def simple_element(name, el_type, dictionary):
  if name in dictionary:
    element = dictionary[name]
    if el_type == element.simple_type:
      return element
    else:
      print ("Conflicting type definitions for: " + name)
      print (el_type)
      print (element.simple_type)
      sys.exit(1)
  # If this is the first mention of 'name' then fine just add it
  element = LiteralElement(name, el_type)
  dictionary[name] = element
  return element

def get_cvode_solver_element(element_dictionary): 
  cvode_children = [ simple_element("TFinal", 
                             "positiveDecimalNotZeroButCanBeMinusOne",
                             element_dictionary),
                     simple_element("TInit", "positiveFloatNotZero",
                                    element_dictionary),
                     simple_element("Interval", "positiveFloatNotZero",
                                    element_dictionary),
                     simple_element("MaxTimes", "positiveIntNotZero", 
                                    element_dictionary),
                     simple_element("Atol", "positiveFloatNotZero",
                                    element_dictionary),
                     simple_element("Reltol", "positiveDecimalNotZero",
                                    element_dictionary),
                   ]

  cvode_element = ContainerElement("CVODESolver", cvode_children)
  element_dictionary[cvode_element.name] = cvode_element
  return cvode_element

def get_model_cost_function(element_dictionary):
  # dcs = data set config
  dsc_children = [ simple_element("FileName", "nonEmptyString",
                                  element_dictionary)
                 ]
  dsc_element = ContainerElement("DataSetConfig", dsc_children)
  element_dictionary[dsc_element.name] = dsc_element
  x2_children = [ simple_element("Weight", "positiveFloatNotZero",
                                 element_dictionary),
                  dsc_element
                ]
  x2_element = ContainerElement("X2Cost", x2_children)
  element_dictionary[x2_element.name] = x2_element


  mcf_choices = [ x2_element ]
  mcf_element = OneOfElement("ModelCostFunction", mcf_choices)
  element_dictionary[mcf_element.name] = mcf_element

  return mcf_element

def get_setup_element(element_dictionary):
  setup_type_choices = [ "SobolSelect", "SobolUnselect", "ReadinFile" ]
  setup_type_element = EnumerationElement("SetupType", setup_type_choices)
  element_dictionary[setup_type_element.name] = setup_type_element
 
  mcf_element = get_model_cost_function(element_dictionary)
  mcf_list =  ListOfElements(mcf_element)

  setup_children = [ simple_element("MaxCost", "positiveIntNotZero",
                                    element_dictionary),
                     simple_element("MaxRandomNumberGenerator", 
                                    "positiveIntNotZero",
                                    element_dictionary),
                     setup_type_element,
                     mcf_list,
                     simple_element("CkpointNum", "xs:int",
                                    element_dictionary),
                     simple_element("ParameterInitialFile",
                                    "nonEmptyString",
                                    element_dictionary),
                   ]


  setup_element = ContainerElement("Setup", setup_children)
  element_dictionary[setup_element.name] = setup_element
  return setup_element

#        <xs:element ref="sbsi_ofwconfig:MaxCost" />
#        <xs:element ref="sbsi_ofwconfig:MaxRandomNumberGenerator" />
#        <xs:element ref="sbsi_ofwconfig:SetupType" />
#        <xs:element ref="sbsi_ofwconfig:CkpointNum" minOccurs="0" />
#        <xs:element ref="sbsi_ofwconfig:ModelCostFunction" maxOccurs="unbounded"/>
#        <xs:element ref="sbsi_ofwconfig:ParameterInitialFile" />
 

def get_element_dictionary():
  """we build up a list of elements for which class definitions
     must be output. this avoids us outputting two class definitions
     for the same element because it is a child of more than one
     other element. Note that there is currently no checking that
     we have not already added an element with the same name, but
     there should be"""
  # we keep a mapping from names to 
  element_dictionary = dict()
  cvode_element = get_cvode_solver_element(element_dictionary)

  solver_choices = [ cvode_element ]
  solver_element = OneOfElement("Solver", solver_choices)
  element_dictionary[solver_element.name] = solver_element

  setup_element = get_setup_element(element_dictionary)

  config_children = [ solver_element, setup_element ]
  config_element = ContainerElement("OFWConfig", config_children)
  element_dictionary[config_element.name] = config_element
  return element_dictionary
  

def print_parse_config(module_file):
  parse_config = """def parse_config_file(filename):
  dom = xml.dom.minidom.parse(filename)
  config_element = dom.getElementsByTagName("OFWConfig")[0]
  ofwconfig = OFWConfig()
  ofwconfig.parse_OFWConfig(config_element)
  return ofwconfig
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
  print_exception_definition(module_file)

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
