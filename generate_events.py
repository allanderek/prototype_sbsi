""" A simple python script which should, given time course data,
    produce a series of events to be added to an SBML model, such
    that the events simulate the changing of the population as in
    the time course data."""
import os
import sys
import argparse
import xml.dom.minidom
import timeseries

class EventAssign:
  """A simple class representing an event assignment"""
  def __init__(self, species, value):
    self.species = species
    self.value = value

class Event:
  """This is actually a poor-man's event class. That is why
     it is here, because it is specific to this use case. Ideally
     we should have a more general Event class and we should just
     use that here."""
  def __init__(self, time, event_assigns):
    self.time = time
    self.event_assigns = event_assigns
    self.math_ns = "http://www.w3.org/1998/Math/MathML"

  def format_event(self):
    """A function to return the event as a formatted xml string"""
    # You need to create the document to which you're going to
    # create elements within. 
    document = xml.dom.minidom.Document()
    element = self.create_element(document)
    formatted = element.toprettyxml(indent="  ", encoding="UTF-8")
    return formatted

  def create_event_assigments(self, document):
    """Create and return an xml element within the given document
       which represents the list of event assignments corresponding
       to those event assignments for this event"""
    #loea = listOfEventAssignments
    loea = document.createElement("listOfEventAssignments")
    for event_assign in self.event_assigns:
      e_assign = document.createElement("eventAssignment")
      loea.appendChild(e_assign)
      e_assign.setAttribute("variable", event_assign.species)
      e_math = document.createElementNS(self.math_ns, "math")
      e_assign.appendChild(e_math)
      e_math.setAttribute("xmlns", self.math_ns)
      e_cn = document.createElement("cn")
      e_math.appendChild(e_cn)
      e_cn_text = document.createTextNode(str(event_assign.value)) 
      e_cn.appendChild(e_cn_text)
    return loea
 
  def create_element(self, document):
    """Created an return the xml element representing this event"""
    event = document.createElement("event")
    trigger = document.createElement("trigger")
    event.appendChild(trigger)
    trigger_math = document.createElement("math")
    trigger.appendChild(trigger_math)
    trigger_math.setAttribute("xmlns", self.math_ns)
    trigger_apply = document.createElement("apply")
    trigger_math.appendChild(trigger_apply)
    trigger_gt = document.createElement("gt")
    trigger_apply.appendChild(trigger_gt)
    trigger_csymbol = document.createElement("csymbol")
    trigger_apply.appendChild(trigger_csymbol)
    trigger_csymbol.setAttribute("encoding", "text")
    time_def_url = "http://www.sbml.org/sbml/symbols/time"
    trigger_csymbol.setAttribute("definitionURL", time_def_url)
    t_text = document.createTextNode("t")
    trigger_csymbol.appendChild(t_text)
    trigger_cn = document.createElement("cn")
    trigger_apply.appendChild(trigger_cn)
    trigger_cntext = document.createTextNode(str(self.time))
    trigger_cn.appendChild(trigger_cntext)
    
     
    #loea = listOfEventAssignments
    loea = self.create_event_assigments(document)
    event.appendChild(loea)
    return event
    
  def create_copasi_element(self, document, event_number):
    """Create and return a copasi xml element representing this event"""
    event = document.createElement("Event")
    event_name = "TimeCourseEvent_" + str(event_number)
    event.setAttribute("key", event_name)
    event.setAttribute("name", event_name)
    event.setAttribute("order", str(event_number + 1))
    trigger = document.createElement("TriggerExpression")
    trigger_text = ("<CN=Root,Model=NoName,Reference=Time> gt " +
                    str(self.time))
    text_node = document.createTextNode(trigger_text)
    trigger.appendChild(text_node)
    event.appendChild(trigger)

    list_of_assignments = document.createElement("ListOfAssignments")
    for event_assign in self.event_assigns:
      assignment = document.createElement("Assignment")
      assignment.setAttribute("targetKey", event_assign.species)
      expression = document.createElement("Expression")
      text_node = document.createTextNode(str(event_assign.value))
      expression.appendChild(text_node)
      assignment.appendChild(expression)
      list_of_assignments.appendChild(assignment)
    
    event.appendChild(list_of_assignments)
    return event

      # <Event key="Event_0" name="event_0" order="1">
      #   <TriggerExpression>
      #     &lt;CN=Root,Model=NoName,Reference=Time&gt; gt 0.1
      #   </TriggerExpression>
      #   <ListOfAssignments>
      #     <Assignment targetKey="ModelValue_2">
      #       <Expression>
      #         18.6355
      #       </Expression>
      #     </Assignment>
      #     <Assignment targetKey="ModelValue_6">
      #       <Expression>
      #         10
      #       </Expression>
      #     </Assignment>
      #   </ListOfAssignments>
      # </Event>



def events_from_timecourse(timecourse):
  """Create a list of events given a timecourse"""
  events = []
  column_names = timecourse.get_column_names()
  for row in timecourse.get_rows():
    time = row[0]
    event_assigns = []
    for i in range(len(column_names)):
      name = column_names[i]
      # Because the column_names don't include the time value
      # hence the index into the row will be one greater than
      # the index into the column_names list.
      value = row[i + 1]
      event_assign = EventAssign(name, value)
      event_assigns.append(event_assign)
    event = Event(time, event_assigns)
    events.append(event)
  return events
  


 
def has_xml_or_sbml_ext(filename):
  """Returns true if we believe the file to be an SBML file
     based on the file's extension"""
  extension = os.path.splitext(filename)[1]
  return extension == ".xml" or extension == ".sbml"

def has_copasi_ext(filename):
  """Returns true if we believe the file to be a copasi model file
     based on the file's extension"""
  extension = os.path.splitext(filename)[1]
  return extension == ".cps"

def check_constant_false(model, species):
  """Checks that if the given species is defined in the model as a
     parameter, then it has its constant attribute set to false.
     This is because any variable which has its value updated
     (for example by an event that we've added) must have its constant
     attribute set to false"""
  list_of_parameters_lists = model.getElementsByTagName("listOfParameters")
  if list_of_parameters_lists:
    for list_of_parameters in list_of_parameters_lists:
      parameters = list_of_parameters.getElementsByTagName("parameter")
      for parameter in parameters:
        name = parameter.getAttribute("id")
        if name in species:
          parameter.setAttribute("constant", "false")

def add_events_sbml(filename, events, species, arguments):
  """Parse in a file as an SBML model, add the given events and then
     print out the augmented model"""
  dom = xml.dom.minidom.parse(filename)
  model = dom.getElementsByTagName("model")[0]
  
  # Where loe = listOfEvents
  loe_elements = model.getElementsByTagName("ListOfEvents")
  if not loe_elements:
    loe_element = dom.createElement("ListOfEvents")

    state_template = model.getElementsByTagName("StateTemplate") 
    model.insertBefore(state_template)
  else:
    loe_element = loe_elements[0]

  for event in events:
    loe_element.appendChild(event.create_element(dom)) 

  check_constant_false (model, species)
 
  if arguments.pretty:
    document = dom.toprettyxml(indent="  ", encoding="UTF-8")
  else:
    document = dom.toxml("UTF-8")
  print(document)

def add_events_copasi(filename, events, species, arguments):
  """Parse in a file as a copasi model, add the given events and then
     print out the augmented model"""
  dom = xml.dom.minidom.parse(filename)
  model = dom.getElementsByTagName("COPASI")[0]
  
  # Where loe = listOfEvents
  loe_elements = model.getElementsByTagName("ListOfEvents")
  if not loe_elements:
    loe_element = dom.createElement("ListOfEvents")
    # listOfEvents is the last child of a model so that we can
    # just append it regardless of what other child nodes are in
    # this particular model.
    model.appendChild(loe_element)
  else:
    loe_element = loe_elements[0]

  event_number = 0
  for event in events:
    loe_element.appendChild(event.create_copasi_element(dom, event_number)) 
    event_number += 1

  sbml_references = model.getElementsByTagName("SBMLReference")
  if sbml_references:
    sbml_reference = sbml_references[0]
  else:
    sbml_reference = dom.createElement("SBMLReference")
    model.appendChild(sbml_reference)

  for index in range (0, event_number):
    event_name = "TimeCourseEvent_" + str(index)
    sbml_map = dom.createElement("SBMLMap")
    sbml_map.setAttribute("COPASIkey", event_name)
    sbml_map.setAttribute("SBMLid", event_name)
    sbml_reference.appendChild(sbml_map)
    

  check_constant_false (model, species)
 
  if arguments.pretty:
    document = dom.toprettyxml(indent="  ", encoding="UTF-8")
  else:
    document = dom.toxml("UTF-8")
  print(document)


def run():
  """Perform the banalities of command-line argument processing
     and then actually do the main work """
  description = "Analyse SBML files for invariants"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="an sbml file to check invariants for")
  parser.add_argument('--pretty', action='store_true',
                      help="Pretty print the xml")
   
  arguments = parser.parse_args()

  sbml_files = [ x for x in arguments.filenames 
                   if has_xml_or_sbml_ext(x) ]
  copasi_files = [ x for x in arguments.filenames
                      if has_copasi_ext(x) ]
  timecourse_files = [ x for x in arguments.filenames
                          if not x in sbml_files and not x in copasi_files ]

  events = []
  species = []
  for filename in timecourse_files:
    timecourse = timeseries.get_timecourse_from_file(filename)
    species.extend(timecourse.get_column_names())
    these_events = events_from_timecourse(timecourse)
    events.extend(these_events)

  if not events:
    print ("No events to add, doing nothing:")
    sys.exit(1)

  if not sbml_files and not copasi_files:
    for event in events:
      print (event.format_event())
  else:
    for filename in sbml_files:
      add_events_sbml(filename, events, species, arguments)
    for filename in copasi_files:
      add_events_copasi(filename, events, species, arguments)



if __name__ == "__main__":
  run()
