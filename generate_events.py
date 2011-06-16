""" A simple python script which should, given time course data,
    produce a series of events to be added to an SBML model, such
    that the events simulate the changing of the population as in
    the time course data."""
import os
import sys
import argparse
from subprocess import Popen
import optimiser
import xml.dom.minidom

class EventAssign:
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

  def format_event(self):
    element = self.create_element()
    formatted = element.toprettyxml()
    return formatted

  def create_element(self, document):
    event = document.createElement("event")
    trigger = document.createElement("trigger")
    event.appendChild(trigger)
    math_ns = "http://www.w3.org/1998/Math/MathML"
    trigger_math = document.createElementNS(math_ns, "math")
    trigger.appendChild(trigger_math)
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
    loea = document.createElement("listOfEventAssignments")
    event.appendChild(loea)
    for event_assign in self.event_assigns:
      e_assign = document.createElement("eventAssignment")
      loea.appendChild(e_assign)
      e_assign.setAttribute("variable", event_assign.species)
      e_math = document.createElementNS(math_ns, "math")
      e_assign.appendChild(e_math)
      e_cn = document.createElement("cn")
      e_math.appendChild(e_cn)
      e_cn_text = document.createTextNode(str(event_assign.value)) 
      e_cn.appendChild(e_cn_text)
      
    return event
    

def events_from_timecourse(timecourse):
  events = []
  column_names = timecourse.get_column_names()
  for row in timecourse.get_rows():
    time = row[0]
    event_assigns = []
    for i in range(1, len(column_names)):
      name = column_names[i]
      value = row[i]
      event_assign = EventAssign(name, value)
      event_assigns.append(event_assign)
    event = Event(time, event_assigns)
    events.append(event)
  return events
  

def get_events_for_file(filename):
  csvfile = open(filename, "r")
  timecourse = optimiser.parse_csv(csvfile) 
  csvfile.close()
  events = events_from_timecourse(timecourse)
  return events 
 
def has_xml_or_sbml_ext(filename):
  extension = os.path.splitext(filename)[1]
  return extension == ".xml" or extension == ".sbml"

def add_events(filename, events):
  """Parse in a file as an SBML model, obfuscate it and then print out
     the obfuscated model"""
  dom = xml.dom.minidom.parse(filename)
  model = dom.getElementsByTagName("model")[0]
  
  # Where loe = listOfEvents
  loe_elements = model.getElementsByTagName("listOfEvents")
  if not loe_elements:
    loe_element = dom.createElement("listOfEvents")
    # listOfEvents is the last child of a model so that we can
    # just append it regardless of what other child nodes are in
    # this particular model.
    model.appendChild(loe_element)
  else:
    loe_element = loe_elements[0]

  for event in events:
    loe_element.appendChild(event.create_element(dom)) 
 
  # pretty should be a command-line option, we mostly won't need it. 
  print(dom.toprettyxml())
 

def run():
  """Perform the banalities of command-line argument processing
     and then actually do the main work """
  description = "Analyse SBML files for invariants"
  parser = argparse.ArgumentParser(description=description)
  # Might want to make the type of this 'FileType('r')'
  parser.add_argument('filenames', metavar='F', nargs='+',
                      help="an sbml file to check invariants for")
 
  arguments = parser.parse_args()

  timecourse_files = [ x for x in arguments.filenames 
                         if not has_xml_or_sbml_ext(x) ]

  sbml_files = [ x for x in arguments.filenames 
                   if has_xml_or_sbml_ext(x) ]

 
  events = []
  for filename in timecourse_files:
    events.extend(get_events_for_file(filename))

  if not sbml_files:
    for event in events:
      print (event.format_event())
  else:
    for filename in sbml_files:
      add_events(filename, events)



if __name__ == "__main__":
  run()
