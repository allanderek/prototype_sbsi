"""A simple module to put some simple utility functions in"""
import os.path

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
