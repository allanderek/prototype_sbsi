"""A simple module to put some simple utility functions in"""
import os.path

def change_filename_ext(filename, new_ext):
  """Returns a new file name based on the first but with the extension
     replaced by the given extension. The new_ext should include the
     '.' separator if that is desired"""
  basename = os.path.splitext(filename)[0]
  return basename + new_ext
 
