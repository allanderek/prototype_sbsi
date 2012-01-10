"""A script to add a noise function to time series data, generally useful
   for having a fake/dry-run of an optimisation problem. One can start from
   known parameters, generate timecourse data, use this script to add
   noise to that data and then attempt to see if you can optimise the model
   to obtain the original known parameters. If not then the problem is
   likely under-constrained
"""
import argparse

import random
import modify_timeseries_shell

def run():
  """ The main method
  """ 
  description = "Add noise to a timeseries"
  parent_arg_parser = modify_timeseries_shell.create_arguments_parser()
  parser = argparse.ArgumentParser(description=description,
                                   parents=[parent_arg_parser])
                                    
  def add_noise_function(timecourse):
    """Simply adds the noise function to the timecourse and returns it"""
    timecourse.apply_noise_function(dream_noise_function)
    return timecourse
  modify_timeseries_shell.run(parser, add_noise_function, False)


  parser = argparse.ArgumentParser(description=description)

def dream_noise_function(orig_value):
  """Transforms the given value using the noise function used in the
     DREAM estimation of parameters competition of 2011"""
  guass_1 = random.gauss(0, 1)
  guass_2 = random.gauss(0, 1) 

  noise_value = orig_value + (0.1 * guass_1) + (0.2 * guass_2 * orig_value)
  return max(0, noise_value)

if __name__ == "__main__":
  run()
