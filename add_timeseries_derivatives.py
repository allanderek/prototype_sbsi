""" A script to the derivatives of each species in a timeseries as a
    separate column. Here we use a very simple method to calculate a
    very approximate derivative. Basically we just take the next data
    point and divide by the time between.
"""
import argparse

import timeseries
import modify_timeseries_shell

def add_derivatives_columns(timecourse):
  """Given a timeseries add columns to the end which represent the
     derivatives of all the original columns
  """
  # All the d_ prefixes here stand for 'derivative'
  names = timecourse.get_column_names()
  d_names = [ "Time" ] + [ "d_" + name for name in names ]

  rows = timecourse.get_rows()
  d_rows = []
  # index here goes up to the length of the rows minus 1, because
  # the final row cannot be calculated hence we will somehow fudge this.
  for index in range(0, len(rows) - 1):
    this_row = rows[index]
    next_row = rows[index + 1]
    this_time = this_row[0]
    next_time = next_row[0]
    time_gap = next_time - this_time
    new_d_row = [ this_time ]
    for column_index in range(1, len(this_row)):
      this_value = this_row[column_index]
      next_value = next_row[column_index]
      derivative = (next_value - this_value) / time_gap
      new_d_row.append(derivative)
    d_rows.append(new_d_row)  

  # Remove the last column of the original timeseries such that
  # we can add the two timeseries together
  timecourse.remove_row(len(rows) - 1)

  d_timecourse = timeseries.Timeseries(d_names, d_rows)
  timecourse.add_timeseries(d_timecourse)

  return timecourse

def run():
  """ The main method. 
  """ 
  description = "Add columns to a timeseries representing derivatives"""
  parent_arg_parser = modify_timeseries_shell.create_arguments_parser()
  parser = argparse.ArgumentParser(description=description,
                                   parents=[parent_arg_parser])
                                    
  modify_timeseries_shell.run(parser, add_derivatives_columns, False)

if __name__ == "__main__":
  run()
