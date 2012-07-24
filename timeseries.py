"""A module implementing a class and utilities for dealing with
   time course data and results"""

import sys
import logging
import matplotlib.pyplot as plt

class Timeseries:
  """This class represents a time series, that is the result of a
     numerical evaluation of the given model. It could also be the
     gold standard to which we are comparing the results in order to
     direct the search"""
 
  # The data structure for time series is questionable, in particular
  # it would be nice to at least just use arrays.

  def __init__(self, columns, rows):
    self.columns = columns
    self.rows = rows

  def get_rows(self):
    """Return as a set of rows with the first value being the
       time and the remainder being those of the columns in
       order of the column names"""
    return self.rows

  def remove_row(self, index):
    """Remove a row from the timeseries"""
    self.rows.pop(index)

  def remove_final_row(self):
    """Remove the final row from a timeseries"""
    self.rows.pop(len(self.rows) - 1)

  def get_times(self):
    """Return a list of all times"""
    return [ row[0] for row in self.rows ]

  def get_final_time(self):
    """Return the final time from the timeseries"""
    last_row = self.rows[len(self.rows) - 1]
    return last_row[0]
 
  def get_column_names(self):
    """Returns the set of columns, however does not include the
       first column because that is the 'Time' column"""
    return self.columns[1:]

  def get_column_data(self, column_name, start=None, end=None):
    """Retrieve the data for just a single column"""
    # We should catch the exception (I think ValueError) in case
    # the given column name is not here.
    column_mapping = self.get_column_as_map(column_name, start, end)
    return [ y for (_, y) in column_mapping ]

  def get_column_as_map(self, column_name, start=None, end=None):
    """Retrieve the data for just a single column as a mapping
       from time to value at that time. By mapping it is a list of
       (key, value) pairs, such that it is still in the order of
       the time column order. Additionally the above 'get_column_data'
       depends on this fact"""
    # We should catch the exception (I think ValueError) in case
    # the given column name is not here.
    column_index = self.columns.index(column_name)
    results = []
    for this_row in self.rows:
      row_time = this_row[0]
      if ( (start == None or row_time >= start) and 
           (end == None or row_time <= end) ):
        entry = (row_time, this_row[column_index])
        results.append(entry)
    return results

  def apply_noise_function(self, noise_function):
    """Apply a given noise function to a timeseries. This could actually
       be any function double->double applied to change the data in
       a timecourse file"""
    for row in self.rows:
      for index in range(1, len(row)):
        row[index] = noise_function(row[index])

  def remove_column(self, column_name):
    """Remove the column with the given name from the timeseries"""
    index = self.columns.index(column_name)
    del self.columns[index]
    for row in self.rows:
      del row[index]

  def is_zero_column(self, column_name):
    """Returns true if all of the given column's data is zero"""
    values = self.get_column_data(column_name)
    return all([value == 0.0 for value in values])

  def remove_zero_columns(self):
    """Sometimes timecourse data is given in a format such that for
       a given species if there is no data for that species then the
       column is filled with zeroes. Hence it is occasionally useful
       to remove any column for which all of the values are zero"""
    # Note that this could be done faster from first principles but 
    # is more maintainable using the methods provided above.

    zero_columns = [x for x in self.columns if self.is_zero_column(x) ]
    for column_name in zero_columns:
      self.remove_column(column_name)
   

  def add_timeseries(self, other_timecourse):
    """Add the columns of data from one another time course into this
       one. Currently we require that the times in both timecourses are
       the same. We should at least allow one to be a subset of the others.
       In particular it should be simple enough to allow the additional
       timecourse to have more times, so long as it has all the times in
       this timecourse"""
    these_times = self.get_times()
    those_times = other_timecourse.get_times()
    if these_times == those_times:
      self.columns += other_timecourse.get_column_names()
      for index in range(len(self.rows)):
        other_row = other_timecourse.rows[index]
        self.rows[index] += other_row[1:]
    else:
      print (these_times)
      print ("-------------------")
      print (those_times)
      raise StandardError

  def number_of_data_points (self):
    """Returns the number of data points within a timecourse. Does not
       include the time data points, so this is essentially the number
       of rows multiplied by the number of columns which are not the
       time columns"""
    return len(self.rows) * (len(self.columns) - 1)

  def write_to_file(self, results_file, separator=None, head_comment=None):
    """Format the time series and write it to the given file"""
    if not separator:
      separator = ", "
    if head_comment:
      results_file.write(head_comment)
    prefix = ""
    for column in self.columns:
      results_file.write(prefix)
      results_file.write(column)
      prefix = separator
  
    results_file.write("\n")

    for row in self.rows:
      prefix = ""
      for value in row:
        results_file.write(prefix)
        results_file.write(str(value))
        prefix = separator
      results_file.write("\n")

  def get_best_matching_time_row(self, gold_time):
    """Return the row which has the closest time to the given target
       time. Used for comparing timeseries, we have data in one row
       of one time series and we wish to find the row of this time
       series to which to compare it."""
    best_row_distance = sys.maxint
    best_row = self.rows[0]
    # The search for the best row could be made a bit faster by
    # only starting from the previous best row
    # We could also break as soon as the value is going higher.
    for row in self.rows:
      row_distance = abs(row[0] - gold_time)
      if row_distance < best_row_distance:
        best_row = row  
        best_row_distance = row_distance
    return best_row

  def get_best_matching_rows(self, gold_standard):
    """Return the rows of this time series which best fit the times
       in the rows of the given gold standard time series. We return
       it as a list of pairs matching each gold row to its best matching
       row from this timeseries"""
    best_rows = [ (self.get_best_matching_time_row(gold_row[0]), gold_row)
                  for gold_row in gold_standard.rows ]
    return best_rows

  def re_timealise(self, new_times):
    """This is a fairly bad name, the idea is for this to used to take
       a timeseries straight from a stochastic simulation and create a
       timeseries which has the desired time points. So for each time point
       in the set of desired time points, we create a row with the
       populations as they were at that point. Note that this means it may
       not choose the closest time for any particular time point. Suppose
       we have the time series from a simulation as:
       0.05, 0.11, 0.23, 0.25, 0.27, 0.29
       and the desired times:
       0.1, 0.2, 0.3
       then it will pick out the rows: 0.05, 0.11 and 0.29.
       Essentially for each time point it picks out the latest time point
       which is not higher than the requested time point.
    """
    latest_row = self.rows[0]
    new_rows = []
    index = 1
    index_limit = len(self.rows)

    # Okay so for each time in the set of new times we have to add
    # exactly one row to the new rows.
    for time in new_times:
      # For each new time, check if the current index into the old
      # list of rows points to a row with a time that is less than the
      # new time.
      while self.rows[index][0] < time:
        # If it is, update the 'latest_row' variable to point to that
        # row since we know that it has a time which is before the one
        # we are currently looking for.
        latest_row = self.rows[index]
        # THEN, we should update the index, however we better check
        # that the index is still within the acceptable ranges first,
        # the reason we do this before updating the index is that this
        # may not be the last 'new time' hence we still need the while
        # loop condition to work even though we would be breaking out of
        # this one.
        if index <= index_limit + 1:
          index += 1
        else:
          break
      # Regardless of whether we broke out of the loop because of the
      # index becoming out of range or because of the loop condition
      # we add the latest row to the list of new rows, but we have to
      # substitute in the new time for the row's previous time.
      new_rows.append([time] + latest_row[1:])

    self.rows = new_rows

  def plot_timecourse(self):
    """ Using matplotlib plot the time course """
    plt.figure()

    column_names = self.get_column_names()
    times = self.get_times()
    for name in column_names:
      timeline = self.get_column_data(name)
      plt.plot(times, timeline, label=name)

    plt.xlabel('Time')
    plt.ylabel('Population')
    plt.legend(loc=0)
    plt.show()

def average_timeseries(timecourses):
  """Create a new timeseries with the results being the average of the
     set of given timeseries.
  """
  # TODO: remove the assumptions of the same order of columns and the
  # same times.
  new_rows = []
  minimum = min([ len(t.rows) for t in timecourses ])

  # for t in timecourses:
  #   print (t.get_times())

  for i in range(minimum):
    this_row = []
    rows = [ t.rows[i] for t in timecourses ]
    # This will actually average the time column but they should all
    # be identical anyway.
    for j in range(len(rows[0])):
      values = [ row[j] for row in rows ]
      average = sum(values) / len(values)
      this_row.append(average)
    new_rows.append(this_row)

  return Timeseries(timecourses[0].columns, new_rows)
  
  
  
 

def parse_csv(csv, separator=None):
  """Parse a comma-separated value file into a timeseries"""
  try:
    # If a separator is specified then we try to find it in
    # each line of the file until we do. 
    if separator:
      header_line = csv.next()
      while not separator in header_line: 
        header_line = csv.next()
    else:
      # If no separator is specified then we have the additional
      # task of finding out what the separator is, hence we just
      # try all the ones we can think of (currently only comma and
      # tab). We must be careful not to add too many possibilities here
      # since in that case it may be found in some comment lines above
      # the actual header line.
      separator = ""
      while not separator:
        header_line = csv.next()
        if "," in header_line: 
          separator = ","
        elif "\t" in header_line:
          separator = "\t"
  except StopIteration:
    logging.error("We could not find a separator in the csv file\n")
    if not separator:
      logging.error("We tried: comma and tab")
    sys.exit(1)

  headers = [ x.lstrip().rstrip()
              for x in header_line.split(separator) ]

  rows = [] 

  try:
    while 1: 
      line = csv.next()
      if line:
        values = [ float(x.rstrip().lstrip()) 
                     for x in line.split(separator) ]
        rows.append(values)
  except StopIteration:
    pass

  return Timeseries(headers, rows)
 
def get_timecourse_from_file(filename):
  """Retrieve a time course from the specified file"""
  csvfile = open(filename, "r")
  timecourse = parse_csv(csvfile) 
  csvfile.close()
  return timecourse


