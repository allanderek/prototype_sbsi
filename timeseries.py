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
    return [ y for (x,y) in column_mapping ]

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



  def write_to_file(self, results_file):
    """Format the time series and write it to the given file"""
    results_file.write("# ")
    prefix = ""
    for column in self.columns:
      results_file.write(prefix)
      results_file.write(column)
      prefix = ", "
  
    results_file.write("\n")

    for row in self.rows:
      prefix = ""
      for value in row:
        results_file.write(prefix)
        results_file.write(str(value))
        prefix = ", "
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
        if "," in header_line: separator = ","
        elif "\t" in header_line: separator = "\t"
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
        values = [ float(x.rstrip()) for x in line.split(separator) ]
        rows.append(values)
  except StopIteration:
    pass

  return Timeseries(headers, rows)
 
