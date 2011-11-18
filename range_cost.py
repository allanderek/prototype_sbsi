"""A module implementing a range cost for the optimiser,
   the cost is zero if all values are within a given range and
   the difference to the limits of the range are squared otherwise.
   This is intended to be used in conjuction with other cost functions
   such as the fft.
"""
# This should be generalised into a range cost.
# It should take the following parameters,
# 1. The columns to check, perhaps all is permitted, or as in
#    the Cholesterol model, all but a given few.
# 2. The time range, rather than look at the gold standard times
#    we should specify a range of times within which we check
#    (perhaps also we can specify an interval, to avoid checking all
#     the output).
# 3. The range values, not that this will not be per column, but a
#    global range, the intention is that this is generally used in
#    conjunction with another cost function.
class RangeCost(object):
  """A special cost function only for the Circadian Clock model.
     However in time we do wish to allow for user cost functions
     and I believe writing them in Python is not a bad way to go"""
  def __init__(self, gold_standard):
    # self.factor = len(gold_standard.rows)
    self.gold_standard = gold_standard
    self.upper_limit = 500
    self.lower_limit = 100
    self.ignored_columns = []

  def set_lower_limit(self, new_lower):
    """Set the lower limit of the range"""
    self.lower_limit = new_lower 

  def set_upper_limit(self, new_upper):
    """Set the upper limit of the range"""
    self.upper_limit = new_upper

  def add_ignored_column(self, column):
    """Add the given column to the list of ignored columns"""
    self.ignored_columns.append(column)

  def set_ignored_columns(self, columns):
    """Set the list of ignored columns to the given list"""
    self.ignored_columns = columns

  def compare_timeseries(self, candidate_ts):
    """cost the time series by checking if all indexes are below"""
    # all_rows = candidate_ts.rows
    # last_row = all_rows[len(all_rows) - 1] 
    cost = 0
    matching_rows = candidate_ts.get_best_matching_rows(self.gold_standard)
    for (best_row, gold_row) in matching_rows:
      # Just skip if the gold standard has a zero time point
      if gold_row[0] <= 0.0:
        continue
      # Now that we've found the best row, we basically ignore the
      # gold standard and just check that all values within that
      # row are within upper and lower limits
      for i in range(1, len(candidate_ts.columns)):
        # For each column we must first check that the species in
        # question is not to be ignored
        candidate_column = candidate_ts.columns[i]
        if candidate_column.lstrip().rstrip() not in self.ignored_columns:
          # if not then great we simply check that the value is within
          # the upper and lower limits and if not then add suitably to
          # the cost
          candidate_value = best_row[i]
          if candidate_value < self.lower_limit:
            diff = self.lower_limit - candidate_value
            cost += diff * diff 
          if candidate_value > self.upper_limit:
            diff = candidate_value - self.upper_limit
            cost += diff * diff
    return cost 


