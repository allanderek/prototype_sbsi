"""A unit test module for the solve_model module"""
import solve_model
import unittest

class Arguments:
  """Emulates the arguments class which can be passed into
     functions of the parameters module such as 'check_parameters'
  """
  def __init__(self):
    self.tolerance = 0.01


class TestCheckSolveModel(unittest.TestCase):
  """ A class designed to test the  function of
      the parameters module
  """
  # Okay we can't do much about the unittest.TestCase having
  # some 45 public methods.k
  # pylint: disable-msg=R0904
  def setUp(self):
    # Similarly we can't do much about the name 'setUp' 
    # pylint: disable-msg=C0103
    self.arg_parser    = solve_model.create_arguments_parser(True)

  def test_check_biopepa_cvodes(self):
    """Simply check the solving of a biopepa model, using the cvodes
       solver, by first converting the biopepa to xml"""
    filename = "test_data/mm/mm.biopepa"
    arguments = self.arg_parser.parse_args(["--solver", "cvodes"])

    arguments.stop_time = 1.0

    solver = solve_model.get_solver(filename, arguments)
    solver.initialise_solver()
    timecourse = solver.solve_model(arguments)
   
    # And now we actually need something to assert
    expected_final_time = 1.0
    final_time = timecourse.get_final_time()
    self.assertEqual(expected_final_time, final_time)

    expected_column_names = ["E", "ES", "P", "S" ]
    column_names = timecourse.get_column_names()
    self.assertEqual(expected_column_names, column_names)

if __name__ == '__main__':
  unittest.main()
