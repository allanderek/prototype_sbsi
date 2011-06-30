"""A unit test module for the parameters module"""
import parameters
import unittest

class Arguments:
  """Emulates the arguments class which can be passed into
     functions of the parameters module such as 'check_parameters'
  """
  def __init__(self):
    self.tolerance = 0.01

class TestCheckParameters(unittest.TestCase):
  """ A class designed to test the check_parameter function of
      the parameters module
  """
  # Okay we can't do much about the unittest.TestCase having
  # some 45 public methods.k
  # pylint: disable-msg=R0904
  def setUp(self):
    # Similarly we can't do much about the name 'setUp' 
    # pylint: disable-msg=C0103
    self.param_p1 = parameters.Parameter("p1", 1.0,
                                         0.1, 10.0)
    self.best_params = dict()
    self.best_params["p1"] = 3.0
    self.arguments = Arguments()

  def test_check_param(self):
    """Checks the parameters.check_parameter function"""
    p1_result = parameters.check_parameter(self.param_p1,
                                           self.best_params["p1"], 
                                           self.arguments) 
    self.assertIsNone(p1_result)

if __name__ == '__main__':
  unittest.main()
