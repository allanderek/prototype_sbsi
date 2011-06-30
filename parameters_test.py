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
    self.arguments = Arguments()

  def test_check_param(self):
    """Checks the parameters.check_parameter function"""
    param_p1  = parameters.Parameter("p1", 1.0, 0.1, 10.0)
    p1_result = parameters.check_parameter(param_p1, 3.0,
                                           self.arguments) 
    self.assertIsNone(p1_result)

    param_p2  = parameters.Parameter("p2", 1.0, 0.0, 10.0)
    p2_result = parameters.check_parameter(param_p2, 0.00001,
                                           self.arguments)
    self.assertIsNotNone(p2_result)
    self.assertTrue(p2_result.too_low)

if __name__ == '__main__':
  unittest.main()
