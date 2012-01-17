"""A unit test module for the function_translator module"""
import function_translator
import unittest

class Arguments:
  """Emulates the arguments class which can be passed into
     functions of the parameters module such as 'check_parameters'
  """
  def __init__(self):
    self.tolerance = 0.01


class TestCheckFunctionTranslator(unittest.TestCase):
  """ A class designed to test the  function of function_translator module.
  """
  # Okay we can't do much about the unittest.TestCase having
  # some 45 public methods.k
  # pylint: disable=R0904
  def setUp(self):
    # Similarly we can't do much about the name 'setUp' 
    # pylint: disable=C0103
    pass 

  def test_check_function_parser(self):
    """Simply check the parsing of a function definition file"""
    my_source = """f a b = a + b ; g x y = x + y ;
                """
    parse_results = function_translator.parse_function_list(my_source)

    expected_number_fund_defs = 2
    fun_defs = [ x for x in parse_results
                   if isinstance(x, function_translator.FunctionDefinition)
               ]
    actual_number_fun_defs = len(fun_defs)
    self.assertEqual(expected_number_fund_defs,
                     actual_number_fun_defs)

if __name__ == '__main__':
  unittest.main()
