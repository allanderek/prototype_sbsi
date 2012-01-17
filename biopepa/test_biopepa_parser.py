"""A unit test module for the biopepa_parser module"""
import unittest

import biopepa.biopepa_parser as biopepa_parser

class Arguments:
  """Emulates the arguments class which can be passed into
     functions of the parameters module such as 'check_parameters'
  """
  def __init__(self):
    self.tolerance = 0.01


class TestCheckBiopepaParser(unittest.TestCase):
  """ A class designed to test the  function of biopepa_parser module.
  """
  # Okay we can't do much about the unittest.TestCase having
  # some 45 public methods.k
  # pylint: disable=R0904
  def setUp(self):
    # Similarly we can't do much about the name 'setUp' 
    # pylint: disable=C0103
    pass 

  def test_check_biopepa_parser(self):
    """Simply check the parsing of a Bio-PEPA model"""
    my_source = """
        x = 1.0;
        y = 2.0 + 3.0;
        zz = 1 + 2 * 3;
        xy = x + y;
        ar = f(1) + g(2);
        f = f();
        g = g(1.0);
        h = h(1.0, 2.0);
        i = k(h(1.0));
        A = x >> ;
        B = y >> + x << ;
        A[0] <*> B[x + y]
        """
    biopepa_model = biopepa_parser.parse_model(my_source)
    
    expected_number_var_defs = 9
    actual_number_var_defs = len(biopepa_model.var_defs)
    self.assertEqual(expected_number_var_defs,
                     actual_number_var_defs)


    expected_number_of_comp_defs = 2
    actual_number_of_component_defs = len(biopepa_model.component_defs)
    self.assertEqual(expected_number_of_comp_defs, 
                     actual_number_of_component_defs)

if __name__ == '__main__':
  unittest.main()
