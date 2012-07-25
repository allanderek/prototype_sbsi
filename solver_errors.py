"""A small helper module for solver errors, essentially this is mostly
   here to avoid external_solvers referencing solve_model or alternatively
   putting the solver errors in external_solvers which seems wrong.
"""
class SolverError(Exception):
  """A simple exception to be raised when we recognise that the model
     cannot be solved successfully. This allows callers, such as the
     optimisation routine, to catch this kind of error rather than
     fail completely. Note that this should be used to indicate an
     error in the actual solving of a particular instance of a model,
     not that the model itself cannot be solved. See SolverModelError
     for that.
  """
  def __init__(self, message):
    Exception.__init__(self)
    self.message = message

  def get_message(self):
    """Return the stored message"""
    return self.message

class SolverModelError(SolverError):
  """A simple subclass of SolverError defined above. Essentially,
     for the purposes of implementing many solves of a particular
     model we wish to distinguish between errors which mean that the
     model cannot be solved with the particular parameters, meaning
     for example if you change the parameters it might then be solvable
     and one in which it will never be solvable, for example in the
     Cvodes solver if we fail to be able to convert from SBML to C then
     we will never be able to solve that particular model.
  """
  pass


