""" A simple module implementing timeouts which can be used to put a
    limit on the amount of time spent performing some operation. For
    example we may use this when implementing a web service to ensure
    that some user operation does not take a long time and thus use up
    a lot of the system's resources.
"""

import signal

class TimeoutError(Exception):
  """A private class of exception for raising when a timeout occurs"""
  pass


def timeout(timeout_time, return_default=False, default=None):
  """The main timeout function, this is actually a higher order method
     which takes in a method and returns another method and can therefore
     be used as a decorator (although not necessarily as for example
     the timeout may be dynamically calculated).
  """
  def timeout_function(user_method):
    """The actually higher order function which converts a time ignorant
       method into one which will timeout.
    """
    def time_limited_method(*args):
      """The time limited method we return"""
      def timeout_handler(_signum, _frame):
        raise TimeoutError()

      old_handler = signal.signal(signal.SIGALRM, timeout_handler) 
      signal.alarm(timeout_time) # triger alarm in timeout_time seconds
      if return_default:
        try: 
          retval = user_method(*args)
        except TimeoutError:
          return default
        finally:
          signal.signal(signal.SIGALRM, old_handler) 
      else:
        try: 
          retval = user_method(*args)
        finally:
          signal.signal(signal.SIGALRM, old_handler) 

      signal.alarm(0)
      return retval
    return time_limited_method
  return timeout_function
