"""A class implementing an fft cost function for the optimiser"""

import math
import logging


class FFTcost:
  """A cost function based on the fast fourier transform"""
  def __init__(self, gold_ts):
    self.gold_standard_timeseries = gold_ts
    # This clearly isn't correct, since not all of the states
    # may be periodic
    self.state_names = gold_ts.get_column_names()


  @staticmethod
  def dft_single(values, k):
    """Compute the 'energy' at the given frequency of a signal defined
       by the time course of given values. The energy is not relative
       to the energy of the signal as a whole, since it is not compared
       with energies at other frequencies"""
    inv = -1 
    num_values = len(values)
    energy = 0
    for value_index in xrange(num_values) :
      exponent = inv * 2j * math.pi * k * value_index / num_values
      energy += values[value_index] * math.e**exponent
    energy /= num_values
    return abs(energy)
 
  @staticmethod
  def dft(values, kman):
    """Compute the cost using the dft of the signal given as a time
       course of values. The cost is determined against the energy
       of the signal at the given frequency (which in turn iss related
       to the desired period)"""
    num_values = len(values)

    # Rather than find the energy at all frequencies, just
    # figure out the energy at a few selected points 
    k1_energy = FFTcost.dft_single(values, 1)
    kman_energy = FFTcost.dft_single(values, kman)
    khalf_energy = FFTcost.dft_single(values, num_values / 2)
    kquarter_energy = FFTcost.dft_single(values, num_values / 4)
    energies = [ k1_energy, kman_energy, 
                 khalf_energy, kquarter_energy ]
    # energies = [ FFTcost.dft_single(x, n) for n in range(1, N / 2) ]
    # The cost then is how much of the energy at all the frequencies
    # is attributable to the energy at the frequency we are interested
    # in. This is essentially saying how relatively well represented is
    # the given frequency. Note there may be a corner case here in which
    # N is less than the number of frequencies we are considering but
    # I think that is not interesting for our purposes.
    sum_energies = sum(energies)
    # Must be positive since all energies are positive and the sum of
    # all the energies must be at least as large as the single energy
    # we are considering.
    cost = (sum_energies / kman_energy) - kman_energy
    logging.debug("dft_energy       = " + str(kman_energy))
    logging.debug("total_dft_energy = " + str(sum_energies))
    return cost

  def compare_timeseries(self, candidate_ts):
    """Evaluate the given time series based on this DFT cost"""
    # So the target period is actually 
    # Should automatically check if these are in the output region
    # of the solver since otherwise we will miss up.
    start_sampling = 100.0
    stop_sampling  = 300.0
    target_period  = 20.0
    num_cycles     = (stop_sampling - start_sampling) / target_period
    # 300 / 20 = 15 
    # num_cycles = 10 # 15 
    cost = 0
    for column_name in ["mRNA"]: # self.state_names:
      column_data = candidate_ts.get_column_data(column_name)
      
      dft_cost = self.dft(column_data, num_cycles)
      logging.debug("dft_cost = " + str(dft_cost))
      cost += dft_cost 
    return cost


