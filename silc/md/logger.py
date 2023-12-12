import copy
import numpy as np

class HistogramLogger:
    """
    Implements a Callback functor for methods.
    Logs the state of the collective variable to generate histograms.


    Parameters
    ----------
    period:
        Time steps between logging of collective variables.

    offset:
        Time steps at the beginning of a run used for equilibration.
    """

    def __init__(self, basename, period: int, offset: int = 0):
        self.basename = basename
        self.period = period
        self.offset = offset
        self.counter = 0
        self.data = None
        self.data_cum = None


    def __call__(self, snapshot, state, timestep):
        """
        Implements the logging itself. Interface as expected for Callbacks.
        """
        self.counter += 1
        if self.counter > self.offset and self.counter % self.period == 0:
            hist = state.hist
            print(hist, self.data_cum, self.data)
            if self.data is None:
                self.data = copy.copy(hist)
                self.data_cum = copy.copy(hist)
            else:
                self.data = hist - self.data_cum
                self.data_cum = copy.copy(hist)
            self.save_file()


    def save_file(self):
        np.savetxt("%s-%d.txt" % (self.basename, self.counter), self.data)

