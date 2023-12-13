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


class ABFLogger:
    """
    Implements a Callback functor for methods.
    Logs the ABF related variable.


    Parameters
    ----------
    period_hist_force:
        Time steps between logging of histogram and force.

    period_CV:
        Time steps between logging of collective variables.

    offset:
        Time steps at the beginning of a run used for equilibration.
    """
    def __init__(self, basename, period_hist_force: int, period_CV: int, offset: int = 0):
        self.basename = basename
        self.period_hist_force = period_hist_force
        self.period_CV = period_CV
        self.offset = offset
        self.counter = 0
        self.hist = None
        self.hist_cum = None
        self.force = None
        self.xi = None
        self.first = True

    def __call__(self, snapshot, state, timestep):
        """
        Implements the logging itself. Interface as expected for Callbacks.
        """
        self.counter += 1

        if self.counter > self.offset and self.counter % self.period_hist_force == 0:
            hist = state.hist
            if self.hist is None:
                self.hist = copy.copy(hist)
                self.hist_cum = copy.copy(hist)
            else:
                self.hist = hist - self.hist_cum
                self.dhist_cum = copy.copy(hist)
            self.save_file_hist()

            shape = (*state.Fsum.shape[:-1], 1)
            self.force = state.Fsum / np.maximum(state.hist.reshape(shape), 1)
            self.save_file_force()

        if self.counter > self.offset and self.counter % self.period_CV == 0:
            self.xi = copy.copy(state.xi)
            if self.first:
                self.save_file_CV("w")
                self.first = False
            else:
                self.save_file_CV("a")

    def save_file_hist(self):
        np.savetxt("%s-hist-%d.txt" % (self.basename, self.counter), self.hist)

    def save_file_force(self):
        np.savetxt("%s-force-%d.txt" % (self.basename, self.counter), self.force)

    def save_file_CV(self, mode):
        with open("%s-cv.txt" % self.basename, mode) as f:
            f.write("%d\t%f\n" % (self.counter, self.xi))