import os
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
    def __init__(self, basename, output_path, period: int, offset: int = 0):
        self.basename = basename
        self.output_path = output_path
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
        filename = os.path.join(self.output_path, "%s-%d.txt" % (self.basename, self.counter))
        np.savetxt(filename, self.data)


class ABFLogger:
    """
    Implements a Callback functor for methods.
    Logs the ABF related variable.


    Parameters
    ----------
    period_hist_force:
        Time steps between logging of histogram and force.

    period_cv:
        Time steps between logging of collective variables.

    offset:
        Time steps at the beginning of a run used for equilibration.
    """
    def __init__(self, basename, output_path, period_hist_force: int, period_cv: int, offset: int = 0):
        self.basename = basename
        self.output_path = output_path
        self.period_hist_force = period_hist_force
        self.period_cv = period_cv
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
                self.hist_cum = copy.copy(hist)
            self.save_file_hist()

            shape = (*state.Fsum.shape[:-1], 1)
            self.force = state.Fsum / np.maximum(state.hist.reshape(shape), 1)
            self.save_file_force()

        if self.counter > self.offset and self.counter % self.period_cv == 0:
            self.xi = copy.copy(state.xi)
            if self.first:
                self.save_file_cv("w")
                self.first = False
            else:
                self.save_file_cv("a")

    def save_file_hist(self):
        if len(self.hist.shape) <= 2:
            filename = os.path.join(self.output_path, "%s-hist-%d.txt" % (self.basename, self.counter))
            np.savetxt(filename, self.hist, fmt="%d")

    def save_file_force(self):
        if self.force.shape[-1] == 1:
            filename = os.path.join(self.output_path, "%s-force-%d.txt" % (self.basename, self.counter))
            print(filename)
            np.savetxt(filename, self.force)
        elif self.force.shape[-1] == 2:
            filename = os.path.join(self.output_path, "%s-force-cv1-%d.txt" % (self.basename, self.counter))
            np.savetxt(filename, self.force[:,:,0])
            filename = os.path.join(self.output_path, "%s-force-cv2-%d.txt" % (self.basename, self.counter))
            np.savetxt(filename, self.force[:,:,1])

    def save_file_cv(self, mode):
        filename = os.path.join(self.output_path, "%s-cv.txt" % self.basename)
        with open(filename, mode) as f:
            f.write("%d" % self.counter)
            for i in range(len(self.xi[0])):
                f.write("\t%f" % self.xi[0][i])
            f.write("\n")
