from pysages.colvars.core import TwoPointCV, ThreePointCV
from pysages.colvars.coordinates import barycenter, distance
from jax import numpy as np
from jax.numpy import linalg

class DistancesSum(ThreePointCV):
    """
    Parameters
    ----------

    indices:
        This is a list of three set of indices. The first one corresponds to the bridge,
        and the rest correspond to the cores.


    Example
    -------

    indices1 = [0, 1, 2, 3, 4, 5]    # atom indices of bridge
    indices2 = [10, 11, 12, 13, 14, 15]    # atom indices of one core
    indices3 = [20, 21, 22, 23, 24, 25]    # atom indices of another core
    cv = DistanceSum([indices1, indices2, indices3])
    """

    @property
    def function(self):
        return distance_sum


def distance_sum(r1, r2, r3):
    r1 = barycenter(r1)
    r2 = barycenter(r2)
    r3 = barycenter(r3)
    return distance(r1, r2) + distance(r1, r3)


class DistancesSum2(ThreePointCV):
    """
    Parameters
    ----------

    indices:
        This is a list of three set of indices. The first one corresponds to the bridge,
        and the rest correspond to the cores.


    Example
    -------

    indices1 = [0, 1, 2, 3, 4, 5]    # atom indices of bridge
    indices2 = [10, 11, 12, 13, 14, 15]    # atom indices of one core
    indices3 = [20, 21, 22, 23, 24, 25]    # atom indices of another core
    cv = DistanceSum([indices1, indices2, indices3])
    """

    @property
    def function(self):
        return distance_sum2


def distance_sum2(r1, r2, r3):
    r1 = barycenter(r1)
    r2 = barycenter(r2)
    r3 = barycenter(r3)
    return distance(r1, r2) + distance(r1, r3) + distance(r2, r3)

class DistancesSum3(ThreePointCV):
    """
    Parameters
    ----------

    indices:
        This is a list of three set of indices. The first one corresponds to the bridge,
        and the rest correspond to the cores.


    Example
    -------

    indices1 = [0, 1, 2, 3, 4, 5]    # atom indices of bridge
    indices2 = [10, 11, 12, 13, 14, 15]    # atom indices of one core
    indices3 = [20, 21, 22, 23, 24, 25]    # atom indices of another core
    cv = DistanceSum([indices1, indices2, indices3])
    """

    @property
    def function(self):
        return distance_sum3


def distance_sum3(r1, r2, r3):
    r1 = barycenter(r1)
    r2 = barycenter(r2)
    r3 = barycenter(r3)
    return distance(r1, r2) + distance(r1, r3) - distance(r2, r3)


class CenterOfMassAngle(ThreePointCV):
    """
    Parameters
    ----------

    indices:
        This is a list of three set of indices. The middle one corresponds to the bridge,
        and the rest correspond to the cores.


    Example
    -------

    indices1 = [0, 1, 2, 3, 4, 5]    # atom indices of bridge
    indices2 = [10, 11, 12, 13, 14, 15]    # atom indices of one core
    indices3 = [20, 21, 22, 23, 24, 25]    # atom indices of another core
    cv = CenterOfMassAngle([indices2, indices1, indices3])
    """

    @property
    def function(self):
        return com_angle


def com_angle(r1, r2, r3):
    r1 = barycenter(r1)
    r2 = barycenter(r2)
    r3 = barycenter(r3)
    q = r1 - r2
    r = r3 - r2
    return np.arctan2(linalg.norm(np.cross(q, r)), np.dot(q, r))


class DistancePBC(TwoPointCV):
    """
    Parameters
    ----------

    indcies:
    """

    def __init__(self, indices, box=[1., 1., 1.,]):
        super().__init__(indices)
        self.box = np.asarray(box)

    @property
    def function(self):
        return lambda r1, r2: distance_pbc(r1, r2, self.box)


def periodic(distance, box):
    return np.mod(distance + box * 0.5, box) - 0.5 * box

def distance_pbc(r1, r2, box):
    r1 = barycenter(r1)
    r2 = barycenter(r2)
    dr = r1 - r2
    return linalg.norm(periodic(dr, box))