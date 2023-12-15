from pysages.colvars.core import ThreePointCV
from pysages.colvars.coordinates import barycenter, distance


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