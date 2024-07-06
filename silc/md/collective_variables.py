from pysages.colvars.core import TwoPointCV, ThreePointCV, FourPointCV
from pysages.colvars.coordinates import barycenter, distance
from jax import numpy as np
from jax.numpy import linalg
from jax import vmap

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

class DistancesProduct(FourPointCV):
    @property
    def function(self):
        return distances_product

def distances_product(r1, r2, r3, r4):
    r1 = barycenter(r1)
    r2 = barycenter(r2)
    r3 = barycenter(r3)
    r4 = barycenter(r4)
    return distance(r1,r2) * distance(r3,r4)


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
        This is a list of two sets of indices. Each sets corresponds to a (part of a) molecule.

    box:
        This is a list of three floats, corresponding to simulation box dimension (assuming cuboid box),
        i.e. [Lx, Ly, Lz]

    Example
    -------

    indices1 = [0, 1, 2, 3, 4, 5]    # atom indices of molecule 1
    indices2 = [10, 11, 12, 13, 14, 15]    # atom indices of molecule 2
    box = [10., 10., 10.]    # simulation box dimension
    cv = DistancePBC([indices1, indices2], box)
    """

    def __init__(self, indices, box=[1., 1., 1.,]):
        super().__init__(indices)
        self.box = np.asarray(box)

    @property
    def function(self):
        return lambda r1, r2: distance_pbc(r1, r2, self.box)


def wrap(distance, box):
    return np.mod(distance + box * 0.5, box) - 0.5 * box


def distance_pbc(r1, r2, box):
    r1 = barycenter(r1)
    r2 = barycenter(r2)
    dr = r1 - r2
    return linalg.norm(wrap(dr, box))


class Alignment(TwoPointCV):
    """
    Collective Variable that calculates the alignment between two groups of particles.
    The alignment is defined as the cos^2(theta), where theta is the angle between the axis
    of two groups of particles. The first group should be rod like, and the second is plate like.

    Parameters
    ----------
    indices: list[int], list[tuple(int)]
        Must be a list or tuple of atoms (ints or ranges) or groups of atoms.
        A group is specified as a nested list or tuple of atoms.
    group_length: int, optional
        Specify if a fixed group length is expected.
    """

    @property
    def function(self):
        """
        Returns
        -------
        Callable
            See `pysages.colvars.pairwise.coordination` for details.
        """
        return align_rod_plate

def mono_inertia(p):
    inertia=np.dot(p,p)*np.identity(3)-np.outer(p,p)
    return inertia

def moment_inertia(positions):
    pos_b = barycenter(positions)
    fit_pos=np.add(positions,-pos_b)
    I=vmap(mono_inertia, in_axes=0)(fit_pos).sum(axis=0)
    return I

def alignment(rod, plate, 2rods=False, assym=False):
    S1 = moment_inertia(rod)
    S2 = moment_inertia(plate)
    _, v1 = linalg.eigh(S1)
    _, v2 = linalg.eigh(S2)
    u1=v1[:,0]    # eigenvector corresponds to the smallest principal moment of inertia, i.e. the axis of rod
    u2=v2[:,-1]   # eigenvector corresponds to the largest principal moment of inertia, i.e. the axis of plate
    
    if 2rods:
        u2=v2[:,0] # eigenvector corresponds to the smallest principal moment of inerta, i.e. the axis of the second rod
    
    dotprod = np.dot(u1,u2)**2 # rods are symmetrical with respect to the 180º rotation 
    
    if assym:
        dotprod = np.dot(u1,u2) # rods are assymetrical with respect to the 180º rotation

    return dotprod

'''
class AlignTwoRods(TwoPointCV):
    """
    Collective Variable that calculates the alignment between two groups of particles.
    The alignment is defined as the cos^2(theta), where theta is the angle between the axis
    of two groups of particles. Both groups should be rod like.

    Parameters
    ----------
    indices: list[int], list[tuple(int)]
        Must be a list or tuple of atoms (ints or ranges) or groups of atoms.
        A group is specified as a nested list or tuple of atoms.
    group_length: int, optional
        Specify if a fixed group length is expected.
    """

    @property
    def function(self):
        """
        Returns
        -------
        Callable
            See `pysages.colvars.pairwise.coordination` for details.
        """
        return align_two_rods

def align_two_rods(rod1, rod2):
    S1 = moment_inertia(rod1)
    S2 = moment_inertia(rod2)
    _, v1 = linalg.eigh(S1)
    _, v2 = linalg.eigh(S2)
    u1=v1[:,0]    # eigenvector corresponds to the smallest principal moment of inertia, i.e. the axis of rod 1
    u2=v2[:,0]   # eigenvector corresponds to the smallest principal moment of inertia, i.e. the axis of rod 2
    return np.dot(u1,u2)
'''
