from dataclasses import dataclass, field
from typing import List

import jax.numpy as np
from jax import grad, jit
from jax.numpy import linalg

from parmed import unit as u

from pysages.colvars.funnels import center, kabsch, periodic
from pysages.colvars.orientation import rmsd, RMSD
from pysages.colvars.orientation import kabsch as rmsd_kabsch

from silc.md.collective_variables import distance_pbc, alignment
from silc.md.util import generate_simulation


@dataclass
class Alignment_Constraint:
    '''
    Contrain the alignment between two objects. 
    The alignment is defined as dot_product(u1, u2)**2, where u1 and u2 are the main axis of the two objects.
    The alignment is defined as dot_product(u1, u2), if the two objects are rods and are asymmetrical with respect to the 180º rotation 
    The alignment is constraint by a harmonic potential with force constant k.

    object1_atoms and object2_atoms: List[int]
        A list of atom index of the object 1 (and object 2)
    object1_rod and object2_rod: bool
        Whether the object 1 (and object 2) is a rod. True if it is a rod. False if it is a plate
    two_rods_symmertric: bool
        If two objects are rods, whether they are symetrical with respect to the 180º rotation 
    minval: float
        min value of the alignment
    maxval: float
        max value of the alignment
    k: float
        force constant that constrains the alignment
    '''
    minval: float
    maxval: float
    k: float
    object1_rod: bool = True
    object2_rod: bool = False
    two_rods_symmertric: bool = False
    object1_atoms: List[int] = field(default_factory=list)
    object2_atoms: List[int] = field(default_factory=list)

    

def alignment_force(pos, ids, alignment_constraint: Alignment_Constraint):
    indices_obj1 = np.array(alignment_constraint.object1_atoms)
    indices_obj2 = np.array(alignment_constraint.object2_atoms)
    if alignment_constraint.object1_rod and alignment_constraint.object2_rod:
        rod = pos[ids[indices_obj1]]
        plate = pos[ids[indices_obj2]]
        two_rods = True
        two_plates = False
    elif alignment_constraint.object1_rod and not alignment_constraint.object2_rod:
        rod = pos[ids[indices_obj1]]
        plate = pos[ids[indices_obj2]]
        two_rods = False
        two_plates = False
    elif not alignment_constraint.object1_rod and alignment_constraint.object2_rod:
        rod = pos[ids[indices_obj12]]
        plate = pos[ids[indices_obj1]]
        two_rods = False
        two_plates = False
    elif not alignment_constraint.object1_rod and not alignment_constraint.object2_rod:
        rod = pos[ids[aindices_obj1]]
        plate = pos[ids[indices_obj2]]
        two_rods = False
        two_plates = True
    minval = alignment_constraint.minval
    maxval = alignment_constraint.maxval
    k = alignment_constraint.k
    asymmetric = alignment_constraint.two_rods_symmertric

    val = alignment(rod, plate, two_plates, two_rods, asymmetric)
    F = np.where(val > maxval, val - maxval, np.where(val < minval, minval - val, 0.0))
    return 0.5 * k * F * F


@dataclass
class Distance_Constraint:
    '''
    Contrain the distance between two objects by harmonic spring, with force constant k.

    object1_atoms and object2_atoms: List[int]
        A list of atom index of the object 1 (and object 2)
    minval: float
        min value of the distance
    maxval: float
        max value of the distance
    k: float
        force constant that constrains the distance
    '''
    minval: float
    maxval: float
    k: float
    object1_atoms: List[int] = field(default_factory=list)
    object2_atoms: List[int] = field(default_factory=list)


def distance_force(pos, ids, box, distance_constraint: Distance_Constraint):
    indices_obj1 = np.array(distance_constraint.object1_atoms)
    indices_obj2 = np.array(distance_constraint.object2_atoms)
    r1 = pos[ids[indices_obj1]]
    r2 = pos[ids[indices_obj2]]
    minval = distance_constraint.minval
    maxval = distance_constraint.maxval
    k = distance_constraint.k
    val = distance_pbc(r1, r2, box)
    F = np.where(val > maxval, val - maxval, np.where(val < minval, minval - val, 0.0))
    return 0.5 * k * F * F


@dataclass
class Recenter_Constraint:
    '''
    Constrain an object to the center of the simulation box

    object_atoms: List[int]
        A list of atom index of the object
    k: float
        force constant that constrains the object to the center of the box
    '''
    k: float
    object_atoms: List[int] = field(default_factory=list)


def recenter_force(pos, ids, box, recenter_constraint: Recenter_Constraint):
    indices_obj = np.array(recenter_constraint.object_atoms)
    k = recenter_constraint.k
    box_center = np.asarray(box)*0.5
    return k * distance_pbc(pos[ids[indices_obj]], box_center, np.asarray(box))**2


@dataclass
class RMSD_Constraint:
    '''
    Constrain the RMSD of an object.

    object: List[int]
        A list of atom index of the object
    '''
    k: float
    object_atoms: List[int] = field(default_factory=list)
    reference_atoms: List[int] = field(default_factory=list)
    reference_files: List[str] = field(default_factory=list)
    
    def __post_init__(self):
        sim = generate_simulation(self.reference_files, minimize_steps=0, NPT_steps=0, NVT_steps=0)
        state = sim.context.getState(getPositions=True)
        ref_pos = state.getPositions(asNumpy=True)
        self.references = ref_pos.value_in_unit(u.nanometer)[self.reference_atoms]


def rmsd_force(pos, ids, rmsd_constraint: RMSD_Constraint):
    indices_obj = np.array(rmsd_constraint.object_atoms)
    references = rmsd_constraint.references
    k = rmsd_constraint.k
    rmsd_restrain = RMSD(rmsd_constraint.object_atoms, references)
    references = rmsd_restrain.Q
    rmsd_w = np.ones(len(indices_obj)) / len(indices_obj)
    return k * rmsd(pos[ids[indices_obj]], references, rmsd_w, rmsd_kabsch)**2


@dataclass
class Funnel_Constraint:
    '''
    Constraint the object molecule/atoms in a funnel cone region. Notice that the projection on the axis is not constraints 
    since it is considered to be a collective variable and is controlled by the parameter `restraints` in `Funnel_ABF`.

    guest_atoms: List[int]
        A list of atom index of the guest object
        
    host_atoms: List[int]
        A list of atom index of the host object

    anchor_atoms: List[int]
        A list of atom index of the anchor (actually should be a list with only one atom index)

    host_reference_atoms: List[int]
        A list of atom index of the host in the reference snapshot served as the reference

    host_reference_files: List[str]
        Simulation files for the reference snapshot. Should be [prmtop file for topology, rst7/xml file for position]

    each_atom: bool
        If False, the funnel constraint is only applied to the center of geometry of the guest object.
        If True, the funnel constraint is applied to each atoms of the guest object
    '''
    Zcc: float
    Z_0: float
    R: float
    k: float
    each_atom: bool = False
    A: np.ndarray = field(default_factory=lambda: np.zeros((3,)))
    B: np.ndarray = field(default_factory=lambda: np.zeros((3,)))
    guest_atoms: List[int] = field(default_factory=list)
    host_atoms: List[int] = field(default_factory=list)
    anchor_atoms: List[int] = field(default_factory=list)
    host_reference_atoms: List[int] = field(default_factory=list)
    host_reference_files: List[str] = field(default_factory=list)


    def __post_init__(self):
        sim = generate_simulation(self.host_reference_files, minimize_steps=0, NPT_steps=0, NVT_steps=0)
        state = sim.context.getState(getPositions=True)
        ref_pos = state.getPositions(asNumpy=True)
        self.references = np.asarray(ref_pos.value_in_unit(u.nanometer)[self.host_reference_atoms])


def y_function(x, Z_0, Zcc, R):
    m = (R - Z_0) / Zcc
    return m * x + Z_0


def cone(x, eje, Zcc, Z_0, R, k):
    x_coord = np.dot(x, eje)
    proj = x_coord * eje
    x_perp = x - proj
    F = linalg.norm(x_perp) - y_function(x_coord, Z_0, Zcc, R)
    return np.where(F < 0.0, 0.0, 0.5 * k * F * F)


def cylinder(x, eje, R, k):
    x_perp = x - np.dot(x, eje) * eje
    F = linalg.norm(x_perp) - R
    return np.where(F < 0.0, 0.0, 0.5 * k * F * F)


def funnel(x, A, B, Zcc, Z_0, R, k):
    A_r = A
    B_r = B
    norm_eje = linalg.norm(B_r - A_r)
    eje = (B_r - A_r) / norm_eje
    #    Z_pos = Zcc * eje
    x_fit = x - A_r
    proj = np.dot(x_fit, eje)
    return np.where(
        proj < Zcc,
        cone(x_fit, eje, Zcc, Z_0, R, k),
        cylinder(x_fit, eje, R, k),
    )


def funnel_force(pos, ids, box, funnel_constraint: Funnel_Constraint):
    indices_guest = np.array(funnel_constraint.guest_atoms)
    indices_host = np.array(funnel_constraint.host_atoms)
    indices_anchor = np.array(funnel_constraint.anchor_atoms)
    A = funnel_constraint.A
    B = funnel_constraint.B
    Zcc = funnel_constraint.Zcc
    Z_0 = funnel_constraint.Z_0
    R = funnel_constraint.R
    k = funnel_constraint.k
    each_atom = funnel_constraint.each_atom
    references = funnel_constraint.references

    weights_guest = None
    weights_host = None
    pos_guest = pos[ids[indices_guest]]
    pos_host = pos[ids[indices_host]]
    pos_anchor = pos[ids[indices_anchor]]
    guest_distances = periodic(pos_guest - pos_anchor, np.asarray(box))
    new_pos_guest = pos_anchor + guest_distances
    center_guest = center(new_pos_guest, weights_guest)
    center_host = center(pos_host, weights_host)
    center_ref = center(references, weights_host)

    if not each_atom:
        guest_rot = np.dot(center_guest - center_host, kabsch(pos_host, references, weights_host)) + center_ref
        return funnel(guest_rot, np.asarray(A), np.asarray(B), Zcc, Z_0, R, k)
    else:
        force = 0
        for pos_atom in new_pos_guest:
            guest_atom_rot = np.dot(pos_atom - center_host, kabsch(pos_host, references, weights_host)) + center_ref
            force += funnel(guest_atom_rot, np.asarray(A), np.asarray(B), Zcc, Z_0, R, k)
        return force