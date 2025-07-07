# funnel functions
from functools import partial

import jax.numpy as np
from jax import grad, jit
from jax.numpy import linalg

from pysages.colvars.funnels import center, kabsch, periodic
from silc.md.constraints import Alignment_Constraint, Distance_Constraint, Recenter_Constraint, RMSD_Constraint, Funnel_Constraint
from silc.md.constraints import alignment_force, distance_force, recenter_force, rmsd_force, funnel_force


def intermediate_funnel(
    pos,
    ids,
    box,
    alignment_constraints,
    distance_constraints,
    recenter_constraints,
    rmsd_constraints,
    funnel_constraints,
):
    myfunnel = 0
    for alignment_constraint in alignment_constraints:
        myfunnel += alignment_force(pos, ids, alignment_constraint)
    for distance_constraint in distance_constraints:
        myfunnel += distance_force(pos, ids, box, distance_constraint)
    for recenter_constraint in recenter_constraints:
        myfunnel += recenter_force(pos, ids, box, recenter_constraint)
    for rmsd_constraint in rmsd_constraints:
        myfunnel += rmsd_force(pos, ids, rmsd_constraint)
    for funnel_constraint in funnel_constraints:
        myfunnel += funnel_force(pos, ids, box, funnel_constraint)
    return myfunnel


def log_funnel(
    pos,
    ids,
    box,
):
    return None


def external_funnel(
    data,
    box,
    alignment_constraints,
    distance_constraints,
    recenter_constraints,
    rmsd_constraints,
    funnel_constraints,
):
    pos = data.positions[:, :3]
    ids = data.indices
    bias = grad(intermediate_funnel)(
        pos,
        ids,
        box,
        alignment_constraints,
        distance_constraints,
        recenter_constraints,
        rmsd_constraints,
        funnel_constraints,
    )
    proj = log_funnel(
        pos,
        ids,
        box,
    )
    return bias, proj


def get_funnel_force(
    box,
    alignment_constraints,
    distance_constraints,
    recenter_constraints,
    rmsd_constraints,
    funnel_constraints,
):

    funnel_force = partial(
        external_funnel,
        box=box,
        alignment_constraints=alignment_constraints,
        distance_constraints=distance_constraints,
        recenter_constraints=recenter_constraints,
        rmsd_constraints=rmsd_constraints,
        funnel_constraints=funnel_constraints,
    )
    return jit(funnel_force)
