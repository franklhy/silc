# funnel functions
from functools import partial

import jax.numpy as np
from jax import grad, jit
from jax.numpy import linalg

from pysages.colvars.funnels import center, kabsch, periodic

from silc.md.collective_variables import alignment, wrap, distance_pbc  

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

def alignforce(rod1, rod2, minval, maxval, k):
    val = alignment(rod1, rod2, tworods)
    F = np.where(val > maxval, val - maxval, np.where(val < minval, minval - val, 0.0))
    return 0.5 * k * F * F

def distanceforce(r1, r2, box, minval, maxval, k):
    val = distance_pbc(r1, r2, box)
    F = np.where(val > maxval, val - maxval, np.where(val < minval, minval - val, 0.0))
    return 0.5 * k * F * F

def borderU(x, eje, k, cv_max):
    # x = lig_com-A
    proj_restr = np.dot(x, eje)
    B = proj_restr - cv_max  # upper limit
    return np.where(B < 0.0, 0.0, 0.5 * k * B * B)


def borderL(x, eje, k, cv_min):
    # x = lig_com-A
    proj = np.dot(x, eje)
    B = proj - cv_min  # lower limit
    return np.where(B > 0.0, 0.0, 0.5 * k * B * B)


def rotation_lig(pos_lig, pos_ref, references, weights, com_ref):
    com_prot = center(pos_ref, weights)
    lig_rot = np.dot(pos_lig - com_prot, kabsch(pos_ref, references, weights))
    return lig_rot + com_ref


def funnel(x, A, B, Zcc, Z_0, R, k, k_cv, cv_min, cv_max):
    A_r = A
    B_r = B
    norm_eje = linalg.norm(B_r - A_r)
    eje = (B_r - A_r) / norm_eje
    #    Z_pos = Zcc * eje
    x_fit = x - A_r
    proj = np.dot(x_fit, eje)
    return np.where(
        proj < Zcc,
        cone(x_fit, eje, Zcc, Z_0, R, k) + borderL(x_fit, eje, k_cv, cv_min),
        cylinder(x_fit, eje, R, k) + borderU(x_fit, eje, k_cv, cv_max),
    )


def proj_funnel(x, A, B, Zcc, Z_0, R, k, k_cv, cv_min, cv_max):
    A_r = A
    B_r = B
    norm_eje = linalg.norm(B_r - A_r)
    eje = (B_r - A_r) / norm_eje
    #    Z_pos = Zcc * eje
    x_fit = x - A_r
    proj = np.dot(x_fit, eje)
    perp = x_fit - proj * eje
    return linalg.norm(perp)


def intermediate_funnel(
    pos,
    ids,
    indexes,
    references,
    weights_ligand,
    weights_protein,
    A,
    B,
    Zcc,
    Z_0,
    R,
    k,
    k_cv,
    cv_min,
    cv_max,
    box,
    minalign,
    maxalign,
    k_algn,
    mindis,
    maxdis,
    k_dis
):

    indices_ligand = np.array(indexes[0])
    indices_protein = np.array(indexes[1])
    indices_anchor = np.array(indexes[2])
    indices_tail = np.array(indexes[3])
    indices_ligand2 = np.array(indexes[4])

    
    pos_anchor = pos[ids[indices_anchor]]
    pos_protein = pos[ids[indices_protein]]
    ligand_distances = periodic(pos[ids[indices_ligand]] - pos_anchor, np.asarray(box))
    tail_distances = periodic(pos[ids[indices_tail]] - pos_anchor, np.asarray(box))
    new_lig_pos = pos_anchor + ligand_distances
    pos_ligand = center(new_lig_pos, weights_ligand)
    pos_tail = pos_anchor + tail_distances
    pos_ref = center(np.asarray(references), weights_protein)
    lig_rot = rotation_lig(
        pos_ligand, pos_protein, np.asarray(references), weights_protein, pos_ref
    )
    tail_rot = []
    pos_tail_array = []
    for i in range(len(pos_tail)):
        pos_tail_array.append(
            center(pos_tail[i], None)
        )
        tail_rot.append(
            rotation_lig(
                pos_tail_array[i], pos_protein, np.asarray(references), weights_protein, pos_ref
            )
        )
    myfunnel = funnel(lig_rot, np.asarray(A), np.asarray(B), Zcc, Z_0, R, k, k_cv, cv_min, cv_max)
    
    myfunnel += alignforce(indices_ligand, indices_ligand2, minalign, maxalign, k_algn)
    
    myfunnel += distanceforce(indices_ligand, indices_ligand2, box, mindis, maxdis, k_dis)

    for i in range(len(pos_tail_array)):
        myfunnel += funnel(tail_rot[i], np.asarray(A), np.asarray(B), Zcc, Z_0, R, k, k_cv, cv_min, cv_max)
    return myfunnel


def log_funnel(
    pos,
    ids,
    indexes,
    references,
    weights_ligand,
    weights_protein,
    A,
    B,
    Zcc,
    Z_0,
    R,
    k,
    k_cv,
    cv_min,
    cv_max,
    box,
):
    indices_ligand = np.array(indexes[0])
    indices_protein = np.array(indexes[1])
    indices_anchor = np.array(indexes[2])
    indices_tail = np.array(indexes[3])
    indices_ligand2 = np.array(indexes[4])

    pos_anchor = pos[ids[indices_anchor]]
    pos_protein = pos[ids[indices_protein]]
    ligand_distances = periodic(pos[ids[indices_ligand]] - pos_anchor, np.asarray(box))
    tail_distances = periodic(pos[ids[indices_tail]] - pos_anchor, np.asarray(box))
    new_lig_pos = pos_anchor + ligand_distances
    pos_ligand = center(new_lig_pos, weights_ligand)
    pos_tail = pos_anchor + tail_distances
    pos_ref = center(np.asarray(references), weights_protein)
    lig_rot = rotation_lig(
        pos_ligand, pos_protein, np.asarray(references), weights_protein, pos_ref
    )
    tail_rot = rotation_lig(
        pos_tail, pos_protein, np.asarray(references), weights_protein, pos_ref
    )
    return proj_funnel(lig_rot, np.asarray(A), np.asarray(B), Zcc, Z_0, R, k, k_cv, cv_min, cv_max)


def external_funnel(
    data,
    indexes,
    references,
    weights_ligand,
    weights_protein,
    A,
    B,
    Zcc,
    Z_0,
    R,
    k,
    k_cv,
    cv_min,
    cv_max,
    box,
    minalign,
    maxalign,
    k_algn,
    mindis,
    maxdis,
    k_dis
):
    pos = data.positions[:, :3]
    ids = data.indices
    bias = grad(intermediate_funnel)(
        pos,
        ids,
        indexes,
        references,
        weights_ligand,
        weights_protein,
        A,
        B,
        Zcc,
        Z_0,
        R,
        k,
        k_cv,
        cv_min,
        cv_max,
        box,
        minalign,
        maxalign,
        k_algn,
        mindis,
        maxdis,
        k_dis
    )
    proj = log_funnel(
        pos,
        ids,
        indexes,
        references,
        weights_ligand,
        weights_protein,
        A,
        B,
        Zcc,
        Z_0,
        R,
        k,
        k_cv,
        cv_min,
        cv_max,
        box,
    )
    return bias, proj


def get_funnel_force(
    indices_sys,
    ref_positions,
    A,
    B,
    Zcc,
    Z_0,
    R_cyl,
    k_cone,
    k_cv,
    cv_min,
    cv_max,
    cv_buffer,
    box,
    minalign = 0.0,
    maxalign = 1.0,
    k_algn = 0.,
    mindis = 0.4,
    maxdis = 1.,
    k_dis = 0.,
    w_ligand=None,
    w_protein=None,
):

    funnel_force = partial(
        external_funnel,
        indexes=indices_sys,
        references=ref_positions,
        weights_ligand=w_ligand,
        weights_protein=w_protein,
        A=A,
        B=B,
        Zcc=Zcc,
        Z_0=Z_0,
        R=R_cyl,
        k=k_cone,
        k_cv=k_cv,
        cv_min=cv_min - cv_buffer,
        cv_max=cv_max + cv_buffer,
        box=box,
        minalign = minalign,
        maxalign = maxalign,
        k_algn = k_algn,
        mindis = mindis,
        maxdis = maxdis,
        k_dis = k_dis
    )
    return jit(funnel_force)
