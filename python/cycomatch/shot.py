"""SHOT descriptor (Tombari et al. 2010), replacing the calc_shot mex file.

Follows the PCL SHOT352 implementation: local reference frame from a
distance-weighted covariance, 32 spatial sectors (8 azimuth x 2 elevation x
2 radial) with (nbins + 1) cosine bins each, and quadrilinear interpolation.
"""
import numpy as np
from scipy.spatial import cKDTree


def _local_frames(X, nbr_i, nbr_j, radius):
    n = X.shape[0]
    d = X[nbr_j] - X[nbr_i]
    w = radius - np.linalg.norm(d, axis=1)
    cov = np.zeros((n, 3, 3))
    np.add.at(cov, nbr_i, w[:, None, None] * d[:, :, None] * d[:, None, :])
    wsum = np.bincount(nbr_i, weights=w, minlength=n)
    cov /= np.maximum(wsum, 1e-12)[:, None, None]
    _, evecs = np.linalg.eigh(cov)  # ascending eigenvalues
    x_ax = evecs[:, :, 2]
    z_ax = evecs[:, :, 0]

    # sign disambiguation: point the axis towards the majority of neighbours
    for ax in (x_ax, z_ax):
        s = np.sign(np.einsum("ij,ij->i", d, ax[nbr_i]))
        vote = np.bincount(nbr_i, weights=np.where(s >= 0, 1.0, -1.0), minlength=n)
        ax[vote < 0] *= -1
    y_ax = np.cross(z_ax, x_ax)
    return x_ax, y_ax, z_ax


def shot_descriptor(X, normals, nbins=10, radius=0.05):
    n = X.shape[0]
    tree = cKDTree(X)
    nbrs = tree.query_ball_point(X, radius)
    nbr_i = np.repeat(np.arange(n), [len(l) for l in nbrs])
    nbr_j = np.concatenate([np.asarray(l, dtype=np.int64) for l in nbrs])

    x_ax, y_ax, z_ax = _local_frames(X, nbr_i, nbr_j, radius)

    delta = X[nbr_j] - X[nbr_i]
    dist = np.linalg.norm(delta, axis=1)
    keep = dist > 0
    nbr_i, nbr_j, delta, dist = nbr_i[keep], nbr_j[keep], delta[keep], dist[keep]

    x = np.einsum("ij,ij->i", delta, x_ax[nbr_i])
    y = np.einsum("ij,ij->i", delta, y_ax[nbr_i])
    z = np.einsum("ij,ij->i", delta, z_ax[nbr_i])
    x[np.abs(x) < 1e-30] = 0
    y[np.abs(y) < 1e-30] = 0
    z[np.abs(z) < 1e-30] = 0

    cos_desc = np.clip(np.einsum("ij,ij->i", normals[nbr_j], z_ax[nbr_i]), -1, 1)
    bin_dist = (1.0 + cos_desc) * nbins / 2.0

    # sector index (same bit trick as PCL)
    bit4 = (y > 0) | ((y == 0) & (x < 0))
    bit3 = np.where((x > 0) | ((x == 0) & (y > 0)), ~bit4, bit4)
    desc = ((bit4.astype(int) << 3) + (bit3.astype(int) << 2)) << 1
    same = (x * y > 0) | (x == 0)
    desc += np.where(same, np.where(np.abs(x) >= np.abs(y), 0, 4),
                     np.where(np.abs(x) > np.abs(y), 4, 0))
    desc += (z > 0).astype(int)
    r12, r14, r34 = radius / 2, radius / 4, radius * 3 / 4
    desc += np.where(dist > r12, 2, 0)

    step = np.where(bin_dist < 0, np.ceil(bin_dist - 0.5), np.floor(bin_dist + 0.5)).astype(int)
    bin_dist = bin_dist - step
    nb1 = nbins + 1
    shot = np.zeros(n * 32 * nb1)
    base = nbr_i * 32 * nb1

    def add(sector, b, val):
        np.add.at(shot, base + sector * nb1 + b, val)

    w_int = 1 - np.abs(bin_dist)
    # cosine interpolation
    pos = bin_dist > 0
    add(desc, np.where(pos, (step + 1) % nbins, (step - 1 + nbins) % nbins),
        np.where(pos, bin_dist, -bin_dist))

    # radial interpolation
    outer = dist > r12
    rd_out = (dist - r34) / r12
    rd_in = (dist - r14) / r12
    w_int += np.where(outer,
                      np.where(dist > r34, 1 - rd_out, 1 + rd_out),
                      np.where(dist < r14, 1 + rd_in, 1 - rd_in))
    m = outer & (dist <= r34)
    add(np.where(m, desc - 2, desc), step, np.where(m, -rd_out, 0))
    m = ~outer & (dist >= r14)
    add(np.where(m, desc + 2, desc), step, np.where(m, rd_in, 0))

    # elevation interpolation
    incl = np.arccos(np.clip(z / dist, -1, 1))
    lower = (incl > np.pi / 2) | ((np.abs(incl - np.pi / 2) < 1e-30) & (z <= 0))
    id_low = (incl - 3 * np.pi / 4) / (np.pi / 2)
    id_up = (incl - np.pi / 4) / (np.pi / 2)
    w_int += np.where(lower,
                      np.where(incl > 3 * np.pi / 4, 1 - id_low, 1 + id_low),
                      np.where(incl < np.pi / 4, 1 + id_up, 1 - id_up))
    m = lower & (incl <= 3 * np.pi / 4)
    add(np.where(m, desc + 1, desc), step, np.where(m, -id_low, 0))
    m = ~lower & (incl >= np.pi / 4)
    add(np.where(m, desc - 1, desc), step, np.where(m, id_up, 0))

    # azimuth interpolation
    has_az = (x != 0) | (y != 0)
    az = np.arctan2(y, x)
    sel = desc >> 2
    az_d = (az - (-np.pi / 8 + (np.pi / 4) * sel)) / (np.pi / 4)
    az_d = np.clip(az_d, -0.5, 0.5)
    az_d = np.where(has_az, az_d, 0)
    w_int += np.where(has_az, np.where(az_d > 0, 1 - az_d, 1 + az_d), 0)
    interp = np.where(az_d > 0, (desc + 4) % 32, (desc - 4 + 32) % 32)
    add(interp, step, np.abs(az_d))

    add(desc, step, w_int)

    shot = shot.reshape(n, 32 * nb1)
    shot /= np.maximum(np.linalg.norm(shot, axis=1, keepdims=True), 1e-12)
    return shot
