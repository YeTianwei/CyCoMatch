"""Mesh I/O and preprocessing (port of tools/+MESH/preprocess.m and friends)."""
import os

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh


def read_ply(path):
    from plyfile import PlyData

    ply = PlyData.read(path)
    v = ply["vertex"]
    verts = np.stack([v["x"], v["y"], v["z"]], axis=1).astype(np.float64)
    faces = np.vstack(ply["face"]["vertex_indices"]).astype(np.int64)
    return verts, faces


def _normv(x):
    return np.sqrt(np.sum(x ** 2, axis=1))


def cot_laplacian(X, T):
    """Cotangent Laplacian W and mixed-Voronoi vertex areas (cotLaplacian.m)."""
    nv = X.shape[0]
    L1 = _normv(X[T[:, 1]] - X[T[:, 2]])
    L2 = _normv(X[T[:, 0]] - X[T[:, 2]])
    L3 = _normv(X[T[:, 0]] - X[T[:, 1]])
    EL = np.stack([L1, L2, L3], axis=1)
    A1 = (L2 ** 2 + L3 ** 2 - L1 ** 2) / (2 * L2 * L3)
    A2 = (L1 ** 2 + L3 ** 2 - L2 ** 2) / (2 * L1 * L3)
    A3 = (L1 ** 2 + L2 ** 2 - L3 ** 2) / (2 * L1 * L2)
    ang = np.arccos(np.clip(np.stack([A1, A2, A3], axis=1), -1, 1))

    I = np.concatenate([T[:, 0], T[:, 1], T[:, 2]])
    J = np.concatenate([T[:, 1], T[:, 2], T[:, 0]])
    S = 0.5 / np.tan(np.concatenate([ang[:, 2], ang[:, 0], ang[:, 1]]))
    W = sp.csr_matrix(
        (np.concatenate([-S, -S, S, S]),
         (np.concatenate([I, J, I, J]), np.concatenate([J, I, I, J]))),
        shape=(nv, nv))

    # mixed Voronoi areas, barycentric fallback for obtuse triangles
    cA = 0.5 / np.tan(ang)
    vp1, vp2 = [1, 2, 0], [2, 0, 1]
    At = 0.25 * (EL[:, vp1] ** 2 * cA[:, vp1] + EL[:, vp2] ** 2 * cA[:, vp2])
    # NB: as in the MATLAB code, Ar is |cross| (= twice the triangle area)
    Ar = _normv(np.cross(X[T[:, 0]] - X[T[:, 1]], X[T[:, 0]] - X[T[:, 2]]))
    for c in range(3):
        locs = cA[:, c] < 0
        w = np.full(3, 1 / 8.0)
        w[c] = 1 / 4.0
        At[locs] = Ar[locs, None] * w[None, :]

    A = np.bincount(np.concatenate([T[:, 0], T[:, 1], T[:, 2]]),
                    weights=np.concatenate([At[:, 0], At[:, 1], At[:, 2]]),
                    minlength=nv)
    return W, A


def triangle_areas(X, T):
    return 0.5 * _normv(np.cross(X[T[:, 1]] - X[T[:, 0]], X[T[:, 2]] - X[T[:, 0]]))


def vertex_normals(X, T):
    fn = np.cross(X[T[:, 1]] - X[T[:, 0]], X[T[:, 2]] - X[T[:, 0]])  # area weighted
    vn = np.zeros_like(X)
    for c in range(3):
        np.add.at(vn, T[:, c], fn)
    return vn / np.maximum(np.linalg.norm(vn, axis=1, keepdims=True), 1e-12)


def laplacian_basis(X, T, k):
    W, a = cot_laplacian(X, T)
    M = sp.diags(a)
    evals, evecs = eigsh(W, k=k, M=M, sigma=-1e-6, which="LM")
    order = np.argsort(np.abs(evals))
    return np.abs(evals[order]), evecs[:, order], a


def geodesic_matrix(X, T, method="heat"):
    """All-pairs geodesic distances (compute_geodesic_dist_matrix.m).

    'heat' uses the heat method from potpourri3d (the MATLAB code uses fast
    marching); 'dijkstra' is a cheaper approximation along mesh edges.
    """
    n = X.shape[0]
    if method == "heat":
        import potpourri3d as pp3d

        solver = pp3d.MeshHeatMethodDistanceSolver(X, T)
        D = np.empty((n, n), dtype=np.float32)
        for i in range(n):
            D[:, i] = solver.compute_distance(i)
    elif method == "dijkstra":
        from scipy.sparse.csgraph import dijkstra

        E = np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]])
        w = np.linalg.norm(X[E[:, 0]] - X[E[:, 1]], axis=1)
        G = sp.csr_matrix((w, (E[:, 0], E[:, 1])), shape=(n, n))
        G = G.maximum(G.T)
        D = dijkstra(G, directed=False).astype(np.float32)
    else:
        raise ValueError(method)
    D = 0.5 * (D + D.T)
    np.maximum(D, 0, out=D)
    np.fill_diagonal(D, 0)
    return D


class Shape:
    """Container mirroring the MATLAB `shapes{i}` struct."""

    def __init__(self, verts, faces, name=""):
        self.name = name
        self.verts = verts
        self.faces = faces
        self.nv = verts.shape[0]

    @classmethod
    def load(cls, path, num_eigs=500, geodesic="heat", cache_dir=None,
             shot_bins=10, shot_radius=5.0):
        from .shot import shot_descriptor

        name = os.path.splitext(os.path.basename(path))[0]
        cache = None
        if cache_dir is not None:
            os.makedirs(cache_dir, exist_ok=True)
            cache = os.path.join(cache_dir, f"{name}_k{num_eigs}_{geodesic}.npz")
        verts, faces = read_ply(path)
        S = cls(verts, faces, name)
        if cache is not None and os.path.exists(cache):
            d = np.load(cache)
            S.evals, S.evecs, S.area_v, S.Gamma = d["evals"], d["evecs"], d["area_v"], d["Gamma"]
        else:
            # LB basis and geodesics are computed on the raw mesh, as in demo.m
            S.evals, S.evecs, S.area_v = laplacian_basis(verts, faces, num_eigs)
            S.Gamma = geodesic_matrix(verts, faces, geodesic)
            if cache is not None:
                np.savez(cache, evals=S.evals, evecs=S.evecs, area_v=S.area_v, Gamma=S.Gamma)
        S.sqrt_area = np.sqrt(triangle_areas(verts, faces).sum())

        # surfaceNorm.m: centre and scale the vertices, then compute SHOT
        Xn = verts - verts.mean(axis=0)
        Xn = Xn / np.sqrt(np.sum(Xn ** 2) / S.nv)
        S.verts_normalized = Xn
        area = triangle_areas(Xn, faces).sum()
        S.shot = shot_descriptor(Xn, vertex_normals(Xn, faces),
                                 shot_bins, shot_radius * np.sqrt(area) / 100)
        return S
