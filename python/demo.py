"""Python counterpart of demo.m: match the six FAUST shapes in ../data."""
import argparse
import glob
import os
import time

import numpy as np

from cycomatch import Shape, knnsearch, multi_ccb


def geo_error(shapes, T):
    """Mean geodesic error over all off-diagonal pairs, normalised by sqrt(area).

    The FAUST registrations share their vertex order, so the ground-truth map
    is the identity.
    """
    errs = []
    n = len(shapes)
    for i in range(n):
        for j in range(n):
            if i != j:
                Sj = shapes[j]
                gt = np.arange(shapes[i].nv)
                errs.append(np.mean(Sj.Gamma[T[i][j], gt]) / Sj.sqrt_area)
    return float(np.mean(errs))


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(here, "..", "data"))
    ap.add_argument("--cache", default=os.path.join(here, "cache"))
    ap.add_argument("--num-eigs", type=int, default=500)
    ap.add_argument("--geodesic", choices=["heat", "dijkstra"], default="heat")
    ap.add_argument("--matlab-compat", action="store_true",
                    help="reproduce the old (pre-fix) distortion indexing of LMD_Filter.m")
    ap.add_argument("--plot", default=None, help="save a colour-transfer figure to this path")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.data, "*.ply")))
    shapes = []
    for f in files:
        t = time.time()
        shapes.append(Shape.load(f, num_eigs=args.num_eigs, geodesic=args.geodesic,
                                 cache_dir=args.cache))
        print(f"loaded {os.path.basename(f)} in {time.time() - t:.1f}s")

    n = len(shapes)
    # initial maps from SHOT descriptors
    T0 = [[np.arange(shapes[i].nv) if i == j else knnsearch(shapes[j].shot, shapes[i].shot)
           for j in range(n)] for i in range(n)]
    G = 1 - np.eye(n, dtype=int)
    print(f"SHOT initialisation geodesic error: {geo_error(shapes, T0):.4f}")

    t = time.time()
    C, T = multi_ccb(shapes, T0, G, matlab_compat=args.matlab_compat)
    print(f"Multi-CCB done in {time.time() - t:.1f}s")
    print(f"Multi-CCB geodesic error: {geo_error(shapes, T):.4f}")

    if args.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(12, 7))
        V0 = shapes[0].verts
        rgb = (V0 - V0.min(0)) / (V0.max(0) - V0.min(0))
        for idx, S in enumerate(shapes):
            ax = fig.add_subplot(2, 3, idx + 1, projection="3d")
            ax.scatter(*S.verts.T, c=rgb[T[idx][0]], s=1)  # transfer colours from shape 1
            ax.set_axis_off()
            ax.view_init(elev=90, azim=-90)
        fig.savefig(args.plot, dpi=120, bbox_inches="tight")
        print(f"figure saved to {args.plot}")


if __name__ == "__main__":
    main()
