# CyCoMatch: Python implementation

This folder reimplements the MATLAB code in NumPy and SciPy. It needs no
MATLAB, mex files or GPU. If PyTorch with CUDA is installed, the nearest-neighbour
searches run on the GPU automatically.

```bash
pip install -r requirements.txt
python demo.py                       # match the six FAUST shapes in ../data
python demo.py --plot result.png     # also save a colour-transfer figure
```

On the first run, the Laplace–Beltrami basis and the geodesic matrices are cached in `./cache`.

## Module map

| MATLAB | Python |
| --- | --- |
| `demo.m` | `demo.py` |
| `src/Multi_CCB.m` | `cycomatch/multi_ccb.py` |
| `src/LMD_Filter.m`, `tools/findK.m` | `cycomatch/fmaps.py: lmd_filter, find_k` |
| `src/MWP_reconstruct.m`, `gsp_design_meyer.m` | `cycomatch/fmaps.py: mwp_reconstruct, meyer_filterbank` |
| `src/add_fmap.m` | `cycomatch/fmaps.py: add_fmap` |
| `tools/fmap2pmap_group.m` | `cycomatch/fmaps.py: fmap2pmap_group` |
| `src/admm_optimization.m` | `cycomatch/admm.py` |
| `external/cycle_basis_generator.m` (+ ADMM weight solver, sampling) | `cycomatch/cycles.py` |
| `MESH.preprocess`, `cotLaplacian.m`, `surfaceNorm.m` | `cycomatch/mesh.py` |
| `calc_shot` (mex) | `cycomatch/shot.py` (PCL SHOT352 port) |
| `knnsearch` (gpuArray) | `cycomatch/knn.py` |

## Differences from the MATLAB code

- **Geodesics** are computed with the heat method (`potpourri3d`). MATLAB uses
  fast marching. You can use `--geodesic dijkstra` for a faster, rougher approximation.
- **SHOT** is reimplemented from the PCL algorithm, so the descriptors (and
  therefore the initial maps `T0`) are close to `calc_shot.mex` but not bit-identical.
- **Landmark distortion.** `LMD_Filter.m` reads the target distances as
  `Dist{j}(idx_source, T(k))`, which indexes shape *j* with vertex ids from
  shape *i*. On FAUST every shape has the same vertex order, so this
  indexing quietly uses the ground-truth correspondence. By default the port uses the
  paper's form `d_j(T(m), T(k))`. Pass `--matlab-compat` (or
  `matlab_compat=True`) to reproduce the MATLAB behaviour exactly. On the
  bundled demo, both settings reach the same accuracy (see below).
- `fast_pMap_NNinterp` only fills `NaN` entries. A nearest-neighbour search never
  produces any, so the port omits it.
- In the cycle-consistency ADMM, MATLAB overwrites `Z2{i,j}` for every cycle that
  contains `i -> j`. The port evaluates only the last cycle for each pair, which gives
  the same result with fewer 500×500 products.
- Cycle sub-sampling uses NumPy's RNG (`seed` argument) instead of MATLAB's `rand`.
  With 6 shapes, every candidate cycle is kept, so this has no effect.

## Results on the bundled data

Six FAUST shapes, 30 ordered pairs. The table gives the mean geodesic error
normalised by √area; the ground truth is the identity map.

| Map | Mean error |
| --- | --- |
| SHOT initialisation `T0` | 0.1086 |
| Multi-CCB (default, paper distortion) | 0.0033 |
| Multi-CCB (`--matlab-compat`) | 0.0034 |

A full run takes about 25–30 minutes on a 4-core CPU without a GPU. Most of the
time goes to nearest-neighbour searches in up to 500 dimensions and to the 500×500
products in the ADMM step. The first run also spends about 2 minutes caching the
eigenbasis and geodesics.
