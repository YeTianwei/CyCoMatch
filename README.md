# Multi-Shape Matching with Cycle Consistency Basis via Functional Maps

Matlab and Python implementation for AAAI'25 paper "Multi-Shape Matching with Cycle Consistency Basis via Functional Maps", by Yifan Xia, Tianwei Ye, Huabing Zhou, Zhongyuan Wang and Jiayi Ma.

This paper focuses on multi-shape matching with a more concise cycle consistency formulation. By modeling shapes as directed graphs, we derived cycle consistency basis as constraints. The method employs a two-stage optimization: optimizing pairwise correspondence accuracy via landmarks and refining consistency through cycle consistency basis.

![teaser](./figs/teaser.jpg "teaser")

If you find this project useful, please cite:

```
@inproceedings{xia2025multi,
  title={Multi-Shape Matching with Cycle Consistency Basis via Functional Maps},
  author={Xia, Yifan and Ye, Tianwei and Zhou, Huabing and Wang, Zhongyuan and Ma, Jiayi},
  booktitle={Proceedings of the AAAI Conference on Artificial Intelligence},
  volume={39},
  number={8},
  pages={8575--8583},
  year={2025}
}
```

## Run the exemplary code

### Matlab

Please unzip the gspbox.zip and utils.zip in this directory. Run demo.m directly.

Note that this demo requires Matlab environment and GPU. (GPU accelerates the KNN, and users can convert to CPU.)

### Python

A Python implementation of the full pipeline is provided in [`python/`](./python). It requires neither Matlab nor a GPU; if PyTorch with CUDA is installed, it uses the GPU for the KNN searches automatically.

```bash
cd python
pip install -r requirements.txt
python demo.py                    # match the six FAUST shapes in ./data
python demo.py --plot result.png  # also save a colour-transfer figure
```

On the bundled data, the Python version brings the mean geodesic error (normalised by √area) from 0.1086 with the SHOT initialisation down to 0.0033. See [`python/README.md`](./python/README.md) for the module layout and the differences from the Matlab code.

## Contact

Please let us know (twye2001@gmail.com) if you have any question regarding the algorithms/paper or you find any bugs in the implementation.

## Acknowledgement
This implementation of functional map reconstruction is adapted from [MWP](https://github.com/Qinsong-Li/MWP) and the implementation of cycle consistency basis generation is modified from [NeurIPS19_Cycle](https://github.com/huangqx/NeurIPS19_Cycle). We thank the authors for making their codes publicly available.
