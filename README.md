# Wannier Orbital Overlap Population tools (WOOPs)

A post-processing tool written in python to get Wannier Orbital Overlap Population (WOOP), Wannier Orbital Hamiltonian Population (WOPP)* from [Wannier90](https://github.com/wannier-developers/wannier90) package.


Before getting into things, you might want to check out these two papers:
1.  [arXiv](https://arxiv.org/pdf/2009.01130.pdf)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

For the script to work, you need to have an valid installation of `python` (2.7.x).
Also, `numpy` package are needed, you can install them by pip:
```
pip install numpy
```
or by conda
```
conda install numpy
```
if you use a supercomputer and don't have enough privilege:

1. install anaconda by downloading from [here](https://www.anaconda.com/download/) and upload it to your directory.
2. using queue system to install anaconda by `chmod 755 anaconda*.sh && ./anaconda*.sh`
3. install `numpy` by download them from [here](https://anaconda.org/anaconda/numpy), upload them as well.
4. manually install package by `conda install numpy*.tar.bz2`.

## Useage
File need for WOOPs:

1. `wannier90_u.mat`
2. `wannier90.eig`

Detailed preparation for Wannier90 generated files can be found in Wanneir90's user guide, or in the `example` folder.

Please read `input.py` for more information.


## Running the tests

Go check the description in `example` folder.

## How to cite

For the method please cite the following paper in any publications arising from the use of this code:

  Sudipta Kundu,Satadeep Bhattacharjee, Seung-Cheol Lee and Manish Jain
  *Population Analysis with Wannier Orbitals*,[arXiv](https://arxiv.org/pdf/2009.01130.pdf)

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Chengcheng XIAO** - *Initial work* - [E-mail](iconxicon@me.com)

## License

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/Chengcheng-Xiao/WOBSTER) file for details
