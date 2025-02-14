# Using Cascade in Quantum Key Distribution

This is a public version of the code used in *[Using Cascade in Quantum Key Distribution](https://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.20.064040)* \[[arxiv](https://arxiv.org/pdf/2307.00576)\]. This was built for [v2.0.2](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/releases/tag/v2.0.2) of the Open QKD Security package.

The code computes key rates and plots the figures from the paper. To plot Fig 1 and 2, run `Cascade/qubitBB84/main.m`. To plot Fig 3 and 4, run `Cascade/DecoyBB84/main.m`. For Fig 1 and Fig 2, set `replaceProb=0` and `replaceProb=0.2` respectively, for Fig 2 and 4 set `perturbProb=0.0` and `perturbProb=0.2` respectively..


## Installation instructions
> [!CAUTION]
> This repository is for archival and transparency purposes; we do not guarantee compatibility with other versions of the Open QKD Security package beyond the ones listed above.

### As zip
1. Download the linked version of the code from above and follow all [installation instructions](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/tree/06a2a37bffd8e207a0803ebb3737a13bb7b4f1fd).
2. Also follow the additional Mosek install instructions if you want an exact match.
3. Download the latest release on the side bar and unzip in your preferred directory and add this folder to the Matlab path.


### with git
1. Clone this repository and its exact submodules navigate to your desired directory and run,
```
git clone --recurse-submodules https://github.com/Optical-Quantum-Communication-Theory/Using-Cascade-in-Quantum-Key-Distribution
```
2. Follow all further [installation instructions](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/tree/06a2a37bffd8e207a0803ebb3737a13bb7b4f1fd).
3. Also follow the additional Mosek install instructions if you want an exact match.
