# Trait gradients

Ideas and code for the analysis


## Installation

Installation of `plant` requires a C++11 compatible C compiler (OSX >= 10.10/Yosemite satisfies this, plus we've had success on Ubuntu 12.04 and 14.04).

The following can be installed directly in R

- remake: `devtools::install_github("richfitz/remake")`
- loggr: `devtools::install_github("smbache/loggr")`
- grader: `devtools::install_github("richfitz/grader")`


These packages are installed from github:

1. plant, development branch

```
git clone git@github.com:traitecoevo/plant.git
cd plant
git checkout development
make -j4
make install
```

2. plant.assembly

```
git clone git@github.com:traitecoevo/plant.assembly.git
cd plant.assembly
make install
```
