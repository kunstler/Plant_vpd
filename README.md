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


### TODO

- implement new parametrisation to include aditional effect along
  climatic gradients

1. trade-off between lma and ll change in slope and elevation. Wright
   et al.(2005) reported that the trade-off between LMA and Leaf Turn over Rate, LTR (1/LLS), is changing with aridity with higher aridity sites having higher LTR at a given LMA and a shallower slope.

2. We could explore the effect of changes in elevation of the LLS-LMA
   trade- off along a precipitation gradient as this is one the
   stronger pattern in the Glopnet data. Which refs ??

3. Leaf N. Wright et al. (2003) explore how N can substitute for
   water. Aarea ∝ Narea gs (with gs stomatal conductance). One
   solution would be to assume that the water stress is represented by
   1/gs and Narea could off-set the effect of water stress on Aarea
   (in our case Amax). The cost of leaf N is in its respiration cost
   with higher dark leaf respiration per area for higher Narea. But
   need to compare with Prentice et al. (2014).

4. Frost gradients affect survival but short vessel size allow for an adaptation to lower temperature in term of survival but with a carbon cost.
