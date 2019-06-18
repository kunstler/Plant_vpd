# Trait gradients

Explore how trait-based tradeoff can explain simultanously the variation of mean traits and the large range of coexisting traits along aridity graident using the model plant.


## Plant installation

Installation of `plant` requires a C++11 compatible C compiler (OSX >= 10.10/Yosemite satisfies this, plus we've had success on Ubuntu 12.04 and 14.04).

The following can be installed directly in R

- remake: `devtools::install_github("richfitz/remake")`
- loggr: `devtools::install_github("smbache/loggr")`
- grader: `devtools::install_github("richfitz/grader")`

For grail see https://github.com/traitecoevo/grail

install https://nlopt.readthedocs.io/en/latest/

on unbuntu `sudo apt-get install libnlopt-dev`

```
git clone https://github.com/traitecoevo/grail.git
cd grail
make install
```
- RcppAramadillo `install.packages("RcppArmadillo")`
- digest `install.packages("digest")`


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
git clone https://kunstler@github.com/traitecoevo/plant.assembly.git
cd plant.assembly
make install
```


### Hypothesis to test

- Decide whether we keep the initial parameterisation of plant photosynthesis and the test of the effect of change in the intercept or the slope of the lma *vs.* leaf lifespan tradeoff. 

1. trade-off between lma and ll change in slope and elevation. Wright
   et al.(2005) reported that the trade-off between LMA and Leaf Turn over Rate, LTR (1/LLS), is changing with aridity with higher aridity sites having higher LTR at a given LMA and a shallower slope.

- Or focus on the new parametrisation of the photosynthesis with aridity represented by vpd (no soil water stress effect) to test

1. change in Amax with vpd and only LMA assembly

2. vpd and change in assembly of both lma and Narea (Narea alone also?). Leaf N. Wright et al. (2003) explore how N can substitute for
   water. Aarea ∝ Narea gs (with gs stomatal conductance).

3. Narea also influence leaf lifespan because this reduce leaf mechanical strength.
   
### TODO

1. Fix error with FvC model version

2. Plot partial results

3. run with only LMA or Narea

4. Implement link between Narea and leaflifespan.
