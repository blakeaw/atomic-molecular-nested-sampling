# Nested Sampling of atomic and molecular systems.

This archived repository contains various codes for Nested Sampling simulations of atomic and molecular systems, including the variations of the method we developed:

  * Isobaric Nested Sampling (IBNS): Wilson et al. (2015) Nested sampling of isobaric phase space for the direct evaluation of the isothermal-isobaric partition function of atomic systems, J. Chem. Phys., 143, 154108; https://doi.org/10.1063/1.4933309
  * Coupling Parameter Path Nested Sampling (CPPNS): Wilson et al. (2018) Computing free energies using nested sampling-based approaches, Molecular Simulation, 44:13-14, 1108-1123; https://doi.org/10.1080/08927022.2017.1416113

Some cases include corresponding Metropolis Monte Carlo implementations that were used for comparison/validation. There are also many codes from unpublished work or which encode experimental nested sampling implementations. The codes were developed as a part of my research and training as a Ph.D. Student at UT Dallas (sometime between 2012 and 2016).

## Implementations

Inside the `src` directory there are the following  collections:

* `classes` - a collection of of custom class objects implemented for these simulations.
* `cg-3bead-lipids` - canonical nested sampling of systems with 3-bead coarse-grained lipids from Cooke et al. https://doi.org/10.1103/PhysRevE.72.011506  
* `CPPNS` - implementations of coupling parameter path nested sampling to estimate free energy perturbations and alchemical free energy changes.
* `grand-canonical` - implementations of constant chemical potential nested sampling to estimate the Grand Canonical partition function.   
* `IBNS` - implementations of isobaric nested sampling to estimate the isothermal-isobaric partition function.   
* `np-dimers-dumbbells` - Nested sampling simulations of nanoparticle dumbbells.
* `ns-replica-exchange` - experimental implementations of Nested Sampling replicas incorporating elements of Replica Exchange methods.
* `ns-test-area-method` - experimental implementations of Nested Sampling combined with the Test Area Method.
