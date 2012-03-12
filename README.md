
# Bare Electron Molecular Dynamics

Bare Electron Molecular Dynamics is a program to simulate molecules using the hard electron approximation for atomic representation.

This level of approximation uses only the electrostatic interaction between nuclei and electrons as the only attractive term. This even allows a limited modeling of charge balance to model the stability of different atoms, even metals. 

An artificial repulsion term is used to model the Pauli exclusion principle. Once you explicitly include electrons, you can dispense with the covalent harmonic springs used in typical molecular dynamics simulations. We will investigate different repulsive functions such as step functions and non-interacting electron spin species.

BEMD therefore allows the creation and destruction of chemical bonds. Since we are interested in configurations, rather than kinetics, we use a hard thermometer to scale velocities to a regime where the electrons can find the optimal configuration. The temperature slows down the effect of the electrostatics, and allows a comprehensive stochastic search for the best electron configuration.

To test the feasbility of BEMD, we hope to show that this representation can generate stable configurations of the following molecules:

## First small steps

- H2
- methane CH4
- ethane CH2=CH2
- CO
- CO2
- N2
- H20
- NO
- methanol
- pyrimidine
- benzene
- acetylmide
- amide
- acetic acid
- ribose

## Other moleucles
- B2H6
- alanine dipeptide
- ATP
- water
- CO2
- benzene
- urea
- small-disulfide protein
- oxygen molecule - paramagnetism

## Ultimate molecules

- halo addition
- SN2 reactions
- diels alder reaction 
- proteolysis in solution (pH)
- ATP hydrolysis
- ion channel discrimination
- water wire for hydrogen transfer
- redox reaction - copper zinc
- semiconductor in crystal geometry
- urey-miller experiment

(c) 2011, Ben Porebski, Mike Kuiper and Bosco Ho.