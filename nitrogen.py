import os

import bemd
from vector3d import Vector3d

if not os.path.isdir('nitrogen'):
  os.makedirs('nitrogen')
os.chdir('nitrogen')

el_mass = 18
out_name = 'nitrogen'
particles = []
particles.append(bemd.Particle(+5, 18, 500, 0.05, Vector3d(0.0, 0.0, 0.0)))
for v in bemd.spherical_vecs(2, 0.3):
  particles.append(bemd.Particle(-1, el_mass, 500, 0.3, v))
for v in bemd.spherical_vecs(3, 3):
  particles.append(bemd.Particle(-1, el_mass, 500, 0.9, v))

bemd.run_md(
    particles, out_name, 20000, 2, 0.001, 
    [bemd.coulomb, bemd.lj], n_save=20)
os.chdir('..')

os.system('python vmdtraj.py nitrogen/nitrogen')