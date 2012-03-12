import bemd
from vector3d import Vector3d, pos_distance
import os

def goto_dir(d):
  if not os.path.isdir(d):
    os.makedirs(d)
  os.chdir(d)

goto_dir('carbon')
out_name = 'carbon'

el_mass = 2
particles = []
particles.append(bemd.Particle(+4, 18, 500, 0.05, Vector3d(0.0, 0.0, 0.0)))
for v in bemd.spherical_vecs(4, 3):
  particles.append(bemd.Particle(-1, el_mass, 500, 0.9, v))

bemd.run_md(
  particles, out_name, 8000, 1, 0.001, 
  [bemd.coulomb, bemd.lj], n_save=50)
goto_dir('..')

os.system('python vmdtraj.py carbon/carbon')