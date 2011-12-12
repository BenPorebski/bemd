import bemd
from vector3d import Vector3d, pos_distance
import os


def goto_dir(d):
  if not os.path.isdir(d):
    os.makedirs(d)
  os.chdir(d)


def nitrogen():
  goto_dir('nitrogen')
  el_mass = 2
  out_name = 'nitrogen'
  p = bemd.Particle(
      +6, 18.0, 500, 0.0, Vector3d(0.0, 0.0, 0.0))
  particles = [p]
  bemd.add_spherical(particles, 1, 0.01, 1, el_mass)
  bemd.add_spherical(particles, 2, 0.3, 2, el_mass)
  bemd.add_spherical(particles, 3, 0.9, 3, el_mass)
  bemd.run_md(
      particles, out_name, 20000, 5, 0.001, 
      [bemd.coulomb, bemd.lj], n_save=20)
  goto_dir('..')


def carbon():
  goto_dir('carbon')
  el_mass = 20
  out_name = 'carbon'
  p = bemd.Particle(
      +5, 18.0, 500, 0.0, Vector3d(0.0, 0.0, 0.0))
  particles = [p]
  bemd.add_spherical(particles, 1, 0.01, 0.2, el_mass)
  bemd.add_spherical(particles, 4, 0.9, 2, el_mass)
  bemd.run_md(
      particles, out_name, 8000, 1, 0.001, 
      [bemd.coulomb, bemd.lj], n_save=50)
  goto_dir('..')


nitrogen()
# carbon()
# os.system('python vmdtraj.py carbon/carbon')