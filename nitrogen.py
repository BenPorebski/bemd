import bemd
import vector3d
import os

if not os.path.isdir('nitrogen'):
  os.makedirs('nitrogen')
os.chdir('nitrogen')
el_mass = 18
out_name = 'nitrogen'
p = bemd.Particle(
    +6, 18.0, 500, 0.0, vector3d.Vector3d(0.0, 0.0, 0.0))
particles = [p]
bemd.add_spherical(particles, 1, 0.01, 1, el_mass)
bemd.add_spherical(particles, 2, 0.3, 2, el_mass)
bemd.add_spherical(particles, 3, 0.9, 3, el_mass)
bemd.run_md(
    particles, out_name, 20000, 3, 0.001, 
    [bemd.coulomb, bemd.lj], n_save=20)
os.chdir('..')

os.system('python vmdtraj.py nitrogen/nitrogen')