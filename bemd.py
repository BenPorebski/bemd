import string
import math
import struct
import copy
import types
import tempfile
import os
import sys
import time

import energy
from vector3d import Vector3d, dot, pos_distance


class Particle:
  def __init__(self, q, m, lj_e, lj_r, pos, vel=Vector3d(0.0, 0.0, 0.0)):
    self.pos = pos
    self.vel = vel
    self.mass = m
    self.q = q
    self.lj_e = lj_e
    self.lj_r = lj_r
    self.force = Vector3d()


def write_psf(psf, particles):
  remarks = ['HA HA', 'HA HA']   
  f = open(psf, 'w')
  f.write('PSF\n')
  f.write('\n')
  f.write('       %d !NTITLE\n' % len(remarks))
  for remark in remarks:
    f.write('REMARKS %s\n' % remark)
  f.write('\n')
  f.write('    %d !NATOM\n' % len(particles))
  for i, particle in enumerate(particles):
    if particle.q > 0:
      name = "C"
    else:
      name = "H"
    s = ""
    s += "%8d " % i 
    s += "A " # chain
    s += "%4d " % 1 # i_res
    s += " %3s " % name # resname
    s += "%3s " % name  # atom name
    s += "%3s " % name  # atom type
    s += "   % 10.6f "  % particle.q
    s += "   % 10.4f "  % particle.mass
    s += "          0"
    s += "\n"
    f.write(s)
  f.write("\n       0 !NBOND: bonds\n\n")
  f.write("\n       0 !NTHETA: angles\n\n")
  f.write("\n       0 !NPHI: dihedrals\n\n")
  f.write("\n       0 !NIMPHI: impropers\n\n")
  f.close()


class DCDWrite:
  def __init__(self, dcd_file, n_atom):
    self.n_frame = 0  # Number of frames
    self.time_offset = 0  # time offset for trajectory
    self.n_save = 0
    self.n_fixed_atom = 0 # won't ever deal with fixed atoms
    self.dt = 0.00
    self.n_atom = n_atom
    remarks = ['REMARKS LINE 1', 'REMARKS LINE 2']
    self.remarks = [string.ljust(r, 80)[:80] for r in remarks]
    self.n_remark_line = len(self.remarks)
    self.size_remark = struct.calcsize('i')+80*self.n_remark_line

    # write header
    pack = struct.pack
    self.f = open(dcd_file, 'w+b')
    self.f.write(pack('i4c', 84, 'C', 'O', 'R', 'D'))
    
    self.pos_n_frame = self.f.tell()
    self.f.write(pack('i', self.n_frame))
    
    self.f.write(pack('i', self.time_offset))
    self.f.write(pack('i', self.n_save))
    self.f.write(pack('4i', 0, 0, 0, 0))
    self.f.write(pack('i', 0)) # Why?
    self.f.write(pack('i', self.n_fixed_atom))
    self.f.write(pack('d', self.dt))
    self.f.write(pack('i', 0)) # Why?
    self.f.write(pack('8i', 0, 0, 0, 0, 0, 0, 0, 0))
    self.f.write(pack('i', 84))
    self.f.write(pack('i', self.size_remark))
    self.f.write(pack('i', self.n_remark_line))
    for r in self.remarks:
      self.f.write(pack(*(['80c']+list(r))))
    self.f.write(pack('i', self.size_remark))
    self.f.write(pack('i', struct.calcsize('i')))
    self.f.write(pack('i', self.n_atom))
    self.f.write(pack('i', struct.calcsize('i')))
    self.f.flush()

  def append_with_numpy(self, particles):
    if 'numpy' not in globals():
      import numpy
    pack = struct.pack
    self.n_frame += 1
    self.f.seek(self.pos_n_frame)
    self.f.write(pack('i', self.n_frame))
    self.f.seek(0, 2) # go to end of file

    x = numpy.array([p.pos.x for p in particles])
    y = numpy.array([p.pos.y for p in particles])
    z = numpy.array([p.pos.z for p in particles])
    size_coordinate_block = struct.calcsize(repr(len(x))+'f')

    self.f.write(pack('i', size_coordinate_block))
    self.f.write(x.tostring())
    self.f.write(pack('i', size_coordinate_block))
    self.f.write(pack('i', size_coordinate_block))
    self.f.write(y.tostring())
    self.f.write(pack('i', size_coordinate_block))
    self.f.write(pack('i', size_coordinate_block))
    self.f.write(z.tostring())
    self.f.write(pack('i', size_coordinate_block))
    self.f.flush()
  
  def append(self, particles):
    pack = struct.pack
    self.n_frame += 1
    self.f.seek(self.pos_n_frame)
    self.f.write(pack('i', self.n_frame))
    self.f.seek(0, 2) # go to end of file

    x_coords = [p.pos.x for p in particles]
    y_coords = [p.pos.y for p in particles]
    z_coords = [p.pos.z for p in particles]
    size_coordinate_block = struct.calcsize("%df" % len(x_coords))
    len_str = '%df' % len(x_coords)

    self.f.write(pack('i', size_coordinate_block))
    self.f.write(pack(*([len_str] + x_coords)))
    self.f.write(pack('i', size_coordinate_block))
    self.f.write(pack('i', size_coordinate_block))
    self.f.write(pack(*([len_str] + y_coords)))
    self.f.write(pack('i', size_coordinate_block))
    self.f.write(pack('i', size_coordinate_block))
    self.f.write(pack(*([len_str] + z_coords)))
    self.f.write(pack('i', size_coordinate_block))
    self.f.flush()

  def close(self):
    self.f.close()


def coulomb(p1, p2, r):
  if p1.q*p1.q < 0 and r < 0.05:
    return 0
  return 331*p1.q*p2.q/(r*r)


def coulomb_energy(p1, p2, r):
  if r >= 1.0:
    return -331*p1.q*p2.q/r
  else:
    return 0.0


def rhooke(p1, p2, r):
  k = 1E8
  if p1.lj_r == 0.0 or p2.lj_r == 0.0:
    return 0
  lj_r = 0.5*(p1.lj_r + p2.lj_r)
  if r > lj_r:
    return 0.0
  return k*(lj_r-r)


def lj(p1, p2, r):
  if p1.q*p2.q < 0:
    return 0
  r0 = 0.5*(p1.lj_r + p2.lj_r)
  if r>r0:
    return 0
  e = math.sqrt(p1.lj_e*p2.lj_e)
  x = r0/r
  x_7 = pow(x, 7)
  x_13 = pow(x, 13)
  return 12.0*e*(x_13 - x_7) 


def lj_energy(p1, p2, r):
  r0 = 0.5*(p1.lj_r + p2.lj_r)
  if r>r0:
    return 0
  e = math.sqrt(p1.lj_e*p2.lj_e)
  x = r0/r
  x_6 = pow(x, 6)
  x_12 = pow(x, 12)
  return e*(x_12 - 2.0*x_6) + e


def calculate_forces(particles, center_forces):
  for p in particles:
    p.force.set(0.0, 0.0, 0.0)
  n_particle = len(particles)
  for i in range(n_particle):
    for j in range(i+1, n_particle):
      disp_i_to_j = particles[j].pos - particles[i].pos
      r = disp_i_to_j.length()
      dir_i_to_j = disp_i_to_j.scaled_vec(1.0/r)
      for center_force_fn in center_forces:
        f_mag = center_force_fn(particles[i], particles[j], r)
        f_mag_on_j = dir_i_to_j.scaled_vec(f_mag)
        particles[j].force += f_mag_on_j
        particles[i].force += -f_mag_on_j


def apply_forces(particles, dt):
  for p in particles:
    p.pos += p.vel.scaled_vec(dt) 
    p.vel += p.force.scaled_vec(dt/p.mass)
 

def velocity_scale(atoms, temp, n_degree_of_freedom):
  "Scales velocity of atoms to energy at temp. Vel: angstroms/ps"
  target_energy = energy.mean_energy(temp, n_degree_of_freedom)
  for atom in atoms:
    kin = energy.kinetic_energy([atom])
    scaling_factor = math.sqrt(target_energy / kin)
    atom.vel.scale(scaling_factor)


def run_md(
    particles, out_name, n_step, temp, 
    dt, center_forces, is_elastic=True, n_save=200):
  psf = out_name + '.psf'
  dcd = out_name + '.dcd'
  write_psf(psf, particles)
  dcd_writer = DCDWrite(dcd, len(particles))
  energy.gas_randomize(particles, temp)
  av_energy = energy.mean_energy(temp, 3*len(particles))
  av_vel = math.sqrt(av_energy*2/el_mass)
  av_dist = av_vel*dt
  print "Temp:%.3f, Energy:%.4f, Velocity:%.4f, Dist:%f" % \
      (temp, av_energy, av_vel, av_dist)
  for i_step in range(n_step):
    calculate_forces(particles, center_forces)
    apply_forces(particles, dt)
    velocity_scale(particles, temp, 3*len(particles))
    if i_step % n_save == 0:
      dcd_writer.append(particles)
      if i_step % (n_save*100) == 0:
        print i_step // n_save, "frames"
  dcd_writer.close()


def print_distances(particles):
  dists = []
  n_particle = len(particles)
  for i in range(n_particle):
    for j in range(i+1, n_particle):
      d = pos_distance(particles[i].pos, particles[j].pos)
      dists.append((d, i, j))
  dists.sort(reverse=True)
  for d in dists:
    print "%4.1f %d %d" % d


def generate_sphere_points(n):
  """
  Returns list of 3d coordinates of points on a sphere using the
  Golden Section Spiral algorithm.
  """
  points = []
  inc = math.pi * (3 - math.sqrt(5))
  offset = 2 / float(n)
  for k in range(int(n)):
    y = k * offset - 1 + (offset / 2)
    r = math.sqrt(1 - y*y)
    phi = k * inc
    points.append((math.cos(phi)*r, y, math.sin(phi)*r))
  return points


def add_spherical(particles, n, r, scale_r):
  for point in generate_sphere_points(n):
    point = [scale_r*r*p for p in point]
    particles.append(
        Particle(-1.0, el_mass, 500, r, Vector3d(*point)))


el_mass = 50.0
out_name = 'nitrogen'
r1 = 0.05
r2 = 0.4
r3 = 0.9
r5 = 2
n = 5
particles = [Particle(+n, 18.0, 500, 0.0, Vector3d(0.0, 0.0, 0.0))]
add_spherical(particles, 1, r1, 1)
# add_spherical(particles, 2, r2, 2)
add_spherical(particles, 4, r3, 2)
# particles = [
#    Particle(+1.0, 18.0, 500, 2.0, Vector3d(0.0, 0.0, 0.0)),
#    Particle(-1.0, 18.0, 500, 2.0, Vector3d(0.0, 4.2, 0.0))]
run_md(
    particles, out_name, 10000, 3, 0.001, 
    [coulomb, lj], n_save=20)
print_distances(particles)
os.system('python vmdtraj.py ' + out_name)

