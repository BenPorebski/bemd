import math
import heating
import trajectory
from vector3d import Vector3d


"""
Units: 

mass: Da = g/mol or 1E-3 kg/mol
dist: Angstroms = 1E-10m
time: fs/T = 1/48.88821 1E-15s
T = 48.88821 is needed to scale 0.5*m*v^2 to kcal/mol

charge: e multiples of electron charge
energy: kcal/mol

thus coulomb prefactor is 332.0636 converts to kcal/mol

in other words we want 
force: kcal/mol/Angstrom
"""

COULOMB_FACTOR = 332.0636
COULOMB_CUTOFF = 0.05

class Particle:
  def __init__(
      self, q, m, lj_e, lj_r, pos, 
      vel=Vector3d(0.0, 0.0, 0.0)):
    self.pos = pos # angstrom
    self.vel = vel # angs/ps
    self.mass = m # dalton
    self.q = q # electron charge
    self.lj_e = lj_e # kcal/mol
    self.lj_r = lj_r # angstrom
    self.force = Vector3d()


def add_spherical(particles, n, r, scale_r, el_mass):
  """
  Returns list of 3d coordinates of points on a sphere using the
  Golden Section Spiral algorithm.
  """

  def generate_sphere_points(n):
    points = []
    inc = math.pi * (3 - math.sqrt(5))
    offset = 2 / float(n)
    for k in range(int(n)):
      y = k * offset - 1 + (offset / 2)
      r = math.sqrt(1 - y*y)
      phi = k * inc
      points.append((math.cos(phi)*r, y, math.sin(phi)*r))
    return points

  for point in generate_sphere_points(n):
    point = [scale_r*r*p for p in point]
    particles.append(
        Particle(-1.0, el_mass, 500, r, Vector3d(*point)))


def coulomb(p1, p2, r):
  "output in kcal/mol/angs"
  if p1.q*p1.q < 0 and r < COULOMB_CUTOFF:
    return 0
  return COULOMB_FACTOR*p1.q*p2.q/(r*r)


def coulomb_energy(p1, p2, r):
  "output in kcal/mol"
  if r >= 1.0:
    return -COULOMB_FACTOR*p1.q*p2.q/r
  else:
    return 0.0


def kinetic_energy_in_kcalmol(p):
  "output in kcal/mol"
  return 2.39E-3 * 0.5 * p.mass * p.vel.length()**2


def rhooke(p1, p2, r):
  "output in kcal/mol/angs"
  k = 1E8
  if p1.lj_r == 0.0 or p2.lj_r == 0.0:
    return 0
  lj_r = 0.5*(p1.lj_r + p2.lj_r)
  if r > lj_r:
    return 0.0
  return k*(lj_r-r)


def lj(p1, p2, r):
  "output in kcal/mol/angs"
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
  "output in kcal/mol"
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
 

header = """
============================================
Bare Electron Molecular Dynamics
(c) 2011 Ben Porebski, Mike Kuiper, Bosco Ho
============================================
"""

def run_md(
    particles, out_name, n_step, temp, time_step_ps, 
    center_forces, is_elastic=True, 
    n_save=200):
  """
  dt in ps, temp in K
  """
  print header
  psf = out_name + '.psf'
  dcd = out_name + '.dcd'
  trajectory.write_psf(psf, particles)
  dcd_writer = trajectory.DCDWrite(dcd, len(particles))
  heating.gas_randomize(particles, temp)
  dt = time_step_ps
  for i_step in range(n_step):
    calculate_forces(particles, center_forces)
    apply_forces(particles, dt)
    heating.strong_velocity_scale(
        particles, temp, 3*len(particles))
    if i_step % n_save == 0:
      dcd_writer.append(particles)
      if i_step % (n_save*100) == 0:
        print i_step // n_save, "frames"
  dcd_writer.close()


