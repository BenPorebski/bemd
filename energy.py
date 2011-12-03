import math
import random
import vector3d


def maxwell_velocity(temp, mass):
  "temp in K, mass in a.m.u., output in angstroms/ps"
  factor = 8314.47148 * temp / mass # in (m/s)^2
  convert_to_ang_per_ps = 0.01 # from m/s to angstroms/ps
  return random.gauss(0, math.sqrt(factor)) * convert_to_ang_per_ps


def mean_energy(temp, n_degree_of_freedom):
  boltzmann = 8314.47148  # Da (m/s)^2 K^-1
  convert_to_ang_per_ps = 1.0E-4 # (m/s)^2 to (angstroms/ps)^2
  return n_degree_of_freedom * 0.5 * boltzmann * \
         temp * convert_to_ang_per_ps # Da (angstroms/ps)^2


def random_energy(temp, n_degree_of_freedom):
  average = mean_energy(temp, n_degree_of_freedom);
  std_dev = math.sqrt(average)
  return random.gauss(average, std_dev)

  
def kinetic_energy(atoms):
  en = 0.0
  for a in atoms:
    vel = a.vel.length()
    en += 0.5 * a.mass * vel * vel
  return en


def gas_randomize(atoms, temp):
  "Heats residues uses gas approximation: vel: angstroms/ps"
  for atom in atoms:
    atom.vel.x = maxwell_velocity(temp, atom.mass)
    atom.vel.y = maxwell_velocity(temp, atom.mass)
    atom.vel.z = maxwell_velocity(temp, atom.mass)


def velocity_scale(atoms, temp, n_degree_of_freedom):
  "Scales velocity of atoms to energy at temp. Vel: angstroms/ps"
  target_energy = mean_energy(temp, n_degree_of_freedom)
  kin = kinetic_energy(atoms)
  if vector3d.is_near_zero(kin):
    gas_randomize(atoms, temp)
  else:
    scaling_factor = math.sqrt(target_energy / kin)
    for atom in atoms:
      atom.vel.scale(scaling_factor)


