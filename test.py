import bemd
import vector3d

def print_distances(particles):
  dists = []
  n_particle = len(particles)
  for i in range(n_particle):
    for j in range(i+1, n_particle):
      d = vector3d.pos_distance(particles[i].pos, particles[j].pos)
      dists.append((d, i, j))
  dists.sort(reverse=True)
  for d in dists:
    print "%4.1f %d %d" % d


el_mass = 2
out_name = 'nitrogen'
particles = [bemd.Particle(+4, 18.0, 500, 0.0, Vector3d(0.0, 0.0, 0.0))]
bemd.add_spherical(particles, 1, 0.01, 1)
# add_spherical(particles, 2, 0.4, 2)
bemd.add_spherical(particles, 3, 0.9, 3)
bemd.run_md(
    particles, out_name, 80000, 5, 0.001, 
    [coulomb, lj], n_save=20)
print_distances(particles)
os.system('python vmdtraj.py ' + out_name)

