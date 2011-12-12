import struct


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


