from mayavi import mlab
from mayavi.modules.surface import Surface

fname = 'sphere.vtk'
fig = mlab.figure()
engine = mlab.get_engine()

vtk_file_reader = engine.open(fname)

surface = Surface()
engine.add_filter(surface,vtk_file_reader)

mlab.show()
