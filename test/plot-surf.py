from mayavi import mlab
from mayavi.modules.surface import Surface
from vtk.util import numpy_support as VN
import vtk
from numpy import *

fname = 'vecplot.vtk'
fig = mlab.figure()
engine = mlab.get_engine()

vtk_file_reader = engine.open(fname)

reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(fname)
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()
output = reader.GetOutput()
u = VN.vtk_to_numpy(output.GetPointData().GetArray('vec_comps'))
points = array(reader.GetOutput().GetPoints().GetData())
surface = Surface()
engine.add_filter(surface,vtk_file_reader)
mlab.quiver3d(points[:,0],points[:,1],points[:,2],u[:,0],u[:,1],u[:,2])


mlab.show()
