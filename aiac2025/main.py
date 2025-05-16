from dolfinx.fem import functionspace, Function
from helezon.parameter_utils import Q_volumetric
from helezon.problem import SteadyState
from helezon.boundary_conditions import BoundaryCondition
from helezon.io_utils import xdmf_writer, XDMFReader, vtk_writer
from ufl import TrialFunction
from dolfinx import default_scalar_type
import numpy as np
import params
import sys

thermal_pad = params.thermal_pad
filename = "MIKRO_UYDU_VC_REORIENTED"

if thermal_pad:
    ind = sys.argv.index('-thickness')
    filename += sys.argv[ind+1]

geometry = XDMFReader("MeshDir/"+filename)
geometry.getInfo()

mesh, subdomains, facet_tags = geometry.getAll()

degree = 1

V = functionspace(mesh, ("CG", degree))
u = TrialFunction(V)

equipment_tag = 2
pad_tag = 3

t_ref = 3.85*1E3
x_ref = np.array([0.0, -0.05, 0])
sigma = 0.25
T_external = 37.69
adjuster = T_external *2
def u_panel(x):
    spatial = (x[0]-x_ref[0])**2 + (x[1]-x_ref[1])**2 + (x[2]-x_ref[2])**2
    n =3
    spatial_term = np.exp(-1*spatial/(2*sigma**2))

    return T_external*spatial_term

if thermal_pad:
    boundary_conditions = [
        BoundaryCondition("DirichletGaussian", 659,  u_panel, V, facet_tags, u),
        BoundaryCondition("DirichletGaussian", 697, u_panel, V, facet_tags, u),
    ]
    Q = functionspace(mesh, ("DG", 0))
    kappa = Function(Q)
    kappa.x.array[:] = params.kappa
    pad_cells = subdomains.find(pad_tag)
    kappa.x.array[pad_cells] = np.full_like(pad_cells, params.kappa_pad, dtype=default_scalar_type)

if not thermal_pad:
    boundary_conditions = [
        BoundaryCondition("DirichletGaussian", 619, u_panel, V, facet_tags, u),
        BoundaryCondition("DirichletGaussian", 656, u_panel, V, facet_tags, u),
    ]

    kappa = params.kappa

Q = Q_volumetric(mesh, subdomains, Q_total=params.Q_total, tag=equipment_tag, degree=0)

problem = SteadyState(V, subdomains, boundary_conditions, kappa, u, Q)
T = problem.solution

if thermal_pad:
    ind = sys.argv.index('-thickness')
    xdmf_writer("ResultsDir/T_pad_"+sys.argv[ind+1], mesh, T)
else:
    xdmf_writer("ResultsDir/T_nopad", mesh, T)
# vtk_writer("ResultsDir/T_steady", mesh, T)

from helezon.dolfinx_utils import getFunctionAverage

faceComponent = 688
faceChassis = 696
T_avg_chassis = getFunctionAverage(mesh, facet_tags, faceChassis, T, "intersection")
T_avg_component = getFunctionAverage(mesh, facet_tags, faceComponent, T, "intersection")
diff = T_avg_chassis-T_avg_component
print("difference = ", diff)

np.savetxt("PlotsDir/T_diff_"+sys.argv[ind+1]+".txt", np.array([float(sys.argv[ind+1]), diff]))