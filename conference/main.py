from dolfinx.fem import functionspace
from helezon.parameter_utils import Q_volumetric
from helezon.problem import SteadyState
from helezon.boundary_conditions import BoundaryCondition
from helezon.io_utils import xdmf_writer, XDMFReader, vtk_writer
from ufl import TrialFunction
import numpy as np
import params

geometry = XDMFReader("MeshDir/MIKRO_UYDU_VC2")
geometry.getInfo()

mesh, subdomains, facet_tags = geometry.getAll()

degree = 1

V = functionspace(mesh, ("CG", degree))
u = TrialFunction(V)

g = 1.88*1E-6

t_ref = 3.85*1E3
# Define the boundary conditions
# boundary_conditions = [
#     BoundaryCondition("DirichletValue", 656, T_external, V, facet_tags, u),
#     BoundaryCondition("DirichletValue", 621, T_external, V, facet_tags, u),
# ]

# boundary_conditions = [
#     BoundaryCondition("DirichletValue", 328, 38.8535, V, facet_tags, u),
#     BoundaryCondition("Neumann", 656, g, V, facet_tags, u),
#     BoundaryCondition("Neumann", 621, g, V, facet_tags, u),
# ]

x_ref = np.array([-0.17669, -0.0019525, 0.125864])
sigma = 0.25
T_external = 37.69
adjuster = T_external *2
def u_panel(x):
    spatial = (x[0]-x_ref[0])**2 + (x[1]-x_ref[1])**2 + (x[2]-x_ref[2])**2
    n =3
    spatial_term = np.exp(-1*spatial/(2*sigma**2))

    return T_external*spatial_term


boundary_conditions = [
    # BoundaryCondition("Neumann", 656, g, V, facet_tags, u),
    BoundaryCondition("DirichletGaussian", 621, u_panel, V, facet_tags, u),
    BoundaryCondition("DirichletGaussian", 656, u_panel, V, facet_tags, u),
]



Q = Q_volumetric(mesh, subdomains, Q_total=params.Q_total, tag=4, degree=0)

problem = SteadyState(V, subdomains, boundary_conditions, params.kappa, u, Q)
T = problem.solution

xdmf_writer("ResultsDir/T_steady", mesh, T)
# vtk_writer("ResultsDir/T_steady", mesh, T)
