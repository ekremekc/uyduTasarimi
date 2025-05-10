from dolfinx.fem import functionspace
from helezon.parameter_utils import Q_volumetric
from helezon.problem import SteadyState
from helezon.boundary_conditions import BoundaryCondition
from helezon.io_utils import xdmf_writer, XDMFReader
from ufl import TrialFunction
import params

geometry = XDMFReader("MeshDir/MIKRO_UYDU_VC")
geometry.getInfo()

mesh, subdomains, facet_tags = geometry.getAll()

degree = 1

V = functionspace(mesh, ("Lagrange", degree))
u = TrialFunction(V)

# Define the boundary conditions
boundary_conditions = [
    BoundaryCondition("DirichletValue", 656, 0, V, facet_tags, u),
]

Q = Q_volumetric(mesh, subdomains, Q_total=params.Q_total, tag=4, degree=0)

problem = SteadyState(V, subdomains, boundary_conditions, params.kappa, u, Q)
T = problem.solution

xdmf_writer("ResultsDir/T_steady", mesh, T)
