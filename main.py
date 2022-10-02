# from src.Dirichlet.Initial_deformation.One_domain import loop_sin_deformation
# from src.Dirichlet.Initial_deformation.One_domain import deformation
from src.Dirichlet.Initial_deformation.subdomains import deformation
# from src.Neumann.Initial_deformation.One_domain import sin_defformation
from src.Graph import graph
# from src.Dirichlet.Initial_deformation.subdomains import *
# from src.Dirichlet.Input_pulse.One_domain import *
# from src.Dirichlet.Input_pulse.subdomains import *
# from src.Neumann.Initial_deformation.One_domain import *
# from src.Neumann.Initial_deformation.subdomains import *
# from src.Neumann_Dirichlet.Input_pulse import *
from dolfin import *

num_steps = 100
dt = 1e-4
nx = ny = 100


name = f'two_domains_sin_nsteps_{num_steps}_dt_{dt}'
path = 'src/Dirichlet/Initial_deformation/subdomains/results/'
# loop_sin_deformation.loop_sin_deformation(num_steps, dt, nx, ny, c2)
# deformation.sinusoidal_deformation(num_steps, dt, nx, ny, c2, name, path)
# graph.graph(path, name, num_steps, step = 1, save=True, plot=False)
c0 = 100
c1 = 10
deformation.sin_deformation(dt, num_steps, name, path, c0, c1)
graph.graph(path, name, num_steps, step = 1, save=True, plot=False)