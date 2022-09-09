# from src.Dirichlet.Initial_deformation.One_domain import loop_sin_deformation
from src.Dirichlet.Initial_deformation.One_domain import sin_deformation
# from src.Neumann.Initial_deformation.One_domain import sin_defformation
from src.Graph import graph
# from src.Dirichlet.Initial_deformation.subdomains import *
# from src.Dirichlet.Input_pulse.One_domain import *
# from src.Dirichlet.Input_pulse.subdomains import *
# from src.Neumann.Initial_deformation.One_domain import *
# from src.Neumann.Initial_deformation.subdomains import *
# from src.Neumann_Dirichlet.Input_pulse import *


num_steps = 60
dt = 1e-4
nx = ny = 500
c2 = 1

name = f'nsteps_{num_steps}_dt{dt}_nx{nx}_ny_{ny}_c_{c2}_cos'
path = 'src/Dirichlet/Initial_deformation/One_domain/results/'
# loop_sin_deformation.loop_sin_deformation(num_steps, dt, nx, ny, c2)
sin_deformation.sinusoidal_deformation(num_steps, dt, nx, ny, c2, name, path)
graph.graph(path, name, num_steps, step = 1, save=True, plot=False)
