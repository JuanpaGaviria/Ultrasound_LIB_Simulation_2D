from src.Neumann.Input.subdomain.deformation import deformation
from src.Graph import graph
from dolfin import *

num_steps = 100
dt = 1e-4
nx = ny = 50
c0 = 10
c1 = 100
name = f'c0_{c0}_c1_{c1}_nsteps_{num_steps}_dt_{dt}'
path = 'src/Neumann/Input/subdomain/results/'
deformation(dt, num_steps, name, path, c0, c1)
print('Graficando')
graph.graph(path, name, num_steps, lower_limit = -0.5, upper_limit = 0.5, step = 1, z_limit=True, save=True, plot=False)
