import os
from Dirichlet.Initial_deformation.One_domain import deformation
from src.Graph import graph

def loop_sin_deformation(num_steps, dt, nx, ny, c2):

    for _dt in dt:
        dt = _dt
        name = f'nsteps_{num_steps}_dt{dt}_nx{nx}_ny_{ny}'
        path = 'src/Dirichlet/Initial_deformation/One_domain/results/'
        print(dt)
        deformation.sinusoidal_defformation(num_steps, dt, nx, ny, c2, name, path)
        print("saving images")
        graph.graph(path, name, num_steps, step = 1, save=True, plot=False)
