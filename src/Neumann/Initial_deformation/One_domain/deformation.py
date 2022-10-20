from fenics import *
import numpy as np
import os


def sinusoidal_deformation(num_steps, dt, nx, ny, c2, name, path):

    mesh = RectangleMesh(Point(-0.005,-0.005),Point(0.005,0.005),nx,ny)
    V = FunctionSpace(mesh, 'P', 1)

    xyvals = mesh.coordinates()
    xvals = xyvals[:,0]
    yvals=xyvals[:,1]

    xx = np.linspace(-0.005,0.005)
    yy = np.linspace(-0.005,0.005)

    XX, YY = np.meshgrid(xx,yy)

    #  Initial values
    u_0 = Expression('sin(x[1])', degree=0)
    u_1 = Expression('sin(x[1])', degree=0)
    u_j0 = interpolate(u_0, V)  # past-initial j-1
    u_j = interpolate(u_1, V)  # present uj

    #  Variational problem
    c2 = 1

    u_j1 = TrialFunction(V)  # future uj+1
    v = TestFunction(V)
    a = ((u_j1*v)+(dt**2*c2*(dot(grad(u_j1),grad(v)))))*dx
    L = (2*u_j - u_j0)*v*dx
    u_j1 = Function(V)
    H=[] #list with solution history

    for i in range(num_steps):
        solve(a == L, u_j1)
        uvals = u_j1.compute_vertex_values(mesh)
        H.append(uvals)
        # Update previous solution
        u_j0.assign(u_j)
        u_j.assign(u_j1)
    H=np.array(H)
        
    # Plot solution
    if not os.path.isdir(path + f'{name}'):
        os.mkdir(path + f'{name}')

    #save data
    np.savetxt(path + f'{name}/xvals.txt',xvals)
    np.savetxt(path + f'{name}/yvals.txt',yvals)
    np.savetxt(path + f'{name}/XX.txt',XX)
    np.savetxt(path + f'{name}/YY.txt',YY)
    np.savetxt(path + f'{name}/H.txt',H)

def cosine_deformation(num_steps, dt, nx, ny, c2, name, path):

    mesh = RectangleMesh(Point(-0.005,-0.005),Point(0.005,0.005),nx,ny)
    V = FunctionSpace(mesh, 'P', 1)

    xyvals = mesh.coordinates()
    xvals = xyvals[:,0]
    yvals=xyvals[:,1]

    xx = np.linspace(-0.005,0.005)
    yy = np.linspace(-0.005,0.005)

    XX, YY = np.meshgrid(xx,yy)

    #  Initial values
    u_0 = Expression('cos(x[1])', degree=0)
    u_1 = Expression('cos(x[1])', degree=0)
    u_j0 = interpolate(u_0, V)  # past-initial j-1
    u_j = interpolate(u_1, V)  # present uj

    #  Variational problem

    u_j1 = TrialFunction(V)  # future uj+1
    v = TestFunction(V)
    a = ((u_j1*v)+(dt**2*c2*(dot(grad(u_j1),grad(v)))))*dx
    L = (2*u_j - u_j0)*v*dx
    u_j1 = Function(V)
    H=[] #list with solution history

    for i in range(num_steps):
        solve(a == L, u_j1)
        uvals = u_j1.compute_vertex_values(mesh)
        # plot(u_j1)
        # plt.show()
        H.append(uvals)
        # Update previous solution
        u_j0.assign(u_j)
        u_j.assign(u_j1)
    H=np.array(H)
    
    # Plot solution
    if not os.path.isdir(path + f'{name}'):
        os.mkdir(path + f'{name}')

    #save data
    np.savetxt(path + f'{name}/xvals.txt',xvals)
    np.savetxt(path + f'{name}/yvals.txt',yvals)
    np.savetxt(path + f'{name}/XX.txt',XX)
    np.savetxt(path + f'{name}/YY.txt',YY)
    np.savetxt(path + f'{name}/H.txt',H)