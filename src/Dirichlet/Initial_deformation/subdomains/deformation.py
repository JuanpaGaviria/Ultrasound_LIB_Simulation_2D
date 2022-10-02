from fenics import *
import numpy as np
import os
from dolfin import *


def cosine_deformation_trapped(num_steps, dt, nx, ny, c2, name, path):
    mesh = RectangleMesh(Point(-0.005,-0.005),Point(0.005,0.005),nx,ny)
    V = FunctionSpace(mesh, 'P', 1)

    xyvals = mesh.coordinates()
    xvals = xyvals[:,0]
    yvals=xyvals[:,1]

    xx = np.linspace(-0.005,0.005)
    yy = np.linspace(-0.005,0.005)

    XX, YY = np.meshgrid(xx,yy)

    boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
    class Material_0(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1e-14
            return x[1] <= 0.005/2 + tol
        
    class Material_1(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1e-14
            return x[1] >= 0.005/2 + tol

    material0 = Material_0()
    material1 = Material_1()
    material0.mark(boundary_markers, 0)
    material1.mark(boundary_markers, 1)

    #  Initial values
    u_0 = Expression('cos(x[1])', degree=0)
    u_1 = Expression('cos(x[1])', degree=0)
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

def cosine_deformation(dt, num_steps, name, path, c0, c1):
    tol = 1E-14
    mesh = RectangleMesh(Point(-0.005,-0.005),Point(0.005,0.005),50,50)
    V = FunctionSpace(mesh, 'P', 1)
    xyvals = mesh.coordinates()
    xvals = xyvals[:,0]
    yvals=xyvals[:,1]

    xx = np.linspace(-0.005,0.005)
    yy = np.linspace(-0.005,0.005)

    XX, YY = np.meshgrid(xx,yy)
    class Omega0(SubDomain):
        def inside(slef, x, on_boundary):
            return True if x[0] <= 0 + tol else False
    class Omega1(SubDomain):
        def inside(self, x, on_boundary):
            return True if x[0] >= 0 + tol else False

    subdomains = MeshFunction('size_t', mesh, mesh.topology().dim(), 0)
    subdomain0 = Omega0().mark(subdomains, 0)
    subdomain1 = Omega1().mark(subdomains, 1)

    class K(UserExpression):
        def __init__(self, subdomains, k_0, k_1, **kwargs):
            super().__init__(**kwargs)
            self.subdomains = subdomains
            self.k_0 = k_0
            self.k_1 = k_1
        def eval_cell(self, values, x, cell):
            if self.subdomains[cell.index] == 0:
                values[0] = self.k_0
            else:
                values[0] = self.k_1

    c2 = K(subdomains=subdomains, k_0 = c0, k_1 = c1, degree=0)
    #  Boundary conditions
    u_D = Expression('0', degree=0)

    def boundary(x, on_boundary):
        return on_boundary
        
    bc = DirichletBC(V, u_D, boundary)
    #  Initial values

    u_0 = Expression('cos(x[0])', degree=0)
    u_1 = Expression('cos(x[0])', degree=0)
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
        solve(a == L, u_j1, bc)
        uvals = u_j1.compute_vertex_values(mesh)
        H.append(uvals)
        # Update previous solution
        u_j0.assign(u_j)
        u_j.assign(u_j1)
    H=np.array(H)
        
    # Plot solution
    if not os.path.isdir(path + f'{name}'):
        os.mkdir(path + f'{name}')
    np.savetxt(path + f'{name}/xvals.txt',xvals)
    np.savetxt(path + f'{name}/yvals.txt',yvals)
    np.savetxt(path + f'{name}/XX.txt',XX)
    np.savetxt(path + f'{name}/YY.txt',YY)
    np.savetxt(path + f'{name}/H.txt',H)

def sin_deformation(dt, num_steps, name, path, c0, c1):
    tol = 1E-14
    mesh = RectangleMesh(Point(-0.005,-0.005),Point(0.005,0.005),50,50)
    V = FunctionSpace(mesh, 'P', 1)
    xyvals = mesh.coordinates()
    xvals = xyvals[:,0]
    yvals=xyvals[:,1]

    xx = np.linspace(-0.005,0.005)
    yy = np.linspace(-0.005,0.005)

    XX, YY = np.meshgrid(xx,yy)
    class Omega0(SubDomain):
        def inside(slef, x, on_boundary):
            return True if x[0] <= 0 + tol else False
    class Omega1(SubDomain):
        def inside(self, x, on_boundary):
            return True if x[0] >= 0 + tol else False

    subdomains = MeshFunction('size_t', mesh, mesh.topology().dim(), 0)
    subdomain0 = Omega0().mark(subdomains, 0)
    subdomain1 = Omega1().mark(subdomains, 1)

    class K(UserExpression):
        def __init__(self, subdomains, k_0, k_1, **kwargs):
            super().__init__(**kwargs)
            self.subdomains = subdomains
            self.k_0 = k_0
            self.k_1 = k_1
        def eval_cell(self, values, x, cell):
            if self.subdomains[cell.index] == 0:
                values[0] = self.k_0
            else:
                values[0] = self.k_1

    c2 = K(subdomains=subdomains, k_0 = c0, k_1 = c1, degree=0)
    #  Boundary conditions
    u_D = Expression('0', degree=0)

    def boundary(x, on_boundary):
        return on_boundary
        
    bc = DirichletBC(V, u_D, boundary)
    #  Initial values

    u_0 = Expression('sin(x[0])', degree=0)
    u_1 = Expression('sin(x[0])', degree=0)
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
        solve(a == L, u_j1, bc)
        uvals = u_j1.compute_vertex_values(mesh)
        H.append(uvals)
        # Update previous solution
        u_j0.assign(u_j)
        u_j.assign(u_j1)
    H=np.array(H)
        
    # Plot solution
    if not os.path.isdir(path + f'{name}'):
        os.mkdir(path + f'{name}')
    np.savetxt(path + f'{name}/xvals.txt',xvals)
    np.savetxt(path + f'{name}/yvals.txt',yvals)
    np.savetxt(path + f'{name}/XX.txt',XX)
    np.savetxt(path + f'{name}/YY.txt',YY)
    np.savetxt(path + f'{name}/H.txt',H)


