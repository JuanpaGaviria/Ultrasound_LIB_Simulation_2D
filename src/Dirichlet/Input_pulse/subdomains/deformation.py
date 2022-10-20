from fenics import *
import numpy as np
import os


def deformation(dt, num_steps, name, path, c0, c1):
    tol = 1E-14
    input_values = [1, 0.75, 0.5, 0.25, 0, -0.25, -0.5, -0.75, -1]
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
    input = Expression('0', degree=0)
    def input_boundary(x, on_boundary):
        s1 = on_boundary and near(x[0], 0.005, tol)
        return s1

    input_bc = DirichletBC(V, input, input_boundary)

    def boundary(x, on_boundary):
        s2 = on_boundary and not near(x[0], 0.005, tol)
        return s2

    u_D = Expression('0', degree = 0)
    bc = DirichletBC(V, u_D, boundary)

    bcs = [input_bc, bc]
    #  Initial values

    u_0 = Expression('0', degree=0)
    u_1 = Expression('0', degree=0)
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
        solve(a == L, u_j1, bcs)
        uvals = u_j1.compute_vertex_values(mesh)
        H.append(uvals)
        # Update previous solution
        u_j0.assign(u_j)
        u_j.assign(u_j1)
        if i < len(input_values):
            input = input_values[i]
            print(input)
        else:
            input = 0
        input_bc = DirichletBC(V, input, input_boundary)
        bcs = [input_bc, bc]
    H=np.array(H)
        
    # Plot solution
    if not os.path.isdir(path + f'{name}'):
        os.mkdir(path + f'{name}')
    np.savetxt(path + f'{name}/xvals.txt',xvals)
    np.savetxt(path + f'{name}/yvals.txt',yvals)
    np.savetxt(path + f'{name}/XX.txt',XX)
    np.savetxt(path + f'{name}/YY.txt',YY)
    np.savetxt(path + f'{name}/H.txt',H)