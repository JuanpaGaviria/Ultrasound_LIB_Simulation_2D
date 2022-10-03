import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import os

def graph(path, name, num_steps, lower_limit, upper_limit, step, z_limit, save, plot):
        plt.rcParams["figure.figsize"] = 15, 15

        H = np.loadtxt(path + f'{name}/H.txt')
        XX = np.loadtxt(path + f'{name}/XX.txt')
        YY = np.loadtxt(path + f'{name}/YY.txt')
        xvals = np.loadtxt(path + f'{name}/xvals.txt')
        yvals = np.loadtxt(path + f'{name}/yvals.txt')

        def makeplot(i):
                Z = griddata((xvals, yvals), H[i], (XX, YY),method='linear')
                plt.rcParams["figure.figsize"] = 12.8, 9.6
                fig = plt.figure(figsize=plt.figaspect(0.5))
                ax = fig.add_subplot(1, 2, 1, projection='3d')
                norm = plt.Normalize(Z.min(), Z.max())    
                colors = cm.jet(norm(Z))
                line = ax.plot_surface(XX, YY, Z, facecolors=colors, shade=False)
                line.set_facecolor((0,0,0,0))
                line = ax.set_xlabel('x')
                line = ax.set_ylabel('y')
                line = ax.set_zlabel('Amplitude')
                if z_limit:
                        line = ax.set_zlim(lower_limit, upper_limit)

                ax = fig.add_subplot(1, 2, 2, projection='3d')
                surf = ax.plot_surface(XX, YY, Z, cmap='jet')
                if z_limit:
                        surf = ax.set_zlim(lower_limit, upper_limit)
                surf = ax.set_xlabel('x')
                surf = ax.set_ylabel('y')
                surf = ax.set_zlabel('Amplitude')
                if save == True:
                        if not os.path.isdir(path + f'{name}/images'):
                                os.mkdir(path + f'{name}/images')
                        plt.savefig(path + f'{name}/images/iter' + str(i) +'.png')                       
                if plot == True:
                        plt.show()
                plt.close()
        for n in range(0, num_steps, step):
                makeplot(n)