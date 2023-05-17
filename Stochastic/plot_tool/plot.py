import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation
import glob

files = sorted(glob.glob('../out_stoc_n3/fmgz*.dat'))

data = np.loadtxt(files[0])
Nt = len(files)
Nx = len(data[0])
levels = len(data)
print(Nt,Nx,levels)

z = np.zeros((Nt, Nx, levels))
for i, file in enumerate(files):
	data = np.loadtxt(file)
	for level in range(levels):
		z[i,:,level] = data[level]

i = 5
z = z[:,:Nx//2**(levels-i-1),i]
N = len(z[0])
x = np.linspace(0,1e7,N)

fig,ax = plt.subplots()
line, = ax.plot(x,z[0])
title = ax.set_title("num = {}".format(0))
ax.set_ylim(-5,5)

def update(frame):
    line.set_ydata(z[frame])
    title.set_text("day = {:2f}".format(frame*0.01))
    return line,

ani = FuncAnimation(fig, update, frames=range(Nt), interval=1, blit=True)

# ax.set_ylim(990,1010)
ani.save("../figure/line.gif",writer='imagemagick')