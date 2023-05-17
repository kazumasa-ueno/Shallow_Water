import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation
import glob

files = sorted(glob.glob('./output/fmgz*.dat'))

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

z = z[:,:,-1]
x = np.linspace(0,4e6,Nx)

fig,ax = plt.subplots()
line, = ax.plot(x,z[0])
title = ax.set_title("num = {}".format(0))
ax.set_ylim(-0.05,0.05)

def update(frame):
    line.set_ydata(z[frame])
    title.set_text("day = {:2f}".format(frame*0.01))
    return line,

ani = FuncAnimation(fig, update, frames=range(Nt), interval=1, blit=True)

# ax.set_ylim(990,1010)
ani.save("line.gif",writer='imagemagick')