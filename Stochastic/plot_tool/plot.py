import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation
import glob

files = sorted(glob.glob('../output/z*.dat'))

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

# i = 5

fig,ax = plt.subplots()
# line, = ax.plot(x,z[0])
for i in range(levels):
    eta = z[:,:Nx//2**(levels-i-1),i]
    N = len(eta[0])
    x = np.linspace(0,10000,N)
    ax.plot(x, eta[0], label="N = {}".format(Nx//2**(levels-i-1)))
title = ax.set_title("num = {}".format(0),fontsize=14)
ax.set_ylim(-5,5)
ax.set_xlabel("x [km]")
ax.set_ylabel("$\eta$ [m]")
ax.legend()

waves = ax.get_lines()

def update(frame):
    for i, wave in enumerate(waves):
        wave.set_ydata(z[frame,:Nx//2**(levels-i-1),i])
    title.set_text("day = {:.2f}".format(frame*0.01))
    return waves

ani = FuncAnimation(fig, update, frames=range(Nt), interval=1, blit=True)

# ax.set_ylim(990,1010)
ani.save("../figure/line.gif",writer='imagemagick')