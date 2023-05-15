import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation

Nt = 150
Nx = 512

z = np.loadtxt("./output/z.txt")
z.reshape(Nt,Nx)
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