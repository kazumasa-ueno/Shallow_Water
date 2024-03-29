import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation

z = np.loadtxt("./output/z.txt")
z.reshape(300,1025)
x = np.linspace(0,2e6,1025)

fig,ax = plt.subplots()
line, = ax.plot(x,z[0])
title = ax.set_title("num = {}".format(0))
ax.set_ylim(-1,1)

def update(frame):
    line.set_ydata(z[frame])
    title.set_text("day = {:2f}".format(frame*9000/86400))
    return line,

ani = FuncAnimation(fig, update, frames=range(300), interval=1, blit=True)

# ax.set_ylim(990,1010)
ani.save("line.gif",writer='imagemagick')