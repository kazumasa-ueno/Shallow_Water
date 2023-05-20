import numpy as np
import matplotlib.pyplot as plt
import glob

dilectname = '../out_stoc_n3/'

files = sorted(glob.glob(dilectname+'fmgz*.dat'))

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

# 仮の入力データ
rho = 1.0  # 質量密度 [kg/m^3]
g = 9.81   # 重力加速度 [m/s^2]

fig, (ax1, ax2) = plt.subplots(2, 1,figsize=(6,8))

# for i in range(levels):
for i in [2,3,4,5]:
	h = np.mean(z[:,:,i],axis=0)
	# print(h.shape)

	# 高さデータのフーリエ変換
	h_fft = np.fft.fft(h[:Nx//2**(levels-i-1)])

	# ポテンシャルエネルギー密度のスペクトルを計算
	energy_spectrum = 0.5 * rho * g**2 * np.abs(h_fft**2)

	# 波数kを計算
	N = len(h_fft)
	k = np.fft.fftfreq(N, d=1/N)
	xf = np.abs(k[:N//2])  # xfが正の値のみを含むようにする
# print(k[:Nx//2])

# スペクトルをプロット
	ax1.loglog(k[:N//2], energy_spectrum[:N//2], '-+')

x = np.logspace(0.8,1.3, 10)
y = 1e1*x**-2
ax1.loglog(x, y,'k')
# 目盛りを2の累乗に設定
x_ticks = [2**i for i in range(int(np.log2(xf[1])), int(np.log2(xf[-1]))+1) if 2**i >= xf[1] and 2**i <= xf[-1]]
ax2.set_xlim(1,256)
ax1.set_xticks(x_ticks, x_ticks)
ax1.text(10,1,"$k^{-2}$")
ax1.set_xlim(1,256)
ax1.set_ylim(1e-13,1e6)


files = sorted(glob.glob(dilectname+'z*.dat'))
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


# for i in range(levels):
for i in [2,3,4,5]:
	h = np.mean(z[:,:,i],axis=0)
	# print(h.shape)

	# 高さデータのフーリエ変換
	h_fft = np.fft.fft(h[:Nx//2**(levels-i-1)])

	# ポテンシャルエネルギー密度のスペクトルを計算
	energy_spectrum = 0.5 * rho * g**2 * np.abs(h_fft**2)

	# 波数kを計算
	N = len(h_fft)
	k = np.fft.fftfreq(N, d=1/N)
	xf = np.abs(k[:N//2])  # xfが正の値のみを含むようにする
# print(k[:Nx//2])

# スペクトルをプロット
	ax2.loglog(k[:N//2], energy_spectrum[:N//2], '-+', label="N = {}".format(Nx//2**(levels-i-1)))

ax2.loglog(x, y,'k')
ax2.text(10,1,"$k^{-2}$")

ax2.legend(fontsize=12,loc='lower left')
ax2.set_xticks(x_ticks, x_ticks)
ax2.set_xlabel('Wavenumber (k)',fontsize=12)
ax1.set_ylabel('output1',fontsize=12)
ax2.set_ylabel('output2',fontsize=12)
ax2.set_ylim(1e-13,1e6)

ax1.set_title('Potential Energy Spectrum [J/s]')

plt.tight_layout()
plt.savefig("../figure/spect.png")
