import numpy as np
import matplotlib.pyplot as plt

Nt = 150
Nx = 512

# 仮の入力データ
rho = 1.0  # 質量密度 [kg/m^3]
g = 9.81   # 重力加速度 [m/s^2]
h = np.loadtxt("./output/z.txt")  # 高さ [m] (長さ1024の配列)
h = np.mean(h,axis=0)
print(h.shape)

# 高さデータのフーリエ変換
h_fft = np.fft.fft(h)

# ポテンシャルエネルギー密度のスペクトルを計算
energy_spectrum = 0.5 * rho * g**2 * np.abs(h_fft**2)

# 波数kを計算
N = len(h_fft)
k = np.fft.fftfreq(N, d=1/N)
xf = np.abs(k[:Nx//2])  # xfが正の値のみを含むようにする
# print(k[:Nx//2])

# スペクトルをプロット
plt.loglog(k[:Nx//2], energy_spectrum[:Nx//2])
# plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter())
# plt.gca().xaxis.set_minor_formatter(plt.ScalarFormatter())
# plt.ticklabel_format(axis='x', style='plain')
plt.xlabel('Wavenumber (k)')
plt.ylabel('Potential Energy Spectrum')

# 目盛りを2の累乗に設定
x_ticks = [2**i for i in range(int(np.log2(xf[1])), int(np.log2(xf[-1]))+1) if 2**i >= xf[1] and 2**i <= xf[-1]]
plt.xticks(x_ticks, x_ticks)

plt.tight_layout()
plt.savefig("spect.png")
