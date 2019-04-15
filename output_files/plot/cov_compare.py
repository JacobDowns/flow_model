import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
import seaborn as sns
import matplotlib

current_palette =  sns.color_palette("dark", 12)#sns.color_palette()
matplotlib.rcParams.update({'font.size': 22})
fig = plt.figure(figsize=(10, 6))
ax = plt.subplot(111)

ts = np.loadtxt('output_files/inversions/cov1/prior/sigma_ts.txt')
#x = np.loadtxt('output_files/inversions/cov1/prior/x.txt')
x1 = np.loadtxt('output_files/inversions/cov1/posterior/x_new.txt')
x2 = np.loadtxt('output_files/inversions/cov2/posterior/x_new.txt')
P1 = np.loadtxt('output_files/inversions/cov1/posterior/Pxx_new.txt')
P2 = np.loadtxt('output_files/inversions/cov2/posterior/Pxx_new.txt')
v1 = np.sqrt(np.diag(P1))
v2 = np.sqrt(np.diag(P2))


plt.fill_between(ts, x1 - 2.0*v1, x1 + 2.0*v1, color = current_palette[0], alpha = 0.3)
plt.fill_between(ts, x2 - 2.0*v2, x2 + 2.0*v2, color = current_palette[3], alpha = 0.45)

plt.plot(ts, x1, color = 'k', lw = 6., ms = 10, marker = 'o', alpha = 1.)
plt.plot(ts, x1, color = sns.color_palette()[0], lw = 4., label = "5th Order", ms = 7, marker = 'o', alpha = 1.)

plt.plot(ts, x2, color = 'k', lw = 6., ms = 10, marker = 'o', alpha = 0.8)
plt.plot(ts, x2, color = sns.color_palette()[3], lw = 4., label = "UT", ms = 7, marker = 'o', alpha = 0.8)

plt.grid(color='slategray', linestyle=':', linewidth=1)

plt.xlabel('Age (ka BP)')
plt.ylabel(r'$\Delta P$ (m.w.e. a$^{-1}$)')
plt.xlim([ts.min(), ts.max()])
plt.tight_layout()
ticks = ax.get_xticks()
ax.set_xticklabels([int(abs(tick / 1000.)) for tick in ticks])
plt.legend()
plt.savefig('cov_compare.png', dpi=400)

plt.show()
