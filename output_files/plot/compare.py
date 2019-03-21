import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
import seaborn as sns
import matplotlib

current_palette =  sns.color_palette()
matplotlib.rcParams.update({'font.size': 22})

x1 = np.loadtxt('output_files/inversions/covariance_comparison/posterior/x_new.txt')
x2 = np.loadtxt('output_files/inversions/covariance_comparison1/posterior/x_new.txt')
P1 = np.loadtxt('output_files/inversions/covariance_comparison/posterior/Pxx_new.txt')
P2 = np.loadtxt('output_files/inversions/covariance_comparison1/posterior/Pxx_new.txt')
v1 = np.sqrt(np.diag(P1))
v2 = np.sqrt(np.diag(P2))

plt.plot(x1 - 2.0*v1, 'b')
plt.plot(x1, 'b')
plt.plot(x1 + 2.0*v1, 'b')

plt.plot(x2 - 2.0*v2, 'r')
plt.plot(x2, 'r')
plt.plot(x2 + 2.0*v2, 'r')



plt.show()
quit()

# Center 
L1_obs = np.array([406878, 396313, 321224, 292845, 288562, 279753]) / 1e3
# south
L1_obs  = np.array([424777, 394942, 332430, 303738, 296659, 284686]) / 1e3
# Observation ages
obs_ages = np.array([-11554., -10284., -9024., -8064., -7234., 0.]) / 1e3
# Observation variances
obs_sigmas = np.array([0.4, 0.2, 0.2, 0.3, 0.3, 0.1])/2.

### Animate
###########################################################

# First
fig = plt.figure(figsize=(14.5,18))
ax1 = plt.subplot(3,1,1)
fig.set_tight_layout(True)
index = H_Ls.argmax()
plt.plot(xs*H_Ls[index], B[index][::-1], 'k', linewidth=6)
line1, = ax1.plot(xs*H_Ls[0], H[0], 'k', linewidth=7)
line2, = ax1.plot(xs*H_Ls[0], H[0],  color = current_palette[0], linewidth=3)
plt.xlim([0., H_Ls.max()])
plt.ylim([-250., 3600.])
plt.xlabel('Along Flow Length (km)')
plt.ylabel('Elevation (m)')

# Second
ax2 = plt.subplot(3,1,2)

ax2.plot(H_ts, H_Ls, color = 'k', lw = 6)
ax2.plot(H_ts, H_Ls, color = current_palette[0], lw = 3)
for i in range(len(obs_ages) - 1):
    ax2.plot([obs_ages[i] - 2.0*obs_sigmas[i], obs_ages[i] + 2.0*obs_sigmas[i]], [L1_obs[i], L1_obs[i]], color = 'k', lw = 5, ms = 6, alpha = 1.)
    #ax2.plot([obs_ages[i] - 2.0*obs_sigmas[i], obs_ages[i] + 2.0*obs_sigmas[i]], [L1_obs[i], L1_obs[i]], color = current_palette[0], lw = 2, ms = 6, alpha = 1.)
    ax2.plot(obs_ages[i], L1_obs[i], 'ko', ms = 12)
    #ax2.plot(obs_ages[i], L1_obs[i], color = current_palette[0], ms = 7, marker = 'o')

plt.ylabel('Glacier Length (km)')
plt.xlabel('Age (ka BP)')
plt.xlim([H_ts.min(), H_ts.max()])
plt.ylim([275., 435.])
line3_a, = ax2.plot(H_ts[0], H_Ls[0] / 1e3,  'ko', ms = 18)
line3, = ax2.plot(H_ts[0], H_Ls[0] / 1e3, marker = 'o', color = current_palette[3], ms = 13)



# Third
ax3 = plt.subplot(3,1,3)
dt_years = np.loadtxt('paleo_data/buizert_ages.txt') / 1e3
dt_vals = np.loadtxt('paleo_data/buizert_dts.txt').T

dt_avg = interp1d(dt_years, dt_vals.mean(axis = 0))(H_ts)
ax3.plot(H_ts, dt_avg, linewidth = 4, color = 'k')
line4_a, = ax3.plot(H_ts[0],  dt_avg[0],  'ko', ms = 18)
line4, = ax3.plot(H_ts[0],  dt_avg[0],  marker = 'o', color = current_palette[3], ms = 13)


plt.ylabel(r'$\Delta T$ ($^{\circ}$ C)')
plt.xlabel('Age (ka BP)')
plt.xlim([H_ts.min(), H_ts.max()])

def update(i):
    label = 'Age: {0} (ka BP)'.format(round(-H_ts[i],1))
    print(label)
    # Update the line and the axes (with a new xlabel). Return a tuple of
    # "artists" that have to be redrawn for this frame.
    line1.set_xdata(xs*H_Ls[i])
    line1.set_ydata((B[i] + H[i])[::-1])
    line2.set_xdata(xs*H_Ls[i])
    line2.set_ydata((B[i] + H[i])[::-1])

    line3_a.set_xdata(H_ts[i])
    line3_a.set_ydata(H_Ls[i])
    line3.set_xdata(H_ts[i])
    line3.set_ydata(H_Ls[i])

    line4_a.set_xdata(H_ts[i])
    line4_a.set_ydata(dt_avg[i])
    line4.set_xdata(H_ts[i])
    line4.set_ydata(dt_avg[i])

    #ax1.set_xlabel(label)
    return line1, ax1, ax2

if __name__ == '__main__':
    # FuncAnimation will call the 'update' function for each frame; here
    # animating over 10 frames, with an interval of 200ms between frames.
    anim = FuncAnimation(fig, update, frames=np.arange(0, H.shape[0]), interval=200)
    if len(sys.argv) > 1 and sys.argv[1] == 'save':
        anim.save('example.gif', dpi=100, writer='imagemagick')
    else:
        # plt.show() will just loop the animation forever.
        plt.show()
