import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
import seaborn as sns
import matplotlib

current_palette =  sns.color_palette()
matplotlib.rcParams.update({'font.size': 16})

age = np.loadtxt('output_files/center_opt/age.txt')
adot = np.loadtxt('output_files/center_opt/adot.txt')
L = np.loadtxt('output_files/center_opt/L.txt')
L_interp = interp1d(age, L)
H = np.loadtxt('output_files/center_opt/H.txt')
H_ts = np.arange(age[0], age.max(), 20.) 
H_Ls = L_interp(H_ts) / 1e3
B = np.loadtxt('output_files/center_opt/B.txt')
xs = np.linspace(0., 1., len(H[0]))

H_ts /= 1e3

# Center 
L1_obs = np.array([406878, 396313, 321224, 292845, 288562, 279753]) / 1e3
# south
#L2_obs  = np.array([424777, 394942, 332430, 303738, 296659, 284686]) / 1e3
# Observation ages
obs_ages = np.array([-11554., -10284., -9024., -8064., -7234., 0.]) / 1e3
# Observation variances
obs_sigmas = np.array([0.4, 0.2, 0.2, 0.3, 0.3, 0.1])/2.

### Animate
###########################################################

# First
fig = plt.figure(figsize=(14,14))
ax1 = plt.subplot(2,1,1)
fig.set_tight_layout(True)
index = H_Ls.argmax()
plt.plot(xs*H_Ls[index], B[index][::-1], 'k', linewidth=6)
line1, = ax1.plot(xs*H_Ls[0], H[0], 'k', linewidth=7)
line2, = ax1.plot(xs*H_Ls[0], H[0],  color = current_palette[0], linewidth=3)
plt.xlim([0., H_Ls.max()])
plt.ylim([-250., 3400.])
plt.xlabel('Along Flow Length (km)')
plt.ylabel('Elevation (m)')

# Second
ax2 = plt.subplot(2,1,2)
ax2.plot(H_ts, H_Ls, 'k', lw = 5)
plt.ylabel('Glacier Length')
plt.xlabel('Age (ka BP)')
plt.xlim([H_ts.min(), H_ts.max()])
plt.ylim([250., 435.])
line3, = ax2.plot(H_ts[0], H_Ls[0] / 1e3,  'ro', ms = 10)

for i in range(len(obs_ages) - 1):
    ax2.plot([obs_ages[i] - 2.0*obs_sigmas[i], obs_ages[i] + 2.0*obs_sigmas[i]], [L1_obs[i], L1_obs[i]], 'k', lw = 5, ms = 6, alpha = 1.)
    #plt.plot([obs_ages[i] - 2.0*obs_sigmas[i], obs_ages[i] + 2.0*obs_sigmas[i]], [L1_obs[i], L1_obs[i]], color = current_palette[0], lw = 2, ms = 6, alpha = 1.)
    ax2.plot(obs_ages[i], L1_obs[i], 'ko', ms = 10)
    #plt.plot(obs_ages[i], L1_obs[i], color = current_palette[0], ms = 7, marker = 'o')


def update(i):
    label = 'Age: {0} (ka BP)'.format(round(-H_ts[i],1))
    print(label)
    # Update the line and the axes (with a new xlabel). Return a tuple of
    # "artists" that have to be redrawn for this frame.
    line1.set_xdata(xs*H_Ls[i])
    line1.set_ydata((B[i] + H[i])[::-1])
    line2.set_xdata(xs*H_Ls[i])
    line2.set_ydata((B[i] + H[i])[::-1])

    line3.set_xdata(H_ts[i])
    line3.set_ydata(H_Ls[i])
    #ax1.set_xlabel(label)
    return line1, ax1, ax2

if __name__ == '__main__':
    # FuncAnimation will call the 'update' function for each frame; here
    # animating over 10 frames, with an interval of 200ms between frames.
    anim = FuncAnimation(fig, update, frames=np.arange(0, H.shape[0]), interval=200)
    if len(sys.argv) > 1 and sys.argv[1] == 'save':
        anim.save('example.gif', dpi=120, writer='imagemagick')
    else:
        # plt.show() will just loop the animation forever.
        plt.show()
