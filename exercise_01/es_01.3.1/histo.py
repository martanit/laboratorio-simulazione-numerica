import numpy as np
import matplotlib.pyplot as plt
import scipy.stats


S1st = np.loadtxt("output_es01.2_st1.dat", usecols=1)
S2st = np.loadtxt("output_es01.2_st2.dat", usecols=1)
S10st = np.loadtxt("output_es01.2_st10.dat", usecols=1)
S100st = np.loadtxt("output_es01.2_st100.dat", usecols=1)

S1exp = np.loadtxt("output_es01.2_exp1.dat", usecols=1)
S2exp = np.loadtxt("output_es01.2_exp2.dat", usecols=1)
S10exp = np.loadtxt("output_es01.2_exp10.dat", usecols=1)
S100exp = np.loadtxt("output_es01.2_exp100.dat", usecols=1)

S1cl = np.loadtxt("output_es01.2_cl1.dat", usecols=1)
S2cl = np.loadtxt("output_es01.2_cl2.dat", usecols=1)
S10cl = np.loadtxt("output_es01.2_cl10.dat", usecols=1)
S100cl = np.loadtxt("output_es01.2_cl100.dat", usecols=1)


n_bins=100
fig1 = plt.figure()
ax1 = fig1.add_subplot(1,1,1)
ax1.set_title('Normal distribution')
(mu, sigma) = scipy.stats.norm.fit(S100st)

y=scipy.stats.norm.pdf(n_bins, mu, sigma)

plt.plot(n_bins, y, 'b--')

ax1.hist([S1st, S2st, S10st, S100st], n_bins, density=True, label=['S1', 'S2', 'S10', 'S100'], range=[0,2], histtype='step')
ax1.legend ()
ax1.set_xlabel('Random number')
ax1.set_ylabel('Probability')

fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
ax2.set_title('Exponential distribution')
ax2.hist([S1exp, S2exp, S10exp, S100exp], n_bins, density=True, label=['S1', 'S2', 'S10', 'S100'], range=[0,4], histtype='step')
ax2.legend ()
ax2.set_xlabel('Random number')
ax2.set_ylabel('Probability')

fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
ax3.set_title('Cauchy Lorentz distribution')
ax3.hist([S1cl, S2cl, S10cl, S100cl], n_bins, density=True, label=['S1', 'S2', 'S10', 'S100'], range=[-20,20], histtype='step' )
ax3.legend ()
ax3.set_xlabel('Random number')
ax3.set_ylabel('Probability')

plt.grid(True)

plt.show()

