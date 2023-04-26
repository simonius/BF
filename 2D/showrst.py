# Python helper script that plots the resulting HDF5 file


import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as patches
import h5py
import solvers
import numpy as np

f = h5py.File("ffsout240.h5", "r")
u = f["u"][:, :, :, :]


C, L, K, N = u.shape

kstep, lstep = int(solvers.ffsdomain.kstep), int(solvers.ffsdomain.lstep)+int(solvers.consts.lmrorder)
kwidth, lwidth = K-kstep, L - lstep
nend = 36

p = np.zeros([L, K])

for k in range(K):
        for l in range(L):
                p[l, k] = solvers.fluxes.p(u[:, l, k, nend])

for k in range(kwidth):
        for l in range(lstep):
                u[:, l, kstep+k, nend] = np.nan
                p[l, kstep+k] = np.nan
          
minv = np.nanmin(u[0, 3:L-3, 1:K-3, nend])
maxv = np.nanmax(u[0, 3:L-3, 1:K-3, nend])
levs = np.linspace(minv, maxv, num=30)
plt.contour(u[0, 3:L-3, 1:K-3, nend], levels=levs, colors='k', linewidths=0.3)
plt.contourf(u[0, 3:L-3, 1:K-3, nend], levels=levs, extend='both')
ax = plt.gca()
ax.set_aspect('equal', 'box')
plt.savefig("ffstestrho240.pdf", bbox_inches='tight')
plt.clf()

plt.contour(p[3:L-3, 1:K-3], levels=30, colors='k', linewidths=0.3)
plt.contourf(p[3:L-3, 1:K-3], levels=30)
ax = plt.gca()
ax.set_aspect('equal', 'box')
plt.savefig("ffstestp240.pdf", bbox_inches='tight')

