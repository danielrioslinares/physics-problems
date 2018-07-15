


import matplotlib.pyplot as plt
import numpy as np
import time
import eq

# Features
x = (np.linspace(-400,400,201) * 1e-7).tolist()
e = [8.854187817e-014 * 11.7 for _ in x]
NC = [1.08e10 for _ in x]
NV = [1.08e10 for _ in x]
EC = [0 for _ in x]
EV = [0 for _ in x]
ND = [2e16 if xx >= 0 else 0 for xx in x]
NA = [1e16 if xx < 0 else 0 for xx in x]
gD = [0 for _ in x]
gA = [0 for _ in x]
ED = [0 for _ in x]
EA = [0 for _ in x]
T = 300
MAXCNT = 300000
MINTOL = 1e-007

# Elapsed time
elapsed = time.time()
v = eq.equilibrium__Poisson_1D__BC_VonNeumann__direct_iterative__manual(x,e,NC,NV,EC,EV,ND,NA,gD,gA,ED,EA,T,MAXCNT,MINTOL)
print("DIRECT ITERATIVE",time.time() - elapsed)

# Plot it
plt.figure()
plt.plot(x,v)
plt.show()
