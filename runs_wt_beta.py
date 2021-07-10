from simulation import Simulation, Parameters
import numpy as np
import os
from joblib import Parallel, delayed


p = Parameters()
p.v_slide = 0.35
p.v_growth = 1.6
p.v_shrink = 3.6
p.dt = 0.01
p.duration_n = 8.53
p.duration_r = 3.17
p.midzone_mu = 1.23
p.midzone_sigma = 0.25
p.ase1 = False
p.rearrange_mts = True


main_dir = "runs_wt_beta"

def runSim(label,a,i,par):

    sim = Simulation(par)
    output = sim.run()
    with open("%s/runs_%.1f_%d/result_%02d.csv" % (main_dir,label,a, i), 'w') as out:
        out.write(output)

for i, resc in enumerate(np.arange(1, 120, 3)):

    p.total_rescue = resc
    for a in [4, 8, 12]:
        print(i,a)
        p.alpha = a
        os.makedirs('%s/runs_%.1f_%d' % (main_dir, resc, a), exist_ok=True)
        Parallel(n_jobs=20,verbose=1)(delayed(runSim)(resc,a,i,p) for i in range(200))

