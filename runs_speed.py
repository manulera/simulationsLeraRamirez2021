from simulation import Simulation, Parameters
import numpy as np
import os
from joblib import Parallel, delayed

p = Parameters()
p.v_slide = 0.35
# p.v_growth = 1.6
p.v_shrink = 3.6
p.dt = 0.01
p.duration_n = 8.53
p.duration_r = 3.17
p.midzone_mu = 1.23
p.midzone_sigma = 0.25
p.ase1 = False
p.rearrange_mts = True
p.total_rescue = 55.0
p.redistribute_rescue = True

main_dir = "runs_wt_speed"

def runSim(label,i,par):

    sim = Simulation(par)
    output = sim.run()
    with open("%s/runs_%.4f/result_%02d.csv" % (main_dir,label, i), 'w') as out:
        out.write(output)


for j, v_growth in enumerate(np.arange(0.35, 1.5, 0.125/4.)):
    print(j)
    os.makedirs('%s/runs_%.4f' % (main_dir,v_growth), exist_ok=True)
    p.v_growth = v_growth

    Parallel(n_jobs=20,verbose=1)(delayed(runSim)(v_growth,i,p) for i in range(500))