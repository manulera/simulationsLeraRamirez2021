from simulation import Simulation, Parameters
import numpy as np
import os
from joblib import Parallel, delayed


p = Parameters()
p.v_slide = 0.35
p.v_growth = 1.6
p.v_shrink = 3.6
p.dt = 0.01
p.duration_n = 6.8
p.duration_r = 2.5
p.midzone_mu = 1.23
p.midzone_sigma = 0.25
p.ase1 = True
p.rearrange_mts = False


main_dir = "runs_ase1"

def runSim(label,i,par):

    sim = Simulation(par)
    output = sim.run()
    with open("%s/runs_%.1f/result_%02d.csv" % (main_dir,label, i), 'w') as out:
        out.write(output)

for i, resc in enumerate(np.arange(1, 120, 3)):
    print(i)
    os.makedirs('%s/runs_%.1f' % (main_dir,resc), exist_ok=True)
    p.total_rescue = resc
    Parallel(n_jobs=20,verbose=1)(delayed(runSim)(resc,i,p) for i in range(200))

