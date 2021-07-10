
import numpy as np
import matplotlib.pyplot as plt


append = ['wt','ase1_1','ase1_2','ase1_3','ase1_4']
files = [
    "runs_wt_beta/runs_55.0_4/result_08.csv",
    "runs_ase1/runs_34.0/result_05.csv",
    "runs_ase1/runs_34.0/result_03.csv",
    "runs_ase1/runs_34.0/result_10.csv",
    "runs_ase1/runs_34.0/result_12.csv",
]

for j in range(5):
    data = np.genfromtxt(files[j],delimiter=' ')

    mt_ids = data[:,0]
    sim_time = data[:,1]
    sim_pos = data[:,2]
    event_type = data[:,3]
    orientation = data[:,4]
    time_loss = np.append(0, sim_time[event_type==2])
    nb_mts0 = len(np.unique(mt_ids))
    nb_mts = np.arange(nb_mts0,nb_mts0-len(time_loss),-1)
    # nb_mts =
    plt.figure()

    colors = ['#dd7eae','#dd7eae','#0076b9']
    for i in np.unique(mt_ids):
        logi = mt_ids==i
        this_orientation = orientation[logi][0]
        plt.plot(sim_time[logi],-sim_pos[logi],colors[int(this_orientation)])

    t = np.linspace(0,20)
    plt.plot(t,2 + t*0.35,c='#5e5e5eff',lw=3)
    plt.plot(t,-2 - t*0.35,c='#5e5e5eff',lw=3)

    # Scalebars
    if j==0:
        plt.plot([0.5,6],[-7.5,-7.5],c='black',lw=3)
        plt.plot([0.5,0.5],[-7.5,-5.5],c='black',lw=3)
    if j==1:
        plt.plot([0.5,6],[-5.5,-5.5],c='black',lw=3)
        plt.plot([0.5,0.5],[-5.5,-3.5],c='black',lw=3)

    plt.savefig('../figures/figures/simulation_%s.svg' % append[j])
    plt.figure()
    plt.step(time_loss,nb_mts)

plt.show()


