from simulation import Simulation, Parameters

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
p.total_rescue = 55.0
p.print_linkers = True
sim = Simulation(p)

sim.microtubules[2].lost = 1
sim.microtubules[6].lost = 1
sim.microtubules[7].lost = 1
print(sim.drawArrangement())
sim.performRearrangement(7)
print(sim.drawArrangement())