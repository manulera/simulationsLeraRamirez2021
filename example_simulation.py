from simulation import Simulation, Parameters

# Set up the parameters ----------------
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

# If you want to print extra info
p.print_linkers = True

# Run the simulation ----------------
sim = Simulation(p)

# Print the results
results = sim.run()

with open('results.csv','w') as output_file:
    output_file.write(results)