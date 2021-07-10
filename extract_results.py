import os
import numpy as np
from scipy.interpolate import interp1d

v_sliding = 0.35

def padcat(lis):
    """
    Concatenate lists into a single numpy array, filling the holes with nans
    :param lis:
    :return:
    """
    max_len = max(list(map(len,lis)))

    out = np.empty([len(lis),max_len])
    out[:] = np.nan
    for i,l in enumerate(lis):
        out[i,:len(l)] = l

    return out


for main_dir in ['runs_ase1_random']:

    # For each condition
    folders = os.listdir(main_dir)
    folders.sort()
    condition_dirs = [os.path.join(main_dir,i) for i in folders if os.path.isdir(os.path.join(main_dir,i)) and i[0] != '.']

    for condition_dir in condition_dirs:
        print(condition_dir)
        rescue_positions = list()
        catastrophe_positions = list()
        nb_mts_all = list()
        time_loss_all = list()
        time_total_length = np.linspace(0, 20)
        total_length_all = list()
        half_spindle_length = 2 + time_total_length * v_sliding
        simulations = [os.path.join(condition_dir,i) for i in os.listdir(condition_dir) if i.endswith('.csv')]

        for sim in simulations:
            data = np.genfromtxt(sim,delimiter=' ')

            mt_ids = data[:, 0]
            sim_time = data[:, 1]
            sim_pos = data[:, 2]
            event_type = data[:, 3]
            orientation = data[:, 4]
            sim_pos = sim_pos * orientation

            time_loss = np.append(0, sim_time[event_type == 2])
            nb_mts0 = len(np.unique(mt_ids))
            nb_mts = np.arange(nb_mts0, nb_mts0 - len(time_loss), -1)
            nb_mts_all.append(nb_mts)
            time_loss_all.append(time_loss)
            rescue_positions += list(sim_pos[event_type == 1])
            catastrophe_positions += list(sim_pos[event_type == 0])

            total_length = np.zeros_like(time_total_length)
            for i in np.unique(mt_ids):
                y = sim_pos[mt_ids == i]
                t = sim_time[mt_ids == i]
                interpolator = interp1d(t, y, bounds_error=False, fill_value=np.nan)
                # nans is when the mt is gone
                mt_length = interpolator(time_total_length) + half_spindle_length
                idxs = np.logical_not(np.isnan(mt_length))
                total_length[idxs] += mt_length[idxs]

            total_length_all.append(total_length)

        summary_dir = os.path.join(condition_dir, 'summary')
        os.makedirs(summary_dir, exist_ok=True)
        np.savetxt(os.path.join(summary_dir,'rescue_positions.csv'),rescue_positions)
        np.savetxt(os.path.join(summary_dir, 'catastrophe_positions.csv'), catastrophe_positions)
        np.savetxt(os.path.join(summary_dir, 'number_microtubules.csv'), padcat(nb_mts_all))
        np.savetxt(os.path.join(summary_dir, 'microtubule_loss_time.csv'), padcat(time_loss_all))
        np.savetxt(os.path.join(summary_dir, 'spindle_length.csv'), half_spindle_length*2)
        np.savetxt(os.path.join(summary_dir, 'polymer_length.csv'), padcat(total_length_all))





