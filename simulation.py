from random import random
from math import log, exp
import numpy as np
from scipy.stats import beta

class Parameters:
    """
    A class to store the parameters of the simulation.
    All the magnitudes are in minutes and micrometers
    """

    def __init__(self):

        # Velocity of microtubule sliding (half of spindle elongation velocity)
        self.v_slide = 0.0

        # Velocity of growth of microtubules
        self.v_growth = 0.0

        # Velocity of shrinkage of microtubules
        self.v_shrink = 0.0

        # The dt for the simulation
        self.dt = 0.0

        # We have fit the wild-type data of position of midzone edge to a normal distribution, these are the mu and the
        # sigma of the fit:
        self.midzone_mu = 0.0
        self.midzone_sigma = 0.0

        # Parameters of the fit of duration of growth events:
        self.duration_r = 0.0
        self.duration_n = 0.0

        # The total rescue that will be distributed along links (in wt) or along the length of mts (ase1)
        self.total_rescue = 0.0

        # Boolean settings
        # Whether we rearrange microtubules when one is lost
        self.rearrange_mts = None

        # You can set this to true to print a cartoon of the spindle whenever a microtubule is lost and when there is a
        # rearrangement
        self.print_linkers = False

        # Whether we simulate the ase1 condition or the wt
        self.ase1 = None

        # The alpha/beta parameters of the beta distribution, if we use it to define the rescue probability
        # in the midzone
        self.alpha = 0.0
        self.beta = 2.0

class Simulation:

    def __init__(self, par):
        """
        :param par:
        :type par:Parameters
        """
        # An instance of the Parameters class
        self.par = par

        # A string containing the content of the output file (see Simulation::run)
        self.mt_coordinates = ""

        # The simulation time
        self.t = 0

        # We sample the position of the midzone edge by random sample of the normal distribution that we fitted to the
        # midzone edge data.
        self.midzone_edge = np.random.normal(self.par.midzone_mu, self.par.midzone_sigma, 1)[0]

        # The array of microtubules (see microtubule class)
        #
        #         0 -- 1 -- 2
        #         |    |    |
        #         3 -- 4 -- 5
        #         |    |    |
        #         6 -- 7 -- 8
        #
        self.microtubules = [
            Microtubule(0, 1, 1, self),  # top left
            Microtubule(1, 1, -1, self),  # top center
            Microtubule(2, 1, 1, self),  # top right
            Microtubule(3, 1, -1, self),  # center left
            Microtubule(4, 1, 1, self),  # center center
            Microtubule(5, 1, -1, self),  # center right
            Microtubule(6, 1, 1, self),  # bottom left
            Microtubule(7, 1, -1, self),  # bottom center
            Microtubule(8, 1, 1, self),  # bottom right
        ]

        # NeighbourList - Each position in the list corresponds to a neighbour_index, and contains the neighbour_index
        # of neighbouring microtubules. For example, position 0, contains the indexes of neighbours (1,3) see cartoon

        self.neighbour_list = [
            [1, 3],  # top left
            [0, 2, 4],  # top center
            [1, 5],  # top right
            [0, 4, 6],  # center left
            [1, 3, 5, 7],  # center center
            [2, 4, 8],  # center right
            [3, 7],  # bottom left
            [4, 6, 8],  # bottom center
            [5, 7],  # bottom right
        ]

        # The initial position of the spindle pole with respect to the spindle center, we start with a spindle of 4 um
        self.spindle_edge = 2.

        # Array that contains the probability of rescue in dt for each microtubule (depending on their id)
        # Updates when a microtubule is lost
        self.prob_rescue = []

        # Total rate for each microtubule (depending on their id, useful for the beta distrib)
        # Updates when a microtubule is lost
        self.rate_rescue = []

        self.updateProbRescue()

    def timeToNextCatastrophe(self):
        """
        Get the time of next catastrophe by random sample of the distribution (1-exp(-r*t))^n
        :return:
        """
        prob = random()
        return -(log(1 - prob ** (1 / self.par.duration_n)) / self.par.duration_r)

    def updateProbRescue(self):
        """
        Update the probability of rescues per number of neighbours. This is updated when a microtubule is lost.
        :return:
        """
        if self.par.ase1:
            return

        linkers = [mt.countNeighbours() if not mt.lost else 0 for mt in self.microtubules]

        # Each linker is counted twice since we iterate through all microtubules
        total_linkers = sum(linkers) / 2.

        # We stop the function just to prevent a warning when dividing by zero (the only case when this function is
        # called when total_linkers is zero is when the last microtubule depolymerises)
        if total_linkers==0:
            return

        # The rescue density is homogeneously distributed
        rate_per_linker = self.par.total_rescue / self.midzone_edge / 2. / total_linkers

        if self.par.alpha == 0.0:
            # This is a list with the probability of being rescued, the index corresponds to the number
            # of neighbours
            rescue_per_neighbour = [1 - exp(-rate_per_linker * self.par.dt * neigh) for neigh in range(6)]

            # Here the index corresponds to the id of the microtubule
            self.prob_rescue = [rescue_per_neighbour[nb_neigh] for nb_neigh in linkers]
        else:
            rate_per_neighbour = [rate_per_linker * neigh for neigh in range(6)]
            self.rate_rescue = [rate_per_neighbour[nb_neigh] for nb_neigh in linkers]

        if self.par.print_linkers:
            print("Number of linkers for each mt_id:", linkers)
            print(self.drawArrangement())

    def updateProbRescueAse1(self):
        """
        Update the probability of rescues for ase1 spindles (the total amount of rescue distributed along mts)
        :return:
        """
        # Sum the length of mts:
        total_length = 0
        for mt in self.microtubules:
            if not mt.lost:
                total_length += mt.pos * mt.orientation + self.spindle_edge

        # Divided by two because now it does not distribute on the interface
        self.prob_rescue = [1. - exp(-self.par.total_rescue / total_length / 2. * self.par.dt)]

    def drawArrangement(self):
        """
        Draw a cartoon of the arrangement of microtubules, where the id microtubule is shown if the microtubule has not
        been lost yet, and otherwise an 'x' is shown
        :return:
        """
        arrangement = "n0 -- n1 -- n2\n|    |    |\nn3 -- n4 -- n5\n|    |    |\nn6 -- n7 -- n8\n"
        for mt in self.microtubules:
            neigh = mt.neighbour_id
            replacement = 'x' if mt.lost else str(mt.id)
            arrangement = arrangement.replace("n" + str(neigh), replacement)
        return arrangement

    def swapMicrotubules(self,mt_1_id,mt_2_id):

        mt_1_old_neighbour_id = self.microtubules[mt_1_id].neighbour_id
        mt_2_old_neighbour_id = self.microtubules[mt_2_id].neighbour_id

        self.microtubules[mt_2_id].neighbour_id = mt_1_old_neighbour_id
        self.microtubules[mt_1_id].neighbour_id = mt_2_old_neighbour_id

        if self.par.print_linkers:
            print("Swapped", mt_1_id, mt_2_id)


    def performRearrangement(self,lost_mt_id):
        """
        The rearrangement can happen in two ways:
            1) A mt is lost and there is an existing one in the same orientation which has
             less neighbours, we exchange them -> see verify_arrangement_case1.py
                Example: 7 is lost
                 x -- 1 -- 2         x -- 1 -- 2
                 |    |    |  3<->7  |    |    |
                 3 -- 4 -- 5    ->   x -- 4 -- 5
                 |    |    |         |    |    |
                 6 -- 7 -- 8         6 -- 3 -- 8
            2) A neighbour of the lost microtubule would have more neighbours if it was in another position
                Example: 5 is lost -> see verify_arrangement_case2.py
                 0 -- 1 -- x         0 -- 1 -- x
                 |    |    |  8<->6  |    |    |
                 3 -- 4 -- 5    ->   3 -- 4 -- x
                 |    |    |         |    |    |
                 x -- 7 -- 8         8 -- 7 -- 6
        :param lost_mt_id:
        :return:
        """

        # We construct an array with (orientation, id, nb_neighbours, lost, neighbour_id) for each mt
        the_array = np.array( [[mt.orientation, mt.id, mt.countNeighbours(), mt.lost, mt.neighbour_id] for mt in self.microtubules])

        # Case 1: other existing microtubules in the same orientation that have less neighbours
        case1 = np.logical_and.reduce((
            the_array[:, 0] == the_array[lost_mt_id, 0],  # oriented as the lost mt
            the_array[:, 1] != lost_mt_id,                # different from the lost mt
            the_array[:, 2] < the_array[lost_mt_id, 2],   # with less neighbours
            the_array[:, 3] == 0                          # not lost
        ))
        if np.any(case1):
            # Pick a random one
            picked_mt_id = np.random.choice(the_array[case1, 1])
            self.swapMicrotubules(lost_mt_id,picked_mt_id)
            if self.par.print_linkers:
                print("Swap case 1")
            return

        # Case 2: A neighbour of the lost microtubule would have more neighbours if it was in another empty position
        # We pick a random neighbour among those with the least neighbours
        neighbours_of_lost_mt = np.logical_and(
            np.isin(the_array[:,4],self.neighbour_list[lost_mt_id]),  # Neighbour of the lost microtubule
            the_array[:,3] == 0                                       # Not lost
        )
        if not np.any(neighbours_of_lost_mt):
            return
        minimum_nb_neighbours = np.min(the_array[neighbours_of_lost_mt,2])
        microtubules_to_pick_from = np.logical_and(
            neighbours_of_lost_mt,
            the_array[:, 2] == minimum_nb_neighbours
        )

        case2 = np.logical_and.reduce((
            the_array[:, 0] != the_array[lost_mt_id, 0],  # oriented opposite to the lost mt
            the_array[:, 2] > minimum_nb_neighbours,         # with more neighbours than the minimum
            the_array[:, 3] == 1                          # lost
        ))

        if np.any(case2):
            neighbour_to_swap = np.random.choice(the_array[microtubules_to_pick_from, 1])
            empty_slot_to_swap = np.random.choice(the_array[case2,1])
            self.swapMicrotubules(neighbour_to_swap,empty_slot_to_swap)
            if self.par.print_linkers:
                print("Swap case 2")

    def addLostMicrotubule(self, lost_mt_id):
        """
        Function before a microtubule is lost, to rearrange microtubules if needed, and to recalculate probability of
        rescue
        :param lost_mt_id:
        :return:
        """
        if self.par.print_linkers:
            print(">> mt %d lost" % lost_mt_id)

        self.microtubules[lost_mt_id].lost = True

        if self.par.rearrange_mts and not self.par.ase1:
            self.performRearrangement(lost_mt_id)

        # If all mts are lost, the neighbours are zero, and updateProbRescue will not do anything
        self.updateProbRescue()

    def checkSpindleDisassembly(self):
        """
        A function to check if the spindle disassembles (no contact between microtubules)
        :return:
        """

        # There must be microtubules oriented in both directions, and the lengths of the two longest microtubules of
        # each side have to be at least as long as the spindle

        # We create lists that contain the lengths of microtubules in each orientation
        lengths_plus1 = list()
        lengths_minus1 = list()

        for mt in self.microtubules:
            if mt.lost:
                continue
            mt_len = mt.pos * mt.orientation + self.spindle_edge
            if mt.orientation == 1:
                lengths_plus1.append(mt_len)
            else:
                lengths_minus1.append(mt_len)

        if len(lengths_plus1) and len(lengths_minus1) and (
                max(lengths_plus1) + max(lengths_minus1)) > self.spindle_edge * 2:
            return False
        else:
            return True

    def run(self):
        """
        Run the simulation, the return value is a string containing the information to be printed to a text file. Each
        line has the following information as fields in a csv file:
            1) the id of the microtubule
            2) the time of the simulation
            3) the position of the microtubule at the time of the simulation
            4) what happened at the simulation at this timepoint
                -1: simulation start
                 0: microtubule catastrophe
                 1: microtubule rescue
                 2: microtubule is lost
                 3: simulation end
            5) orientation of the microtubule
        :return:
        """
        self.t = 0

        # We run 20 minutes of simulation time
        while self.t < 20.:
            self.t += self.par.dt
            self.spindle_edge += self.par.dt * self.par.v_slide

            if self.checkSpindleDisassembly():
                # The spindle is lost, stop the simulation
                break

            if self.par.ase1:
                self.updateProbRescueAse1()

            for mt in self.microtubules:
                if not mt.lost:
                    mt.step()

        # Write the final timepoint
        for mt in self.microtubules:
            if not mt.lost:
                self.mt_coordinates += '%u %.2f %.2f 3 %i\n' % (mt.id, self.t, mt.pos, mt.orientation)

        return self.mt_coordinates


class Microtubule:

    def __init__(self, mt_id, starting_point, orientation, sim):
        """

        :param mt_id:
        :param starting_point:
        :param orientation:
        :param sim:
        :type sim: Simulation
        """
        # The id of the microtubule, which also coincides with its position in the list Simulation.microtubules
        self.id = mt_id

        # Initially the microtubule id and the neighbour_id coincide, but they could switch due to rearrangement
        self.neighbour_id = mt_id

        # Reference to the simulation object, to access the parameters
        self.sim = sim

        # Can be 1 (growing towards the right) or -1
        self.orientation = orientation

        # Position of the microtubule plus-end with respect to the spindle center
        self.pos = starting_point * self.orientation

        # The distance that the microtubule travels during dt when it is growing/shrinking
        self.step_grow = self.sim.par.dt * (self.sim.par.v_growth - self.sim.par.v_slide) * orientation
        self.step_shrink = -self.sim.par.dt * (self.sim.par.v_shrink + self.sim.par.v_slide) * orientation

        # Boolean indicating whether the microtubule is growing. They all start growing.
        self.growing = True

        # Time until the next catastrophe
        self.next_catastrophe = self.sim.timeToNextCatastrophe()

        # Whether the microtubule is lost
        self.lost = False

        # When the microtubule is created, we append it to the output of the simulation
        self.sim.mt_coordinates += '%u %.2f %.2f -1 %i\n' % (self.id, self.sim.t, self.pos, self.orientation)

    def countNeighbours(self):

        count = 0
        # The positions of occupied microtubules
        occupied_pos = [mt.neighbour_id for mt in self.sim.microtubules if not mt.lost]
        for neighbour in self.sim.neighbour_list[self.neighbour_id]:
            if neighbour in occupied_pos:
                count += 1

        return count

    def rescueProbBetaDistribution(self):

        # Convert the position of the plus end into an input for the beta function:
        # (-midzone_lenght, +midzone_length) -> (0,1)
        x_beta = (self.pos*self.orientation + self.sim.midzone_edge) / (2 * self.sim.midzone_edge)

        # Multiply the rate of rescue per number of microtubules by the value of the beta distribution, and return the
        # probability of rescue
        prob = beta.pdf(x_beta,self.sim.par.alpha,self.sim.par.beta)
        rate = self.sim.rate_rescue[self.id]*prob
        return 1.-exp(-rate*self.sim.par.dt)


    def rescueProb(self):
        """
        Return the rescue probability of the microtubule depending on the position.
        If the microtubule has depolymerised beyond the spindle pole, returns -1, to indicate that the
        microtubule is lost.
        Additionally:
            Wild-type:
                If inside the midzone: returns the probability depending on number of neighbours.
                Elsewhere: returns zero (No rescue outside the midzone)
            Ase1:
                The probability of rescue is the same everywhere inside the spindle
        :return:
        """

        # Test if microtubule depolymerises beyond the pole -> then it is lost
        if self.pos * self.orientation < -self.sim.spindle_edge:
            return -1

        # The constant probability of ase1 spindles
        if self.sim.par.ase1:
            return self.sim.prob_rescue[0]

        # The microtubule is inside the midzone
        if abs(self.pos) < self.sim.midzone_edge:
            if self.sim.par.alpha == 0:
                return self.sim.prob_rescue[self.id]
            else:
                return self.rescueProbBetaDistribution()
        else:
            return 0

    def step(self):

        # For microtubules that are growing, we substract dt from the time to next catastrophe, until the value is less
        # than zero, at which point we switch to shrinking (catastrophe).
        if self.growing:
            # The microtubule grows
            self.pos += self.step_grow
            self.next_catastrophe -= self.sim.par.dt

            # Catastrophe also occurs if the microtubules hits the pole
            if self.next_catastrophe < 0 or (self.pos*self.orientation) > self.sim.spindle_edge:
                self.growing = False
                # We print the catastrophe event to the simulation output
                self.sim.mt_coordinates += '%u %.2f %.2f 0 %i\n' % (self.id, self.sim.t, self.pos, self.orientation)
        else:

            # The microtubule shrinks
            self.pos += self.step_shrink

            # Test rescue

            # This function returns the probability of rescue during dt, or -1 if the microtubule depolymerises beyond
            # the spb (microtubule is lost)
            prob = self.rescueProb()

            if prob < 0:
                # Manage the consequences of losing the microtubule
                self.sim.addLostMicrotubule(self.id)
                # We print the loss event to the simulation output
                self.sim.mt_coordinates += '%u %.2f %.2f 2 %i\n' % (self.id, self.sim.t, self.pos, self.orientation)
            elif prob > random():
                self.next_catastrophe = self.sim.timeToNextCatastrophe()
                self.growing = True
                # We print the rescue event to the simulation output
                self.sim.mt_coordinates += '%u %.2f %.2f 1 %i\n' % (self.id, self.sim.t, self.pos, self.orientation)

