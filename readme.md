# Simulations of Lera-Ramírez et al. 2021

This is a simulation of the microtubule dynamics in the anaphase spindle, used in Lera Ramirez 2021 (unpublished).

## How to run the simulation

You will need python 3, with the following libraries:

```
numpy
scipy
matplotlib

# If you want to run the simulations in parallel you also need
joblib
```

## Structure of the simulation

All the positions are recorded with respect to the spindle center. All the magnitudes are represented in micrometers and minutes. There are three main classes:

### Parameters

Contains the parameters of the simulation, see the `Parameters::__init__` function for what they are.

### Simulation class

It's initialised with an instance of `Parameters`, see `Simulation::init`. When initialized:

* Sets the position of spindle poles at -2 and 2 (initial spindle length of 4 um).
* Creates 9 microtubules oriented like in the EM checkerboard, of length 3 um each (position of plus end -1 or +1 with respect to center).
* All microtubules are initially in the growing state.

### Microtubule class

It's initialised with an id, initial position of the plus end, orientation (+1/-1) and with a reference to the `Simulation` instance, see `Microtubule::init`.

### Simulation algorithm

See `Simulation::run`

### Managing rearrangements

This can be a bit confusing. Essentially every microtubule has an `id` and a `grid_position` that initially match. The `grid_position` corresponds to the position on the ZY axis:

```
Grid positions:
	0 -- 1 -- 2
	|    |    |
	3 -- 4 -- 5
	|    |    |
	6 -- 7 -- 8
```

`Simulation::neighbour_list` is a list used to know the `grid_position` of all neighbours of a microtubule.
Each position in the list corresponds to a `grid_position`, and contains the `grid_position` of neighbouring microtubules.

For instance, the `grid_position` of neighbours of a microtubule with `grid_position=4` are `[1,3,5,7]`. This is why `Simulation::neighbour_list` is:

```python
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
```

When a microtubule is lost, there can be a rearrangement to maximize the number of neighbours. See the call to `Simulation::performRearrangement` in `Simulation::addLostMicrotubule`. In a rearrangement, the `grid_position` of two microtubules are swapped.

## How to run the simulation

See `example_simulation.py`
