# gravity-simulator
A newtonian brute force gravity simulator of the solar system written in Python/Vpython. 

# Usage
To run the default configuration (dt=0.01, integrator=euler, infinite loop, the sun and 9 planets) run ``python gravity.py`` to run the simulation.

## Units
The units used in the simulation is Astronomical Units (AU) for length, Solar masses (M☉) for mass and Days for time. 

## config.json

### Planets
Every planet is configured in ``config.json`` by default. The current attributes that a planet has is:
- ``name`` - the name that will be displayed over the planet.
- ``mass`` - the mass of the planet measured in solar masses.
- ``radius`` - the radius of the planet measured in AU.
- ``position`` - the initial position vector. It is a 3D-vector and is given as a array/list: [x, y, z].
- ``velocity`` - the initial velocity vector. It is a 3D-vector and is given as a array/list: [x, y, z].
- ``color`` - the color the planet and its trail will have. Given as a array/list of RGB values in the range(0, 255): [red, blue, green].
- ``trail`` - set to ``true`` if the planet should show a trail after itself. If not set to ``false``.
- ``scale`` - set to ``true`` if you want to scale the visual radius (the one used in calculations is still the same as ``radius``) of the planet by ``scale_factor`` to make it easier to see. For big planets (eg the sun in a solar system) set it to ``false``.

### Simulation Settings


# TODO
- Make documentation
- Runge-Kutta Integrator
- Benchmark: compare execution time and accuracy 

## Contribution
Feel free to submit pull request for improvements or additions. It makes it easier if you open an Issue first so we can discuss a good implementation of your feature.
Add your name to the ``CONTRIBUTORS.md`` file. Your name and what you added. 
