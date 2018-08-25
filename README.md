# Gravity-simulator
A newtonian brute force gravity simulator of the solar system written in Python/Vpython. It is probably not too accurate so don't plan your rocket launch using it. It can be fun though to play around making your own solar systems (or throwing in a black hole in our own).

# Usage
To run the default configuration (dt=0.01, integrator=euler, infinite loop, the sun and 9 planets) run ``python gravity.py`` to run the simulation. It will open up a browser window and it will show the simulation there. To stop the simulation you close the browser tab, CTRL-C does not work. You can navigate the scene with these controls:
- Zoom: scroll while your mouse pointer is above the simulation window.
- Pan: SHIFT + LEFT MOUSE BUTTON.
- Rotate: hold RIGHT MOUSE BUTTON
- Resize window: drag in the lower left corner of the simulation. 

You can also run the simulation in a headless mode ie. no gui. To use headless mode you run ``python gravity.py``. It has the same arguments as ``gravity.py``. In my benchmarking it was around 20-30% faster than the gui version. 

## Options
There are multiple command line arguments you can pass to the simulation:
- `-t` or `--time`, the length of the simulation measured in days. If it is set to ``0`` it will go on forever. Default: 0.
- `--dt`, the timestep the simulation will use measured in days. Default: 0.01.
- `--scale`, to make the simulation look better you can scale all planets with the property `"scale" = true` by a scale factor you specify here. Default: 1000.
- `--rate`, the maximum amount of timesteps per second. Default: 100000.
- `--configfile`, the path to the json file containing your configurations. Default: ``config.json``.
- `--useconfig`, if this is checked then the settings in your configuration file will be used over those you pass in the command line. Not checked by default.
- `--integrator`, the integrator you want to use. Default: ``"euler"``.
- `--endPos`, prints the end position of all planets if checked. Not checked by default.
- `--checkEndPos`, when simulation is ended, compare the end positions of all bodies with their ``"end_position"`` in the configuration file and prints both the individual error and the sum. Used for measuring accuracy of integrators. 

## Units
The units used in the simulation is Astronomical Units (AU) for length, Solar masses (M☉) for mass and Days for time. 

## config.json
The config.json file (you can choose another name if you like) is where you store the initial position, initial velocities, mass, name etc of your planets and bodies and optionally parameters for the simulation which can be used by adding ``--useconfig`` to the run-command. It is structured like this:

``` javascript
[
    [
        // list of planet objects 
        {
            "name": "Earth",
            "mass": 3e-06,
            "radius": 4.247333333333333e-05,
            "position": [-5.111702950987252E-01, -8.734341386147972E-01, 3.902531498407046E-05],
            "velocity": [1.457401965494037E-02, -8.749957786090569E-03, -3.393201214360642E-07],
            "end_position": [1,1,0],
            "color": [0, 255, 0],
            "trail": true,
            "scale": true
        },
        {
            "name": "Mercury",
            "mass": 1.65e-07,
            "radius": 4e-05,
            "position": [3.600062387312980E-01, -8.310671431721671E-02, -3.981766501010686E-02],
            "velocity": [8.732371820239134E-04, 2.867508157942580E-02, 2.263026727476856E-03],
            "end_position": [0,0,0],
            "color": [255, 0, 0],
            "trail": true,
            "scale": true
        }
        // more planets here
    ],
    [
        "dt": 0.01,
        "integrator": "euler"
        // more settings here
    ]
]
```

#### Planets
Every planet is configured in ``config.json`` by default. The current attributes that a planet has is:
- ``name`` - the name that will be displayed over the planet.
- ``mass`` - the mass of the planet measured in solar masses.
- ``radius`` - the radius of the planet measured in AU.
- ``position`` - the initial position vector. It is a 3D-vector and is given as a array/list: [x, y, z].
- ``velocity`` - the initial velocity vector. It is a 3D-vector and is given as a array/list: [x, y, z].
- ``end_position`` - used for benchmarking. Is the actual position of the planet when the simulation has ended.  
- ``color`` - the color the planet and its trail will have. Given as a array/list of RGB values in the range(0, 255): [red, blue, green].
- ``trail`` - set to ``true`` if the planet should show a trail after itself. If not set to ``false``.
- ``scale`` - set to ``true`` if you want to scale the visual radius (the one used in calculations is still the same as ``radius``) of the planet by ``scale_factor`` to make it easier to see. For big planets (eg the sun in a solar system) set it to ``false``.

#### Simulation Settings
If you want to share your simulation with someone else to let them run it on their computers it is handy to set your simulation settings in your ``config.json`` file because then all you have to do is to send them this file and they are ready to go. By default it is turned off, to use it, add ``--useconfig`` to your run-command like this: ``python gravity.py --useconfig``. Settings in your ``config.json`` file will have higher priority than those you pass through the command line. 


## Integrators
TODO

## Benchmark
The benchmark was done with this command (change {dt} and {integrator} with all the combinations you want to test):

``python headless.py -t 365 --configfile 22-07-2018-365-days.json --checkEndPos --dt {dt} --integrator {integrator}``

This will test how accurate the simulation is when it simulates one year (365 days). The data in the json file is from [Horizons](https://ssd.jpl.nasa.gov/horizons.cgi#top). 

### Results
|dt   | Euler  | Verlet | RK4     | FR     | PEFRL  |
|-----|-------:|-------:|--------:|-------:|-------:|
|10   |1.483818|1.521156|42.248875|1.153413|0.549528|
|5    |0.797699|0.414852|12.610780|0.133126|0.013813|
|2    |0.243673|0.166049|6.488631 |0.097902|0.093413|
|1    |0.095696|0.026149|3.461871 |0.009192|0.008896|
|0,1  |0.016178|0.009059|0.347950 |0.008889|0.008889|
|0,01 |0.009603|0.008890|0.037238 |0.008889|0.008889|
|0,001|0.008959|0.008889|0.010172 |0.008889|0.008889|

dt is the timestep measured in days and the numbers under the integrators are the total summed error between the simulated positions of the planets at the end of the simulation and the actual position gathered from [Horizons](https://ssd.jpl.nasa.gov/horizons.cgi#top).

# TODO
- Make documentation
- Benchmark: compare execution time and accuracy 

## Contribution
Feel free to submit pull request for improvements or additions. It makes it easier if you open an Issue first so we can discuss a good implementation of your feature.
Add your name to the ``CONTRIBUTORS.md`` file. Your name and what you added. 
