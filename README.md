# gravity-simulator
A newtonian brute force gravity simulator of the solar system written in Python/Vpython. 

# Usage
To run the default configuration (dt=0.01, integrator=euler, infinite loop, the sun and 9 planets) run ``python gravity.py`` to run the simulation. It will open up a browser window and it will show the simulation there. You can navigate the scene with these controls:
- Zoom: scroll while your mouse pointer is above the simulation window.
- Pan: SHIFT + LEFT MOUSE BUTTON.
- Rotate: hold RIGHT MOUSE BUTTON
- Resize window: drag in the lower left corner of the simulation. 

## Options
There are multiple command line arguments you can pass to the simulation:
- `-t` or `--time`, the length of the simulation measured in days.
- `--dt`, the timestep the simulation will use measured in days.
- `--scale`, to make the simulation look better you can scale all planets with the property `"scale" = true` by a scale factor you specify here.
- `--rate`,
- `--configfile`,
- `--useconfig`,
- `--integrator`,
- `--endPos`,
- `--checkEndPos`,

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
- ``color`` - the color the planet and its trail will have. Given as a array/list of RGB values in the range(0, 255): [red, blue, green].
- ``trail`` - set to ``true`` if the planet should show a trail after itself. If not set to ``false``.
- ``scale`` - set to ``true`` if you want to scale the visual radius (the one used in calculations is still the same as ``radius``) of the planet by ``scale_factor`` to make it easier to see. For big planets (eg the sun in a solar system) set it to ``false``.

#### Simulation Settings
If you want to share your simulation with someone else to let them run it on their computers it is handy to set your simulation settings in your ``config.json`` file because then all you have to do is to send them this file and they are ready to go. By default it is turned off, to use it, add ``--useconfig`` to your run-command like this: ``python gravity.py --useconfig``. Settings in your ``config.json`` file will have higher priority than those you pass through the command line. 

# TODO
- Make documentation
- Benchmark: compare execution time and accuracy 
- Implement Forest Ruth (FR) and Position
Extended Forest-Ruth Like (PEFRL) integrators

## Contribution
Feel free to submit pull request for improvements or additions. It makes it easier if you open an Issue first so we can discuss a good implementation of your feature.
Add your name to the ``CONTRIBUTORS.md`` file. Your name and what you added. 
