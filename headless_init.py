### General imports
import json
from tools import conVec, vector
import integrators
import random

### Body class ###
class Body(object):
    __slots__ = (
        "mass",
        "GM",
        "velocity",
        "position",
        "temp_position",
        "k",
        "xv",
        "radius",
        "acc",
        "name",
        "index",
        "old_positions",
        "old_acc",
        "temp_velocity",
        "energy"
    )

    def __init__(
        self,
        mass=1,
        GM=1,
        radius=1,
        velocity=vector(0, 0, 0),
        position=vector(0, 0, 0),
        name="Body",
        index=0,
    ):
        self.mass = mass
        self.GM = GM
        self.velocity = velocity
        self.position = position
        self.temp_position = position
        self.k = []
        self.xv = conVec(0, 0)
        self.radius = radius
        self.acc = vector(
            0, 0, 0
        )  # gravitational_acc(self.position, params) # approximate the initial acceleration for Verlet
        self.name = name
        self.index = index
        self.old_positions = [vector(0, 0, 0), vector(0, 0, 0), self.position]


### Parameters class ###
class parameters(object):
    __slots__ = (
        "G",
        "AU",
        "M",
        "Theta",
        "Epsilon",
        "Lambda",
        "Chi",
        "time",
        "dt",
        "scale",
        "rate",
        "config",
        "useconfig",
        "integrator",
        "printEndPos",
        "checkEndPos",
        "args",
        "current_time",
        "bodies",
        "body_pairs",
        "comets",
        "all_bodies",
        "start_time",
        "end_time",
        "adaptive",
        "next_dt",
        "old_times",
    )

    def __init__(self, args, G, AU, M, Theta, Epsilon, Lambda, Chi):

        self.G = G
        self.AU = AU
        self.M = M
        self.Theta = Theta
        self.Epsilon = Epsilon
        self.Lambda = Lambda
        self.Chi = Chi

        self.args = args
        self.printEndPos = args.printEndPos
        self.checkEndPos = args.checkEndPos

        self.current_time = 0
        self.old_times = [0, 0, 0]
        self.bodies = []
        self.comets = []
        self.all_bodies = []
        self.body_pairs = []

        # load config from
        with open(args.configfile, "r") as configfile:
            config = json.load(configfile)
        self.config = config
        planet_config = config[0]
        try:
            config = config[1]
        except:
            config = None
        # load settings from configfile
        if args.useconfig:

            if "dt" in config:
                self.dt = config["dt"]
            else:
                self.dt = args.dt

            if "time" in config:
                self.end_time = config["time"]
            else:
                self.end_time = args.time

            if "integrator" in config:
                self.integrator = config["integrator"]
            else:
                self.integrator = args.integrator

        else:
            self.dt = args.dt
            self.end_time = args.time
            self.integrator = args.integrator
        self.next_dt = self.dt
        # check which integrator was chosen
        if self.integrator.lower() == "euler":
            self.integrator = integrators.Euler
            self.adaptive = False
        elif self.integrator.lower() == "rk4":
            self.integrator = integrators.Runge_Kutta
            self.adaptive = False
        elif self.integrator.lower() == "verlet":
            self.integrator = integrators.Verlet
            self.adaptive = False
        elif self.integrator.lower() == "fr":
            self.integrator = integrators.Forest_Ruth
            self.adaptive = False
        elif self.integrator.lower() == "pefrl":
            self.integrator = integrators.PEFRL
            self.adaptive = False
        elif self.integrator.lower() == "rk45":
            self.integrator = integrators.Runge_Kutta_45
            self.adaptive = True
        elif self.integrator.lower() == "nystrom2":
            self.integrator = integrators.Nystrom2
            self.adaptive = False
        elif self.integrator.lower() == "nystrom3":
            self.integrator = integrators.Nystrom3
            self.adaptive = False
        elif self.integrator.lower() == "nystrom4":
            self.integrator = integrators.Nystrom4
            self.adaptive = False
        elif self.integrator.lower() == "nystrom5":
            self.integrator = integrators.Nystrom5
            self.adaptive = False
        elif self.integrator.lower() == "nystrom6":
            self.integrator = integrators.Nystrom6
            self.adaptive = False
        elif self.integrator.lower() == "nystromsimos":
            self.integrator = integrators.NystromSimos
            self.adaptive = False
        elif self.integrator.lower() == "rkn64":
            self.integrator = integrators.RKN64
            self.adaptive = True
        elif self.integrator.lower() == "rkn12":
            self.integrator = integrators.RKN12
            self.adaptive = True
        elif self.integrator.lower() == "nystrom11":
            self.integrator = integrators.Nystrom11
            self.adaptive = False
        elif self.integrator.lower() == "nystrom8":
            self.integrator = integrators.Nystrom8
            self.adaptive = False
        elif self.integrator.lower() == "nystrom12":
            self.integrator = integrators.Nystrom12
            self.adaptive = False
        elif self.integrator.lower() == "sympny10":
            self.integrator = integrators.SympNystrom10
            self.adaptive = False
        elif self.integrator.lower() == "yoshida6":
            self.integrator = integrators.Yoshida6
            self.adaptive = False
        elif self.integrator.lower() == "yoshida8":
            self.integrator = integrators.Yoshida8
            self.adaptive = False
        elif self.integrator.lower() == "kahanli8":
            self.integrator = integrators.KahanLi8
            self.adaptive = False
        elif self.integrator.lower() == "mcate5":
            self.integrator = integrators.McAte5
            self.adaptive = False
        elif self.integrator.lower() == "kuraev":
            self.integrator = integrators.Kuraev
            self.adaptive = False
        elif self.integrator.lower() == "beeman":
            self.integrator = integrators.Beeman
            self.adaptive = False
        elif self.integrator.lower() == "opverlet":
            self.integrator = integrators.Optimized_Verlet
            self.adaptive = False
        elif self.integrator.lower() == "bh":
            self.integrator = integrators.PEFRL_BH
            self.adaptive = False
        else:
            raise Exception("No valid integrator was provided")

        # load planets from configfile
        for i, body in enumerate(planet_config):
            # check to see if all needed fields are filled
            if not "mass" in body and not "gm" in body:
                raise Exception(f"Body nr. {i} don't have any mass properties")
            else:
                if "gm" in body:
                    gm = body["gm"]
                else:
                    gm = body["mass"] * self.G
                if "mass" in body:
                    mass = body["mass"]
                else:
                    mass = body["gm"] / self.G

            if not "radius" in body:
                raise Exception(f"Body nr. {i} don't have any radius property")
            else:
                radius = body["radius"]

            if not "position" in body:
                raise Exception(f"Body nr. {i} don't have any position property")
            else:
                position = body["position"]

            if not "velocity" in body:
                raise Exception(f"Body nr. {i} don't have any velocity property")
            else:
                velocity = body["velocity"]

            if not "name" in body:
                print(
                    f"WARNING: Body nr. {i} don't have a name property, default 'Body{i}'' will be used."
                )
                name = f"Body{i}"
            else:
                name = body["name"]

            if not "comet" in body:
                print(
                    f"WARNING: Body nr. {i} don't have a comet property, default 'false' will be used."
                )
                comet = False
            else:
                comet = body["comet"]

            # add planet to corresponding list
            if comet:
                self.comets.append(
                    Body(
                        name=name,
                        mass=mass,
                        GM=gm,
                        radius=radius,
                        position=vector(position[0], position[1], position[2]),
                        velocity=vector(velocity[0], velocity[1], velocity[2]),
                        index=len(self.bodies) + len(self.comets),
                    )
                )
            else:
                self.bodies.append(
                    Body(
                        name=name,
                        mass=mass,
                        GM=gm,
                        radius=radius,
                        position=vector(position[0], position[1], position[2]),
                        velocity=vector(velocity[0], velocity[1], velocity[2]),
                        index=len(self.bodies) + len(self.comets),
                    )
                )
        """
        # add N random comets
        for i in range(10000):
            self.bodies.append(Body(
                name = "",
                mass=1.1e-16,
                GM=1.1e-16*self.G,
                radius=6.67e-06,
                position=vector(random.uniform(-50, 50), random.uniform(-50, 50), random.uniform(-0.01, 0.01)),
                velocity=vector(random.uniform(-1e-2, 1e-2), random.uniform(-1e-2, 1e-2), random.uniform(-1e-4, 1e-4)),
                index=len(self.bodies) + len(self.comets),
            ))
        """
        self.all_bodies = self.bodies + self.comets

        # initialize acceleration for verlet
        for body in self.all_bodies:
            body.acc = integrators.gravitational_acc(body.position, self)
        # initialize for leapfrogs
        for body in self.all_bodies:
            body.temp_position = body.position - body.velocity * self.dt - 0.5 * body.acc * self.dt ** 2
        for body in self.all_bodies:
            body.old_acc = integrators.gravitational_acc_temp(body.temp_position, self)
        ### END parameters.__init__() ###
