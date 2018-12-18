### General imports
import json
from tools import conVec, vector
import integrators

### Body class ###
class Body(object):
    __slots__ = ('mass', 'GM', 'velocity', 'position', 'temp_position', 'k', 'xv', 'radius', 'acc', 'name', 'index')
    def __init__(self, mass=1, GM=1, radius=1, velocity=vector(0,0,0), position=vector(0,0,0), name="Body", index=0):
        self.mass = mass
        self.GM = GM
        self.velocity = velocity
        self.position = position
        self.temp_position = vector(0,0,0)
        self.k = []
        self.xv = conVec(0,0)
        self.radius = radius
        self.acc = vector(0,0,0) #gravitational_acc(self.position, params) # approximate the initial acceleration for Verlet
        self.name = name
        self.index = index


### Parameters class ###
class parameters(object):
    __slots__ = ('G', 'AU', 'M', 'Theta', 'Epsilon', 'Lambda', 'Chi', 'time', 'dt', 'scale', 'rate', 'config', 'useconfig', 'integrator', 'printEndPos', 'checkEndPos', 'args',
                     'current_time', 'bodies', 'body_pairs', 'comets', 'all_bodies', 'start_time', 'end_time')
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
        self.bodies = []
        self.comets = []
        self.all_bodies = []
        self.body_pairs = []

        # load config from
        with open(args.configfile, "r") as configfile:
            config = json.load(configfile)
        self.config = config
        planet_config = config[0]
        config = config[1] 
        # load settings from configfile
        if args.useconfig:
            
            if 'dt' in config:
                self.dt = config['dt']
            else:
                self.dt = args.dt
            
            if 'time' in config:
                self.end_time = config['time']
            else:
                self.end_time = args.time
            
            if 'integrator' in config:
                self.integrator = config['integrator']
            else:
                self.integrator = args.integrator

        else:
            self.dt = args.dt
            self.end_time = args.time
            self.integrator = args.integrator
        
        # check which integrator was chosen
        if self.integrator.lower() == "euler":
            self.integrator = integrators.Euler
        elif self.integrator.lower() == "rk4":
            self.integrator = integrators.Runge_Kutta
        elif self.integrator.lower() == "verlet":
            self.integrator = integrators.Verlet
        elif self.integrator.lower() == "fr":
            self.integrator = integrators.Forest_Ruth
        elif self.integrator.lower() == "pefrl":
            self.integrator = integrators.PEFRL
        else:
            raise Exception('No valid integrator was provided')
        
        # load planets from configfile
        for i, body in enumerate(planet_config):
            # check to see if all needed fields are filled
            if not 'mass' in body and not 'gm' in body:
                raise Exception(f"Body nr. {i} don't have any mass properties")
            else:
                if 'gm' in body:
                    gm = body['gm']
                else:
                    gm = body['mass'] * self.G
                if 'mass' in body:
                    mass = body['mass']
                else:
                    mass = body['gm']/self.G

            if not 'radius' in body:
                raise Exception(f"Body nr. {i} don't have any radius property")
            else:
                radius = body['radius']

            if not 'position' in body:
                raise Exception(f"Body nr. {i} don't have any position property")
            else:
                position = body['position']

            if not 'velocity' in body:
                raise Exception(f"Body nr. {i} don't have any velocity property")
            else:
                velocity = body['velocity']

            if not 'name' in body:
                print(f"WARNING: Body nr. {i} don't have a name property, default 'Body{i}'' will be used.")
                name = f"Body{i}"
            else:
                name = body['name']

            if not 'comet' in body:
                print(f"WARNING: Body nr. {i} don't have a comet property, default 'false' will be used.")
                comet = False
            else:
                comet = body['comet']
            
            # add planet to corresponding list
            if comet:
                self.comets.append(Body(
                    name = name,
                    mass=mass,
                    GM=gm,
                    radius=radius,
                    position=vector(position[0], position[1], position[2]),
                    velocity=vector(velocity[0], velocity[1], velocity[2]),
                    index=len(self.bodies) + len(self.comets)
                ))
            else:
                self.bodies.append(Body(
                    name = name,
                    mass=mass,
                    GM=gm,
                    radius=radius,
                    position=vector(position[0], position[1], position[2]),
                    velocity=vector(velocity[0], velocity[1], velocity[2]),
                    index=len(self.bodies) + len(self.comets)
                ))
            
            self.all_bodies = self.bodies + self.comets

            # initialize acceleration for verlet
            for body in self.all_bodies:
                body.acc = integrators.gravitational_acc(body.position, self)
            ### END parameters.__init__() ###