from vpython import color, label, sphere, vector
from tools import conVec
import json
import integrators


def color_to_vector(color_list):
    return vector(color_list[0]/255, color_list[1]/255, color_list[2]/255)

class Body():
    __slots__ = ('mass', 'GM', 'velocity', 'position', 'temp_position', 'k', 'xv', 'radius', 'acc', 'name', 'index', 'color', 'label', 'sphere', 'scale')
    def __init__(self, mass=1, GM=1, radius=1, velocity=vector(0,0,0), position=vector(0,0,0), color=color.white, trail=True, name="Body", scale=True, index=0):
        self.mass = mass
        self.GM = GM
        self.velocity = velocity
        self.position = position
        self.temp_position = vector(0,0,0)
        self.k = []
        self.xv = conVec(0,0)
        #self.sum_forces = vector(0,0,0)
        self.color = color
        self.radius = radius
        self.acc = vector(0,0,0) #gravitational_acc(self.position) # approximate the initial acceleration for Verlet
        self.name = name
        self.label = label(pos=self.position, text=self.name, height=10)
        self.index = index
        self.sphere = sphere(pos=self.position, color=self.color, radius=self.radius, make_trail=trail, retain=200, index=self.index)
        self.scale = scale


class parameters(object):
    __slots__ = ('G', 'AU', 'M', 'Theta', 'Epsilon', 'Lambda', 'Chi', 'time', 'dt', 'scale', 'rate', 'config', 'useconfig', 'integrator', 'printEndPos', 'checkEndPos', 'scale_factor','args',
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
            
            if 'scale_factor' in config:
                self.scale_factor = config['scale_factor']
            else:
                self.scale_factor = args.scale

        else:
            self.dt = args.dt
            self.end_time = args.time
            self.integrator = args.integrator
            self.scale_factor = args.scale
        
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
            
            if not 'trail' in body:
                print(f"WARNING: Body nr. {i} don't have a trail property, default 'true' will be used.")
                trail = True
            else:
                trail = body['trail']

            if not 'color' in body:
                print(f"WARNING: Body nr. {i} don't have a color property, default [255, 255, 255] will be used.")
                color = [255, 255, 255]
            else:
                color = body['color']
            
            if not 'scale' in body:
                print(f"WARNING: Body nr. {i} don't have a scale property, default 'true' will be used.")
                scale = True
            else:
                scale = body['scale']

            # add planet to corresponding list
            if comet:
                self.comets.append(Body(
                    name = name,
                    mass=mass,
                    GM=gm,
                    radius=radius,
                    position=vector(position[0], position[1], position[2]),
                    velocity=vector(velocity[0], velocity[1], velocity[2]),
                    index=len(self.bodies) + len(self.comets),
                    trail=trail,
                    color=color_to_vector(color),
                    scale=scale,
                ))
            else:
                self.bodies.append(Body(
                    name = name,
                    mass=mass,
                    GM=gm,
                    radius=radius,
                    position=vector(position[0], position[1], position[2]),
                    velocity=vector(velocity[0], velocity[1], velocity[2]),
                    index=len(self.bodies) + len(self.comets),
                    trail=trail,
                    color=color_to_vector(color),
                    scale=scale,
                ))
            
        self.all_bodies = self.bodies + self.comets
        # initialize acceleration for verlet
        for body in self.all_bodies:
            body.acc = integrators.gravitational_acc(body.position, self)
            if body.scale:
                body.sphere.radius *= self.scale_factor
            ### END parameters.__init__() ###
