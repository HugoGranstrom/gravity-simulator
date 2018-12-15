from vpython import *
import argparse
import json
from collections import namedtuple
import time
from math import pi;
from itertools import combinations



# UNITS:
# Mass: solar mass
# Length: Astronomical unit
# Time: days
# G = 4pi^2*AU^3/(M * 365.25) => G = 4*pi^2/365.25^2
#G = 6.67e-11
# AU^3/D^2 = 1/448485856027460.06  km^3/s^2 = 2.2297247205467538e-15 * km^3/s^2
#scale_factor = 1000
#dt = 0.01

### Constants ###

#G = 2.9592e-04
# G = 4*pi**2/365.25**2
# AU = 1.5e11
# M = 2e30
# Theta = 1/(2-2**(1/3))
# Epsilon = 0.1786178958448091
# Lambda = -0.2123418310626054
# Chi = -0.6626458266981849E-01

### containerVector ###
conVec = namedtuple("conVec","x y")

### Functions ###
# gravitational acceleration for Euler and Verlet
def gravitational_acc(position, p):
    sum_acc = vector(0,0,0)
    # calculate the gravitational acceleration from all other bodies
    for body in p.bodies:
        # distance to the other body
        r_vec = body.position - position
        r = r_vec.mag
        # skip if body is itself
        if r < body.radius:
            continue
        # the magnitude of the force
        acc = body.GM / r**2
        # the unit vector for the force
        dir = r_vec.hat
        # the force vector
        acc = acc * dir
        # add force vector to the sum of forces
        sum_acc += acc
    return sum_acc

# gravitational acceleration for Runge-Kutta
def gravitational_acc_runge(xv, p):
    sum_acc = vector(0,0,0)
    # calculate the gravitational acceleration from all other bodies
    for body in p.bodies:
        # distance to the other body
        r_vec = body.temp_position - xv.x
        r = r_vec.mag
        # skip if body is itself
        if r < body.radius:
            continue
        # the magnitude of the force
        acc = body.GM / r**2
        # the unit vector for the force
        dir = r_vec.hat
        # the force vector
        acc = acc * dir
        # add force vector to the sum of forces
        sum_acc += acc
    return conVec(xv.y, sum_acc)

### Integrators ###

# Euler Integrator
def Euler(p):
    # acceleration and velocity calculations
    for body in p.all_bodies:
        body.acc = gravitational_acc(body.position, p)
        body.velocity += p.dt * body.acc
    # position calculations
    for body in p.all_bodies:
        body.position += p.dt * body.velocity
        
        # update position of graphics
        body.sphere.pos = body.position
        body.label.pos = body.position
        
# Runge-Kutta 4 (RK4) Integrator
def Runge_Kutta(p):
    
    # calculate k1
    for body in p.all_bodies:
        body.k = [0] # zero at beginning to push indices one step up for readability
        body.temp_position = body.position
        body.xv = conVec(body.temp_position, body.velocity)
        temp_k1 = gravitational_acc_runge(body.xv, p)
        k1 = conVec(temp_k1.x * p.dt, temp_k1.y * p.dt)
        body.k.append(k1)
        
    # move temp_pos according to k1
    for body in p.all_bodies:
        body.temp_position = body.position + body.k[1].x/2
    
    # calculate k2
    for body in p.all_bodies:
        temp_xv = conVec(body.xv.x + body.k[1].x/2, body.xv.y + body.k[1].y/2)
        temp_k2 = gravitational_acc_runge(temp_xv, p)
        k2 = conVec(temp_k2.x * p.dt, temp_k2.y * p.dt)
        body.k.append(k2)
    
    # move temp_pos according to k2
    for body in p.all_bodies:
        body.temp_position = body.position + body.k[2].x/2

    # calculate k3
    for body in p.all_bodies:
        temp_xv = conVec(body.xv.x + body.k[2].x/2, body.xv.y + body.k[2].y/2)
        temp_k3 = gravitational_acc_runge(temp_xv, p)
        k3 = conVec(temp_k3.x * p.dt, temp_k3.y * p.dt)
        body.k.append(k3)
    
    # move temp_pos according to k3
    for body in p.all_bodies:
        body.temp_position = body.position + body.k[3].x

    # calculate k4
    for body in p.all_bodies:
        temp_xv = conVec(body.xv.x + body.k[3].x, body.xv.y + body.k[3].y)
        temp_k4 = gravitational_acc_runge(temp_xv, p)
        k4 = conVec(temp_k4.x * p.dt, temp_k4.y * p.dt)
        body.k.append(k4)

    # calculate weighted sum and move bodies
    for body in p.all_bodies:
        body.position += 1/6 * (body.k[1].x + 2*body.k[2].x + 2*body.k[3].x + body.k[4].x)
        body.velocity += 1/6 * (body.k[1].y + 2*body.k[2].y + 2*body.k[3].y + body.k[4].y)

        body.sphere.pos = body.position
        body.label.pos = body.position

# Velocity-Verlet Integrator
def Verlet(p):
    # position calculations
    for body in p.all_bodies:
        body.position += body.velocity * p.dt + body.acc/2 * p.dt**2
        body.sphere.pos = body.position
        body.label.pos = body.position
    # acceleration and velocity calculations
    for body in p.all_bodies:
        temp_acc = gravitational_acc(body.position, p)
        body.velocity += p.dt/2*(body.acc + temp_acc)
        body.acc = temp_acc

def Forest_Ruth(p):
    for body in p.all_bodies:
        body.position += p.Theta*p.dt/2*body.velocity
    for body in p.all_bodies:
        body.velocity += p.Theta*p.dt*gravitational_acc(body.position, p)
    for body in p.all_bodies:
        body.position += (1-p.Theta)*p.dt/2*body.velocity
    for body in p.all_bodies:
        body.velocity += (1-2*p.Theta)*p.dt*gravitational_acc(body.position, p)
    for body in p.all_bodies:
        body.position += (1-p.Theta)*p.dt/2*body.velocity
    for body in p.all_bodies:
        body.velocity += p.Theta*p.dt*gravitational_acc(body.position, p)
    for body in p.all_bodies:
        body.position += p.Theta*p.dt/2*body.velocity
        body.sphere.pos = body.position
        body.label.pos = body.position


def PEFRL(p):
    for body in p.all_bodies:
        body.position += p.Epsilon*p.dt*body.velocity
    for body in p.all_bodies:
        body.velocity += (1-2*p.Lambda)*p.dt/2*gravitational_acc(body.position, p)
    for body in p.all_bodies:
        body.position += p.Chi*p.dt*body.velocity
    for body in p.all_bodies:
        body.velocity += p.Lambda*p.dt*gravitational_acc(body.position, p)
    for body in p.all_bodies:
        body.position += (1-2*(p.Chi+p.Epsilon))*p.dt*body.velocity
    for body in p.all_bodies:
        body.velocity += p.Lambda*p.dt*gravitational_acc(body.position, p)
    for body in p.all_bodies:
        body.position += p.Chi*p.dt*body.velocity
    for body in p.all_bodies:
        body.velocity += (1-2*p.Lambda)*p.dt/2*gravitational_acc(body.position, p)
    for body in p.all_bodies:
        body.position += p.Epsilon*p.dt*body.velocity
        body.sphere.pos = body.position
        body.label.pos = body.position

### End Integrators ###



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


def color_to_vector(color_list):
    return vector(color_list[0]/255, color_list[1]/255, color_list[2]/255)

def onClick(e):
        obj = scene.mouse.pick  # get the sphere
        if(obj != None):
            body = p.bodies[obj.index]  # each sphere has a index attribute wich points to the planets position in the bodies list
            # info about the planets can noe be retrieved from the class
            info_label.text = '<b><i>'+body.name+'</i></b>\n<i>Mass:</i> '+str(body.mass)+' M☉\n' + '<i>Position:</i> '+str(body.position)+'\n'
            # TODO add more info about planets & convert units
            scene.camera.follow(obj)
        else:
            info_label.text = ''
            #scene.camera.follow()



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
                self.scale_factor = args.scale_factor

        else:
            self.dt = args.dt
            self.end_time = args.time
            self.integrator = args.integrator
            self.scale_factor = args.scale
        
        # check which integrator was chosen
        if self.integrator.lower() == "euler":
            self.integrator = Euler
        elif self.integrator.lower() == "rk4":
            self.integrator = Runge_Kutta
        elif self.integrator.lower() == "verlet":
            self.integrator = Verlet
        elif self.integrator.lower() == "fr":
            self.integrator = Forest_Ruth
        elif self.integrator.lower() == "pefrl":
            self.integrator = PEFRL
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
            
            self.body_pairs = combinations(self.bodies, 2)
            self.all_bodies = self.bodies + self.comets

            ### END parameters.__init__() ###

# argument parsing
parser = argparse.ArgumentParser(description="A newtonian gravity simulator")
parser.add_argument("-t", "--time", type=float, default=0, dest="time",
            help="The amount of time the simulation will simulate measured in days. Type 0 for infinite time. (Default: 0)")
parser.add_argument("--dt", type=float, default=0.1, dest="dt",
            help="The timestep to use. (Default: 0.1)")
parser.add_argument("--scale", type=float, default=1000, dest="scale",
            help="The number to scale the radiuses of the planets to make them visible. Does only affect the visuals not collisions. (Default: 1000)")
parser.add_argument("--rate", type=int, default=100000, dest="rate",
            help="Number of timesteps per second (Default: 100000)")
parser.add_argument("--configfile", type=str, default="config.json", dest="configfile",
            help="Path to the config file containing the bodies. (Default: config.json)")
parser.add_argument("--useconfig", action="store_true", default=False, dest="useconfig",
            help="Use this flag if you want to use the settings in the configfile instead of defaults and cmd arguments. (Default: False)")
parser.add_argument("--integrator", type=str, default="verlet", dest="integrator", 
            help="The integrator to be used. Options: euler, verlet, rk4, fr, pefrl (Default: verlet)")
parser.add_argument("--endPos", action="store_true", default=False, dest="printEndPos",
            help="When flagged the end position of all bodies will be printed (Deafult: False)")
parser.add_argument("--checkEndPos", action="store_true", default=False, dest="checkEndPos",
            help="When flagged the end position of all bodies is compared to their real end positions, which are given as 'end_position' in config.json (Default: False)")

# parse cmd arguments
args = parser.parse_args()

p = parameters(
    args=args,
    G=4*pi**2/365.25**2,
    AU = 1.5e11,
    M = 2e30,
    Theta = 1/(2-2**(1/3)),
    Epsilon = 0.1786178958448091,
    Lambda = -0.2123418310626054,
    Chi = -0.6626458266981849E-01,
)

# initialize acceleration for verlet and scale spheres
for body in p.all_bodies:
    body.acc = gravitational_acc(body.position, p)
    if body.scale:
        body.sphere.radius *= p.scale_factor

time_label = label(pos=vector(20, 350, 0), pixel_pos=True, align='left', text="Time: " + str(p.current_time/365) + " years")

### info box ###
info_label = label(pos=vector(20, scene.height/3, 0), pixel_pos=True, box=False, align='left', text="")
scene.bind('click', onClick)
### ###

# loop over every body and run its update method every timestep
p.start_time = time.time()
if p.end_time > 0:
    for epoch in range(int(p.end_time/p.dt)):
        rate(p.args.rate)

        p.integrator(p)

        p.current_time = epoch*p.dt
        time_label.text = "Time: {:.2f} years".format(p.current_time/365)
    # print body positions for benchmarking
    if p.printEndPos:
        for body in p.bodies:
            print(f"{body.name}: {body.position}")
    if args.checkEndPos:
        error_sum = 0
        for body in p.bodies:
            try:
                end_pos = p.config[0][body.index]["end_position"]
                end_pos = vector(end_pos[0], end_pos[1], end_pos[2])
                error = mag(body.position - end_pos) # the magnitude of the error
                error_sum += error
                print(f"{body.name}: {error} AU")
            except:
                pass
        print(f"Total error: {error_sum}")
        print(f"dt: {p.dt}")
        print(f"Integrator: {args.integrator}")



else:
    while True:
        rate(p.args.rate)

        p.integrator(p)

        time_label.text = "Time: {:.2f} years".format(p.current_time/365)
        p.current_time += p.dt

print(f"Execution time: {time.time()-p.start_time} seconds")
