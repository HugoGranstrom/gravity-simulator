from vpython import *
import argparse
import json
from collections import namedtuple

# argument parsing
parser = argparse.ArgumentParser(description="A newtonian gravity simulator")
parser.add_argument("-t", "--time", type=float, default=0, dest="time",
                    help="The amount of time the simulation will simulate measured in days. Type 0 for infinite time. (Default: 0)")
parser.add_argument("--dt", type=float, default=0.01, dest="dt",
                    help="The timestep to use. (Default: 0.01)")
parser.add_argument("--scale", type=float, default=1000, dest="scale",
                    help="The number to scale the radiuses of the planets to make them visible. Does only affect the visuals not collisions. (Default: 1000)")
parser.add_argument("--rate", type=int, default=100000, dest="rate",
                    help="Number of timesteps per second (Default: 100000)")
parser.add_argument("--configfile", type=str, default="config.json", dest="configfile",
                    help="Path to the config file containing the bodies. (Default: config.json)")
parser.add_argument("--useconfig", action="store_true", default=False, dest="useconfig",
                    help="Use this flag if you want to use the settings in the configfile instead of defaults and cmd arguments. (Default: False)")
parser.add_argument("--integrator", type=str, default="euler", dest="integrator", 
                    help="The integrator to be used. Options: euler, verlet, rk4, fr, pefrl (Default: euler)")
parser.add_argument("--endPos", action="store_true", default=False, dest="printEndPos",
                    help="When flagged the end position of all bodies will be printed (Deafult: False)")
parser.add_argument("--checkEndPos", action="store_true", default=False, dest="checkEndPos",
                    help="When flagged the end position of all bodies is compared to their real end positions, which are given as 'end_position' in config.json (Default: False)")


### containerVector ###
conVec = namedtuple("conVec","x y")

### Functions ###
# gravitational acceleration for Euler and Verlet
def gravitational_acc(position):
    sum_acc = vector(0,0,0)
    # calculate the gravitational acceleration from all other bodies
    for body in bodies:
        # distance to the other body
        r = mag(position-body.position)
        # skip if body is itself
        if r < body.radius:
            continue
        # the magnitude of the force
        acc = G * body.mass / r**2
        # the unit vector for the force
        dir = norm(body.position - position)
        # the force vector
        acc = acc * dir
        # add force vector to the sum of forces
        sum_acc += acc
    return sum_acc

# gravitational acceleration for Runge-Kutta
def gravitational_acc_runge(xv):
    sum_acc = vector(0,0,0)
    # calculate the gravitational acceleration from all other bodies
    for body in bodies:
        # distance to the other body
        r = mag(xv.x-body.temp_position)
        # skip if body is itself
        if r < body.radius:
            continue
        # the magnitude of the force
        acc = G * body.mass / r**2
        # the unit vector for the force
        dir = norm(body.temp_position - xv.x)
        # the force vector
        acc = acc * dir
        # add force vector to the sum of forces
        sum_acc += acc
    return conVec(xv.y, sum_acc)

### Integrators ###

# Euler Integrator
def Euler():
    # acceleration and velocity calculations
    for body in bodies:
        body.acc = gravitational_acc(body.position)
        body.velocity += dt * body.acc

    # position calculations
    for body in bodies:
        body.position += dt * body.velocity
        
        # update position of graphics
        body.sphere.pos = body.position
        body.label.pos = body.position

# Runge-Kutta 4 (RK4) Integrator
def Runge_Kutta():
    
    # calculate k1
    for body in bodies:
        body.k = [0] # zero at beginning to push indices one step up for readability
        body.temp_position = body.position
        body.xv = conVec(body.temp_position, body.velocity)
        temp_k1 = gravitational_acc_runge(body.xv)
        k1 = conVec(temp_k1.x * dt, temp_k1.y * dt)
        body.k.append(k1)
        
    # move temp_pos according to k1
    for body in bodies:
        body.temp_position = body.position + body.k[1].x/2
    
    # calculate k2
    for body in bodies:
        temp_xv = conVec(body.xv.x + body.k[1].x/2, body.xv.y + body.k[1].y/2)
        temp_k2 = gravitational_acc_runge(temp_xv)
        k2 = conVec(temp_k2.x * dt, temp_k2.y * dt)
        body.k.append(k2)
    
    # move temp_pos according to k2
    for body in bodies:
        body.temp_position = body.position + body.k[2].x/2

    # calculate k3
    for body in bodies:
        temp_xv = conVec(body.xv.x + body.k[2].x/2, body.xv.y + body.k[2].y/2)
        temp_k3 = gravitational_acc_runge(temp_xv)
        k3 = conVec(temp_k3.x * dt, temp_k3.y * dt)
        body.k.append(k3)
    
    # move temp_pos according to k3
    for body in bodies:
        body.temp_position = body.position + body.k[3].x

    # calculate k4
    for body in bodies:
        temp_xv = conVec(body.xv.x + body.k[3].x, body.xv.y + body.k[3].y)
        temp_k4 = gravitational_acc_runge(temp_xv)
        k4 = conVec(temp_k4.x * dt, temp_k4.y * dt)
        body.k.append(k4)

    # calculate weighted sum and move bodies
    for body in bodies:
        body.position += 1/6 * (body.k[1].x + 2*body.k[2].x + 2*body.k[3].x + body.k[4].x)
        body.velocity += 1/6 * (body.k[1].y + 2*body.k[2].y + 2*body.k[3].y + body.k[4].y)

        body.sphere.pos = body.position
        body.label.pos = body.position

# Velocity-Verlet Integrator
def Verlet():
    # position calculations
    for body in bodies:
        body.position += body.velocity * dt + body.acc/2 * dt**2
        body.sphere.pos = body.position
        body.label.pos = body.position
    # acceleration and velocity calculations
    for body in bodies:
        temp_acc = gravitational_acc(body.position)
        body.velocity += dt/2*(body.acc + temp_acc)
        body.acc = temp_acc

def Forest_Ruth():
    for body in bodies:
        body.position += Theta*dt/2*body.velocity
    for body in bodies:
        body.velocity += Theta*dt*gravitational_acc(body.position)
    for body in bodies:
        body.position += (1-Theta)*dt/2*body.velocity
    for body in bodies:
        body.velocity += (1-2*Theta)*dt*gravitational_acc(body.position)
    for body in bodies:
        body.position += (1-Theta)*dt/2*body.velocity
    for body in bodies:
        body.velocity += Theta*dt*gravitational_acc(body.position)
    for body in bodies:
        body.position += Theta*dt/2*body.velocity
        body.sphere.pos = body.position
        body.label.pos = body.position


def PEFRL():
    for body in bodies:
        body.position += Epsilon*dt*body.velocity
    for body in bodies:
        body.velocity += (1-2*Lambda)*dt/2*gravitational_acc(body.position)
    for body in bodies:
        body.position += Chi*dt*body.velocity
    for body in bodies:
        body.velocity += Lambda*dt*gravitational_acc(body.position)
    for body in bodies:
        body.position += (1-2*(Chi+Epsilon))*dt*body.velocity
    for body in bodies:
        body.velocity += Lambda*dt*gravitational_acc(body.position)
    for body in bodies:
        body.position += Chi*dt*body.velocity
    for body in bodies:
        body.velocity += (1-2*Lambda)*dt/2*gravitational_acc(body.position)
    for body in bodies:
        body.position += Epsilon*dt*body.velocity
        body.sphere.pos = body.position
        body.label.pos = body.position

### End Integrators ###

# parse cmd arguments
args = parser.parse_args()

# load config json file
with open(args.configfile, "r") as configfile:
    config = json.load(configfile)

# use configurations from config json file
if args.useconfig:
    try:
        dt = config[1]["dt"]
    except:
        dt = args.dt
    try:
        scale_factor = config[1]["scale_factor"]
    except:
        scale_factor = args.scale
    try:
        end_time = config[1]["time"]
    except:
        end_time = args.time
    try:
        integrator = config[1]["integrator"]
    except:
        integrator = args.integrator
# use argument configurations
else:
    dt = args.dt
    scale_factor = args.scale
    end_time = args.time
    integrator = args.integrator

# check which integrator was chosen
if integrator.lower() == "euler":
    integrator = Euler
elif integrator.lower() == "rk4":
    integrator = Runge_Kutta
elif integrator.lower() == "verlet":
    integrator = Verlet
elif integrator.lower() == "fr":
    integrator = Forest_Ruth
elif integrator.lower() == "pefrl":
    integrator = PEFRL

# UNITS:
# Mass: solar mass
# Length: Astronomical unit
# Time: days
# G = 4pi^2*AU^3/(M * 365.25) => G = 4*pi^2/365.25^2
#G = 6.67e-11
#scale_factor = 1000
#dt = 0.01

### Constants ###

G = 2.9592e-04
AU = 1.5e11
M = 2e30
Theta = 1/(2-2**(1/3))
Epsilon = 0.1786178958448091
Lambda = -0.2123418310626054
Chi = -0.6626458266981849E-01



# list of all the bodies in the simulation
bodies = []

class Body():
    def __init__(self, mass=1, radius=1, velocity=vector(0,0,0), position=vector(0,0,0), color=color.white, trail=True, name="Body", scale=True, index=0):
        self.mass = mass
        self.velocity = velocity
        self.position = position
        self.temp_position = vector(0,0,0)
        self.k = []
        self.xv = conVec(0,0)
        self.sum_forces = vector(0,0,0)
        self.color = color
        self.radius = radius
        self.acc = gravitational_acc(self.position) # approximate the initial acceleration for Verlet
        self.name = name
        self.label = label(pos=self.position, text=self.name, height=10)
        self.index = index
        if scale:
            self.sphere = sphere(pos=self.position, color=self.color, radius=self.radius*scale_factor, make_trail=trail, retain=200, index=self.index)
        else:
            self.sphere = sphere(pos=self.position, color=self.color, radius=self.radius, make_trail=trail, retain=200, index=self.index)
        #bodies.append(self) # uncomment if you want automatic adding to bodies list


def color_to_vector(color_list):
    return vector(color_list[0]/255, color_list[1]/255, color_list[2]/255)


for body in config[0]:
    bodies.append(Body(
        name = body["name"],
        mass = body["mass"],
        radius = body["radius"],
        position = vector(body["position"][0], body["position"][1], body["position"][2]),
        velocity = vector(body["velocity"][0], body["velocity"][1], body["velocity"][2]),
        trail = body["trail"],
        color = color_to_vector(body["color"]),
        scale = body["scale"],
        index=len(bodies)
    ))

# create all bodies in the simulation

""" sun = Body(mass=1,
           radius=7e8 * 5 / scale_factor / AU,
           color=color.yellow,
           trail=False,
           name="Sun",
           )


earth = Body(mass=6e24/M,
             radius=6.371e6/AU,
             #velocity=vector(0,3e4,0),
             #position=vector(AU,0,0),
             position=vector(-5.111702950987252E-01, -8.734341386147972E-01, 3.902531498407046E-05),
             velocity=vector(1.457401965494037E-02, -8.749957786090569E-03, -3.393201214360642E-07),
             color=color.green,
             name="Earth"
             )


mercury = Body(
               mass=earth.mass*0.055,
               radius=6e6/AU,
               #velocity=vector(0,4.7e4,0),
               #position=vector(0.4*AU,0,0),
               position=vector(3.600062387312980E-01, -8.310671431721671E-02, -3.981766501010686E-02),
               velocity=vector(8.732371820239134E-04, 2.867508157942580E-02, 2.263026727476856E-03),
               color=color.red,
               name="Mercury",
               )


venus = Body(
               mass=earth.mass*0.815,
               radius=6e6/AU,
               #velocity=vector(0,3.5e4,0),
               #position=vector(0.7*AU,0,0),
               position=vector(-5.460148756311848E-01, 4.654289630909307E-01, 3.789319798488837E-02),
               velocity=vector(-1.319751648139675E-02, -1.549708277964608E-02, 5.490020542624818E-04),
               color=color.white,
               name="Venus",
               )


mars = Body(
               mass=earth.mass*0.107,
               radius=6e6/AU,
               velocity=vector(1.444719742599419E-02,-2.365918534978303E-04,-3.594885612448260E-04),
               position=vector(-1.508529480814324E-01,-1.460121856503524,-2.689190873994556E-02),
               color=color.red,
               name="Mars",
               )


jupiter = Body(
               mass=earth.mass*318,
               radius=7e7/AU,
               velocity=vector(5.611682808441865E-03,-4.596785105938998E-03,-1.064356940327842E-04),
               position=vector(-3.545075313382027,-4.081361865858232,9.627457319753692E-02),
               color=color.blue,
               name="Jupiter",
               )


saturn = Body(
               mass=earth.mass*95,
               radius=6e7/AU,
               velocity=vector(5.262021976694793E-03,4.141890616120753E-04,-2.169327374705523E-04),
               position=vector(7.842529344684837E-01,-1.003393486265119E+01,1.431896871358062E-01),
               color=color.white,
               name="Saturn",
               )


uranus = Body(
               mass=earth.mass*14,
               radius=2.5e7/AU,
               velocity=vector(-1.905201349392400E-03,3.265505721711341E-03,3.669073443400500E-05),
               position=vector(1.746114323983198E+01,9.517430938519276,-1.907513002050031E-01),
               color=color.blue,
               name="Uranus",
               )


neptune = Body(
               mass=earth.mass*17,
               radius=2.4e7/AU,
               velocity=vector(8.427417626787077E-04,3.035037625808767E-03,-8.199842541642128E-05),
               position=vector(2.880079206580985E+01,-8.173900363488711,-4.954784189728160E-01),
               color=color.yellow,
               name="Neptune",
               )


pluto = Body(
               mass=earth.mass*0.0022,
               radius=1e6/AU,
               #velocity=vector(0,4.67e3,0),
               #position=vector(39*AU,0,0),
               position=vector(1.120198708794019E+01, -3.164123744663468E+01, 1.446313453325374E-01),
               velocity=vector(3.029567845289497E-03, 3.743167934314588E-04, -9.260696937062970E-04),
               color=color.yellow,
               name="Pluto",
               )
 """
# black_hole = Body(
#                   mass=4e6,
#                   radius=2.5e7/AU,
#                   color=color.white,
#                   position=vector(15, -20, 0),
#                   velocity=vector(0, 0.5, 0),
# )

# hole = Body(#mass=2e30*4.3e6,
           # mass=2e30*10,
           # radius=7e8 * 20,
           # color=color.yellow,
           # position=vector(5.2*AU,-39*AU,0),
           # velocity=vector(0,100e3,0),
           # trail=False,
           # )

current_time = 0

time_label = label(pos=vector(20, 350, 0), pixel_pos=True, align='left', text="Time: " + str(current_time/365) + " years")

### info box ###
info_label = label(pos=vector(20, scene.height/3, 0), pixel_pos=True, box=False, align='left', text="")

def onClick(e):
    obj = scene.mouse.pick  # get the sphere
    if(obj != None):
        body = bodies[obj.index]  # each sphere has a index attribute wich points to the planets position in the bodies list
        # info about the planets can noe be retrieved from the class
        info_label.text = '<b><i>'+body.name+'</i></b>\n<i>Mass:</i> '+str(body.mass)+' M☉\n'
        # TODO add more info about planets & convert units
    else:
        info_label.text = ''


scene.bind('click', onClick)
### ###

# loop over every body and run its update method every timestep

if end_time > 0:
    for epoch in range(int(end_time/dt)):
        rate(args.rate)

        integrator()

        current_time = epoch*dt
        time_label.text = "Time: {:.2f} years".format(current_time/365)
    # print body positions for benchmarking
    if args.printEndPos:
        for body in bodies:
            print(f"{body.name}: {body.position}")
    if args.checkEndPos:
        error_sum = 0
        for body in bodies:
            end_pos = config[0][body.index]["end_position"]
            end_pos = vector(end_pos[0], end_pos[1], end_pos[2])
            error = mag(body.position - end_pos) # the magnitude of the error
            error_sum += error
            print(f"{body.name}: {error} AU")
        print(f"Total error: {error_sum}")
        print(f"dt: {dt}")
        print(f"Integrator: {args.integrator}")

        

else:
    while True:
        rate(args.rate)
        
        integrator()
        
        time_label.text = "Time: {:.2f} years".format(current_time/365)
        current_time += dt
