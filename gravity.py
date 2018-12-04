from vpython import *
import argparse
import json
from collections import namedtuple
import time
from math import pi;

### Functions ###
# gravitational acceleration for Euler and Verlet
def gravitational_acc(position):
    sum_acc = vector(0,0,0)
    # calculate the gravitational acceleration from all other bodies
    for body in bodies:
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
def gravitational_acc_runge(xv):
    sum_acc = vector(0,0,0)
    # calculate the gravitational acceleration from all other bodies
    for body in bodies:
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
def Euler():
    # acceleration and velocity calculations
    for body in all_bodies:
        body.acc = gravitational_acc(body.position)
        body.velocity += dt * body.acc
    # position calculations
    for body in all_bodies:
        body.position += dt * body.velocity
        
        # update position of graphics
        body.sphere.pos = body.position
        body.label.pos = body.position
        
# Runge-Kutta 4 (RK4) Integrator
def Runge_Kutta():
    
    # calculate k1
    for body in all_bodies:
        body.k = [0] # zero at beginning to push indices one step up for readability
        body.temp_position = body.position
        body.xv = conVec(body.temp_position, body.velocity)
        temp_k1 = gravitational_acc_runge(body.xv)
        k1 = conVec(temp_k1.x * dt, temp_k1.y * dt)
        body.k.append(k1)
        
    # move temp_pos according to k1
    for body in all_bodies:
        body.temp_position = body.position + body.k[1].x/2
    
    # calculate k2
    for body in all_bodies:
        temp_xv = conVec(body.xv.x + body.k[1].x/2, body.xv.y + body.k[1].y/2)
        temp_k2 = gravitational_acc_runge(temp_xv)
        k2 = conVec(temp_k2.x * dt, temp_k2.y * dt)
        body.k.append(k2)
    
    # move temp_pos according to k2
    for body in all_bodies:
        body.temp_position = body.position + body.k[2].x/2

    # calculate k3
    for body in all_bodies:
        temp_xv = conVec(body.xv.x + body.k[2].x/2, body.xv.y + body.k[2].y/2)
        temp_k3 = gravitational_acc_runge(temp_xv)
        k3 = conVec(temp_k3.x * dt, temp_k3.y * dt)
        body.k.append(k3)
    
    # move temp_pos according to k3
    for body in all_bodies:
        body.temp_position = body.position + body.k[3].x

    # calculate k4
    for body in all_bodies:
        temp_xv = conVec(body.xv.x + body.k[3].x, body.xv.y + body.k[3].y)
        temp_k4 = gravitational_acc_runge(temp_xv)
        k4 = conVec(temp_k4.x * dt, temp_k4.y * dt)
        body.k.append(k4)

    # calculate weighted sum and move bodies
    for body in all_bodies:
        body.position += 1/6 * (body.k[1].x + 2*body.k[2].x + 2*body.k[3].x + body.k[4].x)
        body.velocity += 1/6 * (body.k[1].y + 2*body.k[2].y + 2*body.k[3].y + body.k[4].y)

        body.sphere.pos = body.position
        body.label.pos = body.position

# Velocity-Verlet Integrator
def Verlet():
    # position calculations
    for body in all_bodies:
        body.position += body.velocity * dt + body.acc/2 * dt**2
        body.sphere.pos = body.position
        body.label.pos = body.position
    # acceleration and velocity calculations
    for body in all_bodies:
        temp_acc = gravitational_acc(body.position)
        body.velocity += dt/2*(body.acc + temp_acc)
        body.acc = temp_acc

def Forest_Ruth():
    for body in all_bodies:
        body.position += Theta*dt/2*body.velocity
    for body in all_bodies:
        body.velocity += Theta*dt*gravitational_acc(body.position)
    for body in all_bodies:
        body.position += (1-Theta)*dt/2*body.velocity
    for body in all_bodies:
        body.velocity += (1-2*Theta)*dt*gravitational_acc(body.position)
    for body in all_bodies:
        body.position += (1-Theta)*dt/2*body.velocity
    for body in all_bodies:
        body.velocity += Theta*dt*gravitational_acc(body.position)
    for body in all_bodies:
        body.position += Theta*dt/2*body.velocity
        body.sphere.pos = body.position
        body.label.pos = body.position


def PEFRL():
    for body in all_bodies:
        body.position += Epsilon*dt*body.velocity
    for body in all_bodies:
        body.velocity += (1-2*Lambda)*dt/2*gravitational_acc(body.position)
    for body in all_bodies:
        body.position += Chi*dt*body.velocity
    for body in all_bodies:
        body.velocity += Lambda*dt*gravitational_acc(body.position)
    for body in all_bodies:
        body.position += (1-2*(Chi+Epsilon))*dt*body.velocity
    for body in all_bodies:
        body.velocity += Lambda*dt*gravitational_acc(body.position)
    for body in all_bodies:
        body.position += Chi*dt*body.velocity
    for body in all_bodies:
        body.velocity += (1-2*Lambda)*dt/2*gravitational_acc(body.position)
    for body in all_bodies:
        body.position += Epsilon*dt*body.velocity
        body.sphere.pos = body.position
        body.label.pos = body.position

### End Integrators ###



class Body():
    def __init__(self, mass=1, GM=1, radius=1, velocity=vector(0,0,0), position=vector(0,0,0), color=color.white, trail=True, name="Body", scale=True, index=0):
        self.mass = mass
        self.GM = GM
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

def onClick(e):
        obj = scene.mouse.pick  # get the sphere
        if(obj != None):
            body = bodies[obj.index]  # each sphere has a index attribute wich points to the planets position in the bodies list
            # info about the planets can noe be retrieved from the class
            info_label.text = '<b><i>'+body.name+'</i></b>\n<i>Mass:</i> '+str(body.mass)+' M☉\n' + '<i>Position:</i> '+str(body.position)+'\n'
            # TODO add more info about planets & convert units
            scene.camera.follow(obj)
        else:
            info_label.text = ''
            #scene.camera.follow()

def run():
    # define globals
    global args
    global G
    global AU
    global M
    global Theta
    global Epsilon
    global Lambda
    global Chi
    global conVec
    global dt
    global scale_factor
    global end_time
    global integrator
    global bodies
    global comets
    global all_bodies
    global info_label

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

    # parse cmd arguments
    args = parser.parse_args()


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
    G = 4*pi**2/365.25**2
    AU = 1.5e11
    M = 2e30
    Theta = 1/(2-2**(1/3))
    Epsilon = 0.1786178958448091
    Lambda = -0.2123418310626054
    Chi = -0.6626458266981849E-01

    ### containerVector ###
    conVec = namedtuple("conVec","x y")



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


    current_time = 0

    # list of all the bodies in the simulation
    bodies = []
    comets = []
    all_bodies = []

    for body in config[0]:
        try:
            GM = body["gm"]
        except:
            GM = body["mass"] * G
        if not body["comet"]:
            bodies.append(Body(
                name = body["name"],
                mass = body["mass"],
                GM = GM,
                radius = body["radius"],
                position = vector(body["position"][0], body["position"][1], body["position"][2]),
                velocity = vector(body["velocity"][0], body["velocity"][1], body["velocity"][2]),
                trail = body["trail"],
                color = color_to_vector(body["color"]),
                scale = body["scale"],
                index=len(bodies)
            ))
        else:
            comets.append(Body(
                name = body["name"],
                mass = body["mass"],
                GM = GM,
                radius = body["radius"],
                position = vector(body["position"][0], body["position"][1], body["position"][2]),
                velocity = vector(body["velocity"][0], body["velocity"][1], body["velocity"][2]),
                trail = body["trail"],
                color = color_to_vector(body["color"]),
                scale = body["scale"],
                index=len(bodies)
            ))

    all_bodies = bodies + comets



    time_label = label(pos=vector(20, 350, 0), pixel_pos=True, align='left', text="Time: " + str(current_time/365) + " years")

    ### info box ###
    info_label = label(pos=vector(20, scene.height/3, 0), pixel_pos=True, box=False, align='left', text="")
    scene.bind('click', onClick)
    ### ###

    # loop over every body and run its update method every timestep
    start_time = time.time()
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

    print(f"Execution time: {time.time()-start_time} seconds")

if __name__ == "__main__":
    run()