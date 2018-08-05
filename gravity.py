from vpython import *
import argparse
import json

# argument parsing
parser = argparse.ArgumentParser(description="A newtonian gravity simulator")
parser.add_argument("-t", "--time", type=float, default=0, dest="time",
                    help="The amount of time the simulation will simulate measured in days. Type 0 for infinite time. (Default: 0)")
parser.add_argument("--dt", type=float, default=0.01, dest="dt",
                    help="The timestep to use. (Default: 0.01)")
parser.add_argument("--scale", type=float, default=1000, dest="scale",
                    help="The number to scale the radiuses of the planets to make them visible. Does only affect the visuals not collisions. (Default: 1000)")
parser.add_argument("--rate", type=int, default=10000, dest="rate",
                    help="Number of timesteps per second (Default: 10000)")
parser.add_argument("--configfile", type=str, default="config.json", dest="configfile",
                    help="Path to the config file containing the bodies. (Default: config.json)")
parser.add_argument("--useconfig", action="store_true", default=False, dest="useconfig",
                    help="Use this flag if you want to use the settings in the configfile instead of defaults and cmd arguments. (Default: False)")
# TODO Start
parser.add_argument("--integrator", type=str, default="euler", dest="integrator",
                    help="The integrator to be used. Options: euler, verlet, runge (Default: euler)")

# TODO End

args = parser.parse_args()

with open(args.configfile, "r") as configfile:
    config = json.load(configfile)

if args.useconfig:
    dt = config[1]["dt"]
    scale_factor = config[1]["scale_factor"]
    end_time = config[1]["time"]
else:
    dt = args.dt
    scale_factor = args.scale
    end_time = args.time


# UNITS:
# Mass: solar mass
# Length: Astronomical unit
# Time: days
# G = 4pi^2*AU^3/(M * 365.25) => G = 4*pi^2/365.25^2
#G = 6.67e-11
#scale_factor = 1000
#dt = 0.01

G = 2.9592e-04
time = 0

AU = 1.5e11
M = 2e30

# list of all the bodies in the simulation
bodies = []

class Body():
    def __init__(self, mass=1, radius=1, velocity=vector(0,0,0), position=vector(0,0,0), color=color.white, trail=True, name="Body", scale=True, index=0):
        self.mass = mass
        self.velocity = velocity
        self.position = position
        self.color = color
        self.radius = radius
        self.forces = []
        self.acc = vector(0,0,0)
        self.sum_force = vector(0,0,0)
        self.name = name
        self.label = label(pos=self.position, text=self.name, height=10)
        self.index = index
        if scale:
            self.sphere = sphere(pos=self.position, color=self.color, radius=self.radius*scale_factor, make_trail=trail, retain=200, index=self.index)
        else:
            self.sphere = sphere(pos=self.position, color=self.color, radius=self.radius, make_trail=trail, retain=200, index=self.index)
        #bodies.append(self) # uncomment if you want automatic adding to bodies list
    def update(self):
        self.forces = []
        self.sum_force = vector(0,0,0)
        self.gravitational_force()
        # add other forces here

        # sum all forces (Newtons Second Law)
        for force in self.forces:
            self.sum_force += force

        # Euler method
        self.acc = self.sum_force/self.mass
        self.velocity += dt * self.acc

    def move(self):
        self.position += dt * self.velocity

        # update position of graphics
        self.sphere.pos = self.position
        self.label.pos = self.position

    def gravitational_force(self):
        # calculate the gravitational force from all other bodies
        for body in bodies:
            # distance to the other body
            r = mag(self.position-body.position)
            # skip if body is itself
            if r < self.radius+body.radius:
                continue
            # the magnitude of the force
            force = G * self.mass * body.mass / r**2
            # the unit vector for the force
            dir = norm(body.position - self.position)
            # the force vector
            force = force * dir
            # add force vector to list of forces
            self.forces.append(force)

    # this function is an alternative to the Euler integration
    def updateVerlet(self):
        self.forces = []
        self.sum_force = vector(0,0,0)
        self.gravitational_force()
        # add other forces here

        # sum the forces
        for force in self.forces:
            self.sum_force += force


        self.position += self.velocity * dt + self.acc/2 * dt**2
        self.velocity += dt/2*(self.acc + self.sum_force/self.mass)
        self.acc = self.sum_force/self.mass

        self.sphere.pos = self.position
        self.label.pos = self.position


def color_to_vector(color_list):
    return vector(color_list[0]/255, color_list[1]/255, color_list[2]/255)


for body in config[0]:
    print(body["name"])
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

time_label = label(pos=vector(20, 350, 0), pixel_pos=True, align='left', text="Time: " + str(time/365) + " years")

### info box ###
info_label = label(pos=vector(20, scene.height/3, 0), pixel_pos=True, box=False, align='left', text="")

def onClick(e):
    obj = scene.mouse.pick  # get the sphere
    if(obj != None):
        body = bodies[obj.index]  # each sphere has a index attribute wich points to the planets position in the bodies list
        # info about the planets can noe be retrieved from the class
        info_label.text = '<b><i>'+body.name+'</i></b>\n<i>Mass:</i> '+str(body.mass)+' M☉\n'
        # TODO add more info about planets & convert units


scene.bind('click', onClick)
### ###

# loop over every body and run its update method every timestep
if end_time > 0:
    while time < end_time:
        rate(args.rate)
        # calculate the change in velocity for all bodies...
        for body in bodies:
            body.update()
        # ... then move all bodies
        for body in bodies:
            body.move()
        time_label.text = "Time: {:.2f} years".format(time/365)
        time += dt

else:
    while True:
        rate(args.rate)
        # calculate the change in velocity for all bodies...
        for body in bodies:
            body.update()
        # ... then move all bodies
        for body in bodies:
            body.move()
        time_label.text = "Time: {:.2f} years".format(time/365)
        time += dt
