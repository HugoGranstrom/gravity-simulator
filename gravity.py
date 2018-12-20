from vpython import *
import argparse
import time
from math import pi
import integrators
from gravity_init import parameters
def run():


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

    time_label = label(pos=vector(20, 350, 0), pixel_pos=True, align='left', text="Time: " + str(p.current_time/365) + " years")

    ### info box ###
    info_label = label(pos=vector(20, scene.height/3, 0), pixel_pos=True, box=False, align='left', text="")
    scene.bind('click', onClick)
    ### ###

    # loop over every body and run its update method every timestep
    p.start_time = time.time()
    if p.end_time != 0:
        for epoch in range(int(p.end_time/p.dt)):
            rate(p.args.rate)

            p.integrator(p)
            for body in p.all_bodies:
                body.sphere.pos = body.position
                body.label.pos = body.position
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
            for body in p.all_bodies:
                body.sphere.pos = body.position
                body.label.pos = body.position
            time_label.text = "Time: {:.2f} years".format(p.current_time/365)
            p.current_time += p.dt

    print(f"Execution time: {time.time()-p.start_time} seconds")

if __name__ == "__main__":
    run()
