def run():
    import argparse
    import time
    from math import pi
    from headless_init import parameters
    from tools import vector

    # argument parsing
    parser = argparse.ArgumentParser(description="A newtonian gravity simulator")
    parser.add_argument(
        "-t",
        "--time",
        type=float,
        default=0,
        dest="time",
        help="The amount of time the simulation will simulate measured in days. Type 0 for infinite time. (Default: 0)",
    )
    parser.add_argument(
        "--dt",
        type=float,
        default=0.1,
        dest="dt",
        help="The timestep to use. (Default: 0.1)",
    )
    parser.add_argument(
        "--scale",
        type=float,
        default=1000,
        dest="scale",
        help="The number to scale the radiuses of the planets to make them visible. Does only affect the visuals not collisions. (Default: 1000)",
    )
    parser.add_argument(
        "--rate",
        type=int,
        default=100000,
        dest="rate",
        help="Number of timesteps per second (Default: 100000)",
    )
    parser.add_argument(
        "--configfile",
        type=str,
        default="config.json",
        dest="configfile",
        help="Path to the config file containing the bodies. (Default: config.json)",
    )
    parser.add_argument(
        "--useconfig",
        action="store_true",
        default=False,
        dest="useconfig",
        help="Use this flag if you want to use the settings in the configfile instead of defaults and cmd arguments. (Default: False)",
    )
    parser.add_argument(
        "--integrator",
        type=str,
        default="verlet",
        dest="integrator",
        help="The integrator to be used. Options: euler, verlet, rk4, fr, pefrl (Default: verlet)",
    )
    parser.add_argument(
        "--endPos",
        action="store_true",
        default=False,
        dest="printEndPos",
        help="When flagged the end position of all bodies will be printed (Deafult: False)",
    )
    parser.add_argument(
        "--checkEndPos",
        action="store_true",
        default=False,
        dest="checkEndPos",
        help="When flagged the end position of all bodies is compared to their real end positions, which are given as 'end_position' in config.json (Default: False)",
    )

    # parse cmd arguments
    args = parser.parse_args()

    p = parameters(
        args=args,
        G=4 * pi ** 2 / 365.2422 ** 2,
        AU=1.5e11,
        M=2e30,
        Theta=1 / (2 - 2 ** (1 / 3)),
        Epsilon=0.1786178958448091,
        Lambda=-0.2123418310626054,
        Chi=-0.6626458266981849E-01,
    )

    # loop over every body and run its update method every timestep
    p.start_time = time.time()
    if p.end_time != 0:
        if p.adaptive:
            while p.current_time < p.end_time:
                p.integrator(p)

                p.old_times[0] = p.old_times[1]
                p.old_times[1] = p.old_times[2]
                p.old_times[2] = p.current_time

                p.current_time += p.dt
                print(p.current_time)
                p.dt = p.next_dt

            # do interpolation here
            for body in p.all_bodies:
                x = p.end_time
                x1 = p.old_times[0]
                x2 = p.old_times[1]
                x3 = p.old_times[2]
                x4 = p.current_time
                y1 = body.old_positions[0]
                y2 = body.old_positions[1]
                y3 = body.old_positions[2]
                y4 = body.position
                # Lagrange Interpolation Polynomial
                body.position = (
                    (x - x2)
                    * (x - x3)
                    * (x - x4)
                    / ((x1 - x2) * (x1 - x3) * (x1 - x4))
                    * y1
                    + (x - x1)
                    * (x - x3)
                    * (x - x4)
                    / ((x2 - x1) * (x2 - x3) * (x2 - x4))
                    * y2
                    + (x - x1)
                    * (x - x2)
                    * (x - x4)
                    / ((x3 - x1) * (x3 - x2) * (x3 - x4))
                    * y3
                    + (x - x1)
                    * (x - x2)
                    * (x - x3)
                    / ((x4 - x1) * (x4 - x2) * (x4 - x3))
                    * y4
                )

        else:
            for epoch in range(int(p.end_time / p.dt)):
                p.integrator(p)
                p.current_time = epoch * p.dt
        # print body positions for benchmarking
        if p.printEndPos:
            for body in p.bodies:
                print(f"{body.name}: {body.position}")
        if p.checkEndPos:
            error_sum = 0
            n_error_counter = 0
            for body in p.bodies:
                try:
                    end_pos = p.config[0][body.index]["end_position"]
                    end_pos = vector(end_pos[0], end_pos[1], end_pos[2])
                    error = (body.position - end_pos).mag  # the magnitude of the error
                    error_sum += error
                    n_error_counter += 1
                    print(f"{body.name}: {error} AU")
                except:
                    print(f"{body.name} has no endPos")
            print(f"Total error: {error_sum}")
            print(f"Average error: {error_sum/n_error_counter}")
            print(f"dt: {p.dt}")
            print(f"Integrator: {p.args.integrator}")

    else:
        while True:
            p.integrator(p)
            p.current_time += p.dt

    print(f"Execution time: {time.time()-p.start_time} seconds")

    # UNITS:
    # Mass: solar mass
    # Length: Astronomical unit
    # Time: days
    # G = 4pi^2*AU^3/(M * 365.25) => G = 4*pi^2/365.25^2
    # G = 6.67e-11
    # AU^3/D^2 = 1/448485856027460.06  km^3/s^2 = 2.2297247205467538e-15 * km^3/s^2
    # scale_factor = 1000
    # dt = 0.01

    ### Constants ###

    # G = 0.00029592338593516714
    # G = 4*pi**2/365.25**2
    # AU = 1.5e11
    # M = 2e30
    # Theta = 1/(2-2**(1/3))
    # Epsilon = 0.1786178958448091
    # Lambda = -0.2123418310626054
    # Chi = -0.6626458266981849E-01


if __name__ == "__main__":
    run()
