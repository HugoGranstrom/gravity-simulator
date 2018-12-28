from tools import conVec, vector

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


# Runge-Kutta 4 (RK4) Integrator
def Runge_Kutta(p):

    # calculate k1
    for body in p.all_bodies:
        body.k = [0]  # zero at beginning to push indices one step up for readability
        body.temp_position = body.position
        body.xv = conVec(body.temp_position, body.velocity)
        temp_k1 = gravitational_acc_runge(body.xv, p)
        k1 = conVec(temp_k1.x * p.dt, temp_k1.y * p.dt)
        body.k.append(k1)

    # move temp_pos according to k1
    for body in p.all_bodies:
        body.temp_position = body.position + body.k[1].x / 2

    # calculate k2
    for body in p.all_bodies:
        temp_xv = conVec(body.xv.x + body.k[1].x / 2, body.xv.y + body.k[1].y / 2)
        temp_k2 = gravitational_acc_runge(temp_xv, p)
        k2 = conVec(temp_k2.x * p.dt, temp_k2.y * p.dt)
        body.k.append(k2)

    # move temp_pos according to k2
    for body in p.all_bodies:
        body.temp_position = body.position + body.k[2].x / 2

    # calculate k3
    for body in p.all_bodies:
        temp_xv = conVec(body.xv.x + body.k[2].x / 2, body.xv.y + body.k[2].y / 2)
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
        body.position += (
            1 / 6 * (body.k[1].x + 2 * body.k[2].x + 2 * body.k[3].x + body.k[4].x)
        )
        body.velocity += (
            1 / 6 * (body.k[1].y + 2 * body.k[2].y + 2 * body.k[3].y + body.k[4].y)
        )


def Runge_Kutta_45(p):
    tol = 1e-4
    dt_min = 1e-2
    dt_max = 1
    while True:
        # calculate k1
        for body in p.all_bodies:
            body.k = [
                0
            ]  # zero at beginning to push indices one step up for readability
            body.temp_position = body.position
            body.xv = conVec(body.temp_position, body.velocity)
            temp_k1 = gravitational_acc_runge(body.xv, p)
            k1 = conVec(temp_k1.x * p.dt, temp_k1.y * p.dt)
            body.k.append(k1)
        # move temp_pos according to k1
        for body in p.all_bodies:
            body.temp_position = body.position + body.k[1].x / 4

        # calculate k2
        for body in p.all_bodies:
            temp_xv = conVec(body.xv.x + body.k[1].x / 4, body.xv.y + body.k[1].y / 4)
            temp_k2 = gravitational_acc_runge(temp_xv, p)
            k2 = conVec(temp_k2.x * p.dt, temp_k2.y * p.dt)
            body.k.append(k2)

        # move temp_pos according to k2
        for body in p.all_bodies:
            body.temp_position = (
                body.position + 3 / 32 * body.k[1].x + 9 / 32 * body.k[2].x
            )

        # calculate k3
        for body in p.all_bodies:
            temp_xv = conVec(
                body.xv.x + 3 / 32 * body.k[1].x + 9 / 32 * body.k[2].x,
                body.xv.y + 3 / 32 * body.k[1].y + 9 / 32 * body.k[2].y,
            )
            temp_k3 = gravitational_acc_runge(temp_xv, p)
            k3 = conVec(temp_k3.x * p.dt, temp_k3.y * p.dt)
            body.k.append(k3)

        # move temp_pos according to k3
        for body in p.all_bodies:
            body.temp_position = (
                body.position
                + 1932 / 2197 * body.k[1].x
                - 7200 / 2197 * body.k[2].x
                + 7296 / 2197 * body.k[3].x
            )

        # calculate k4
        for body in p.all_bodies:
            temp_xv = conVec(
                body.xv.x
                + 1932 / 2197 * body.k[1].x
                - 7200 / 2197 * body.k[2].x
                + 7296 / 2197 * body.k[3].x,
                body.xv.y
                + 1932 / 2197 * body.k[1].y
                - 7200 / 2197 * body.k[2].y
                + 7296 / 2197 * body.k[3].y,
            )
            temp_k4 = gravitational_acc_runge(temp_xv, p)
            k4 = conVec(temp_k4.x * p.dt, temp_k4.y * p.dt)
            body.k.append(k4)

        # move temp_pos according to k4
        for body in p.all_bodies:
            body.temp_position = (
                body.position
                + 439 / 216 * body.k[1].x
                - 8 * body.k[2].x
                + 3680 / 513 * body.k[3].x
                - 845 / 4104 * body.k[4].x
            )

        # calculate k5
        for body in p.all_bodies:
            temp_xv = conVec(
                body.xv.x
                + 439 / 216 * body.k[1].x
                - 8 * body.k[2].x
                + 3680 / 513 * body.k[3].x
                - 845 / 4104 * body.k[4].x,
                body.xv.y
                + 439 / 216 * body.k[1].y
                - 8 * body.k[2].y
                + 3680 / 513 * body.k[3].y
                - 845 / 4104 * body.k[4].y,
            )
            temp_k5 = gravitational_acc_runge(temp_xv, p)
            k5 = conVec(temp_k5.x * p.dt, temp_k5.y * p.dt)
            body.k.append(k5)

        # move temp_pos according to k5
        for body in p.all_bodies:
            body.temp_position = (
                body.position
                - 8 / 27 * body.k[1].x
                + 2 * body.k[2].x
                - 3544 / 2565 * body.k[3].x
                + 1859 / 4104 * body.k[4].x
                - 11 / 40 * body.k[5].x
            )

        # calculate k6
        for body in p.all_bodies:
            temp_xv = conVec(
                body.xv.x
                - 8 / 27 * body.k[1].x
                + 2 * body.k[2].x
                - 3544 / 2565 * body.k[3].x
                + 1859 / 4104 * body.k[4].x
                - 11 / 40 * body.k[5].x,
                body.xv.y
                - 8 / 27 * body.k[1].y
                + 2 * body.k[2].y
                - 3544 / 2565 * body.k[3].y
                + 1859 / 4104 * body.k[4].y
                - 11 / 40 * body.k[5].y,
            )
            temp_k6 = gravitational_acc_runge(temp_xv, p)
            k6 = conVec(temp_k6.x * p.dt, temp_k6.y * p.dt)
            body.k.append(k6)

        total_error = 0
        for body in p.all_bodies:
            body.temp_position = (
                body.position
                + 25 / 216 * body.k[1].x
                + 1408 / 2565 * body.k[3].x
                + 2197 / 4101 * body.k[4].x
                - 1 / 5 * body.k[5].x
            )
            z = (
                body.position
                + 16 / 135 * body.k[1].x
                + 6656 / 12825 * body.k[3].x
                + 28561 / 56430 * body.k[4].x
                - 9 / 50 * body.k[5].x
                + 2 / 55 * body.k[6].x
            )
            error = (body.temp_position - (z)).mag
            total_error += error
        p.next_dt = p.dt * 0.84 * (tol * p.dt / total_error) ** 0.25
        if p.next_dt < dt_min:
            p.next_dt = dt_min
        elif p.next_dt > dt_max:
            p.next_dt = dt_max
        if total_error < tol:
            break
        else:
            p.dt = p.next_dt

    # calculate weighted sum and move bodies
    for body in p.all_bodies:
        body.old_positions[0] = body.old_positions[1]
        body.old_positions[1] = body.old_positions[2]
        body.old_positions[2] = body.position

        body.position = body.temp_position
        body.velocity += (
            25 / 216 * body.k[1].y
            + 1408 / 2565 * body.k[3].y
            + 2197 / 4101 * body.k[4].y
            - 1 / 5 * body.k[5].y
        )


# Velocity-Verlet Integrator
def Verlet(p):
    # position calculations
    for body in p.all_bodies:
        body.position += body.velocity * p.dt + body.acc / 2 * p.dt ** 2
    # acceleration and velocity calculations
    for body in p.all_bodies:
        temp_acc = gravitational_acc(body.position, p)
        body.velocity += p.dt / 2 * (body.acc + temp_acc)
        body.acc = temp_acc


def Forest_Ruth(p):
    for body in p.all_bodies:
        body.position += p.Theta * p.dt / 2 * body.velocity
    for body in p.all_bodies:
        body.velocity += p.Theta * p.dt * gravitational_acc(body.position, p)
    for body in p.all_bodies:
        body.position += (1 - p.Theta) * p.dt / 2 * body.velocity
    for body in p.all_bodies:
        body.velocity += (1 - 2 * p.Theta) * p.dt * gravitational_acc(body.position, p)
    for body in p.all_bodies:
        body.position += (1 - p.Theta) * p.dt / 2 * body.velocity
    for body in p.all_bodies:
        body.velocity += p.Theta * p.dt * gravitational_acc(body.position, p)
    for body in p.all_bodies:
        body.position += p.Theta * p.dt / 2 * body.velocity


def PEFRL(p):
    for body in p.all_bodies:
        body.position += p.Epsilon * p.dt * body.velocity
    for body in p.all_bodies:
        body.velocity += (
            (1 - 2 * p.Lambda) * p.dt / 2 * gravitational_acc(body.position, p)
        )
    for body in p.all_bodies:
        body.position += p.Chi * p.dt * body.velocity
    for body in p.all_bodies:
        body.velocity += p.Lambda * p.dt * gravitational_acc(body.position, p)
    for body in p.all_bodies:
        body.position += (1 - 2 * (p.Chi + p.Epsilon)) * p.dt * body.velocity
    for body in p.all_bodies:
        body.velocity += p.Lambda * p.dt * gravitational_acc(body.position, p)
    for body in p.all_bodies:
        body.position += p.Chi * p.dt * body.velocity
    for body in p.all_bodies:
        body.velocity += (
            (1 - 2 * p.Lambda) * p.dt / 2 * gravitational_acc(body.position, p)
        )
    for body in p.all_bodies:
        body.position += p.Epsilon * p.dt * body.velocity


### End Integrators ###

### Functions ###
# gravitational acceleration for Euler and Verlet
def gravitational_acc(position, p):
    sum_acc = vector(0, 0, 0)
    # calculate the gravitational acceleration from all other bodies
    for body in p.bodies:
        # distance to the other body
        r_vec = body.position - position
        r = r_vec.mag
        # skip if body is itself
        if r == 0:
            continue
        acc = r_vec * body.GM / r ** 3
        # add force vector to the sum of forces
        sum_acc += acc
    return sum_acc


# gravitational acceleration for Runge-Kutta
def gravitational_acc_runge(xv, p):
    sum_acc = vector(0, 0, 0)
    # calculate the gravitational acceleration from all other bodies
    for body in p.bodies:
        # distance to the other body
        r_vec = body.temp_position - xv.x
        r = r_vec.mag
        # skip if body is itself
        if r < body.radius:
            continue
        # # the magnitude of the force
        # acc = body.GM / r_vec.mag2
        # # the unit vector for the force
        # dir = r_vec.hat
        # # the force vector
        # acc = acc * dir

        acc = r_vec * body.GM / r ** 3

        # add force vector to the sum of forces
        sum_acc += acc
    return conVec(xv.y, sum_acc)
