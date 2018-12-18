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

# Velocity-Verlet Integrator
def Verlet(p):
    # position calculations
    for body in p.all_bodies:
        body.position += body.velocity * p.dt + body.acc/2 * p.dt**2
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

### End Integrators ###

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
        acc = body.GM / r_vec.mag2
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
        acc = body.GM / r_vec.mag2
        # the unit vector for the force
        dir = r_vec.hat
        # the force vector
        acc = acc * dir
        # add force vector to the sum of forces
        sum_acc += acc
    return conVec(xv.y, sum_acc)