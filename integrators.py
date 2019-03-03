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
    dt_min = 1e-3
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
                + 2197 / 4104 * body.k[4].x
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
            + 2197 / 4104 * body.k[4].y
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

def Nystrom2(p):
    for body in p.all_bodies:
        body.k = [0]
        body.k.append(p.dt * gravitational_acc(body.position, p))
    
    for body in p.all_bodies:
        body.temp_position = body.position + p.dt * 2/3 * body.velocity + 2/9 * p.dt * body.k[1]

    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))

    for body in p.all_bodies:
        body.position = body.position + p.dt * body.velocity + p.dt * (body.k[1]/4 + body.k[2]/4)
        body.velocity = body.velocity + body.k[1]/4 + body.k[2] * 3/4

def Nystrom3(p):
    """
    delta2 = 0.6 - 0.06**0.5
    delta3 = 0.6 + 0.06**0.5
    a1 = 0.21 - 0.6 * 0.06*0.5
    b1 = (0.15 + 4 * 0.06**0.5)/25
    b2 = (5.1 + 11 * 0.06**0.5)/25
    alpha1 = 1/9
    alpha2 = (7 + 20 * 0.06**0.5) / 36
    alpha3 = (7 - 20 * 0.06**0.5) / 36
    beta1 = 1/9
    beta2 = (8 + 5 * 0.06**0.5)/18
    beta3 = (8 - 5 * 0.06**0.5)/18
    """
    delta2 = 0.3550510257
    delta3 = 0.8449489743 
    a1 = 0.0630306154 
    b1 = 0.0451918359
    b2 = 0.3117775487
    alpha1 = 0.1111111111
    alpha2 = 0.3305272081
    alpha3 = 0.0583616809
    beta1 = 0.1111111111
    beta2 = 0.5124858262
    beta3 = 0.3764030627

    for body in p.all_bodies:
        body.k = [0]
        body.k.append(p.dt * gravitational_acc(body.position, p))
    
    for body in p.all_bodies:
        body.temp_position = body.position + delta2 * p.dt * body.velocity + a1 * p.dt * body.k[1]
    
    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
        
    
    for body in p.all_bodies:
        body.temp_position = body.position + delta3 * p.dt * body.velocity + p.dt * (b1*body.k[1] + b2 * body.k[2])

    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))

    for body in p.all_bodies:
        body.position = body.position + p.dt * body.velocity + p.dt * (alpha1 * body.k[1] + alpha2 * body.k[2] + alpha3 * body.k[3])
        body.velocity = body.velocity + beta1 * body.k[1] + beta2 * body.k[2] + beta3 * body.k[3]

def Nystrom4(p):
    delta2 = 0.2123405385
    delta3 = 0.5905331358
    delta4 = 0.9114120406
    a1 = 0.02254425214 
    b1 = -0.0011439805
    b2 = 0.1755086728
    c1 = 0.1171541673
    c2 = 0.1393754710
    c3 = 0.1588063156
    alpha1 = 0.0625000001
    alpha2 = 0.2590173402
    alpha3 = 0.1589523623
    alpha4 = 0.0195302974
    beta1 = 0.0625000001
    beta2 = 0.3288443202
    beta3 = 0.3881934687
    beta4 = 0.2204622110

    for body in p.all_bodies:
        body.k = [0]
        body.k.append(p.dt * gravitational_acc(body.position, p))
    
    for body in p.all_bodies:
        body.temp_position = body.position + delta2 * p.dt * body.velocity + a1 * p.dt * body.k[1]
    
    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
    
    # k3
    for body in p.all_bodies:
        body.temp_position = body.position + delta3 * p.dt * body.velocity + p.dt * (b1*body.k[1] + b2 * body.k[2])

    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
    
    # k4
    for body in p.all_bodies:
        body.temp_position = body.position + delta4 * p.dt * body.velocity + p.dt * (c1*body.k[1] + c2 * body.k[2] + c3 * body.k[3])

    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))

    for body in p.all_bodies:
        body.position = body.position + p.dt * body.velocity + p.dt * (alpha1 * body.k[1] + alpha2 * body.k[2] + alpha3 * body.k[3] + alpha4 * body.k[4])
        body.velocity = body.velocity + beta1 * body.k[1] + beta2 * body.k[2] + beta3 * body.k[3] + beta4 * body.k[4]

def Nystrom5(p):
    delta2 = 1/2
    delta3 = 1/3
    delta4 = 2/3
    delta5 = 1
    a1 = 1/8
    b1 = 1/18
    b2 = 0
    c1 = 1/9
    c2 = 0
    c3 = 1/9
    d1 = 0
    d2 = -8/11
    d3 = 9/11
    d4 = 9/22
    alpha1 = 11/120
    alpha2 = -4/15
    alpha3 = 9/20
    alpha4 = 9/40
    alpha5 = 0
    beta1 = 11/120
    beta2 = -8/15
    beta3 = 27/40
    beta4 = 27/40
    beta5 = 11/120

    for body in p.all_bodies:
        body.k = [0]
        body.k.append(p.dt * gravitational_acc(body.position, p))
    
    for body in p.all_bodies:
        body.temp_position = body.position + delta2 * p.dt * body.velocity + a1 * p.dt * body.k[1]
    
    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
    
    # k3
    for body in p.all_bodies:
        body.temp_position = body.position + delta3 * p.dt * body.velocity + p.dt * (b1*body.k[1] + b2 * body.k[2])

    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
    
    # k4
    for body in p.all_bodies:
        body.temp_position = body.position + delta4 * p.dt * body.velocity + p.dt * (c1*body.k[1] + c2 * body.k[2] + c3 * body.k[3])

    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))

    # k5
    for body in p.all_bodies:
        body.temp_position = body.position + delta5 * p.dt * body.velocity + p.dt * (d1*body.k[1] + d2 * body.k[2] + d3 * body.k[3] + d4 * body.k[4])

    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))

    for body in p.all_bodies:
        body.position = body.position + p.dt * body.velocity + p.dt * (alpha1 * body.k[1] + alpha2 * body.k[2] + alpha3 * body.k[3] + alpha4 * body.k[4] + alpha5 * body.k[5])
        body.velocity = body.velocity + beta1 * body.k[1] + beta2 * body.k[2] + beta3 * body.k[3] + beta4 * body.k[4] + beta5 * body.k[5]


def NystromSimos(p):
    delta2 = 0.25475295159
    delta3 = 0.50540316962
    delta4 = 1
    a1 = 0.03244953299
    b1 = 0.03292284493
    b2 = 0.09479333659
    c1 = 0.19014504913
    c2 = 0
    c3 = 0.30985495135 
    alpha1 = 0.10022791553
    alpha2 = 0.18458008559 
    alpha3 = 0.19318290696
    alpha4 = 0.02200909191 
    beta1 = 0.16323603847
    beta2 = 0.01892388531
    beta3 = 0.65237178035
    beta4 = 0.16546832871

    for body in p.all_bodies:
        body.k = [0]
        body.k.append(p.dt * gravitational_acc(body.position, p))
    
    for body in p.all_bodies:
        body.temp_position = body.position + delta2 * p.dt * body.velocity + a1 * p.dt * body.k[1]
    
    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
    
    # k3
    for body in p.all_bodies:
        body.temp_position = body.position + delta3 * p.dt * body.velocity + p.dt * (b1*body.k[1] + b2 * body.k[2])

    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
    
    # k4
    for body in p.all_bodies:
        body.temp_position = body.position + delta4 * p.dt * body.velocity + p.dt * (c1*body.k[1] + c2 * body.k[2] + c3 * body.k[3])

    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))

    for body in p.all_bodies:
        body.position = body.position + p.dt * body.velocity + p.dt * (alpha1 * body.k[1] + alpha2 * body.k[2] + alpha3 * body.k[3] + alpha4 * body.k[4])
        body.velocity = body.velocity + beta1 * body.k[1] + beta2 * body.k[2] + beta3 * body.k[3] + beta4 * body.k[4]


def Nystrom6(p):
    delta2 = 0.1065417886
    delta3 = 0.2130835772
    delta4 = 0.5926723008
    delta5 = 0.916
    delta6 = 0.972
    a1 =  0.5675576359e-2
    b1 =  0.0756743515e-1
    b2 =  0.1513487029e-1
    c1 =  0.1400361674
    c2 =  -0.2544780570
    c3 =  0.2900721177
    d1 =  -1.0216436141
    d2 =  2.6539701073
    d3 = -1.4861590950
    d4 =  0.2733606017
    e1 =  -20.4083294915
    e2 =  50.3143181086
    e3 =  -32.3044178724
    e4 =  2.9494960939
    e5 =  -0.0786748385
    alpha1 =  0.0627170177
    alpha2 = 0
    alpha3 = 0.2596874616
    alpha4 = 0.1587555586
    alpha5 = 0.0191237845
    alpha6 = -0.0002838224
    beta1 = 0.0627170177
    beta2 = 0
    beta3 = 0.3300064074
    beta4 = 0.3897489881
    beta5 = 0.2276641014
    beta6 = -0.0101365146
    
    for body in p.all_bodies:
        body.k = [0]
        body.k.append(p.dt * gravitational_acc(body.position, p))
    
    for body in p.all_bodies:
        body.temp_position = body.position + delta2 * p.dt * body.velocity + a1 * p.dt * body.k[1]
    
    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
        
    
    for body in p.all_bodies:
        body.temp_position = body.position + delta3 * p.dt * body.velocity + p.dt * (b1*body.k[1] + b2 * body.k[2])

    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))

    # k4
    for body in p.all_bodies:
        body.temp_position = body.position + delta4 * p.dt * body.velocity + p.dt * (c1*body.k[1] + c2 * body.k[2] + c3 * body.k[3])

    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
    
    # k5
    for body in p.all_bodies:
        body.temp_position = body.position + delta5 * p.dt * body.velocity + p.dt * (d1*body.k[1] + d2 * body.k[2] + d3 * body.k[3] + d4 * body.k[4])

    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
    
    # k6
    for body in p.all_bodies:
        body.temp_position = body.position + delta6 * p.dt * body.velocity + p.dt * (e1*body.k[1] + e2 * body.k[2] + e3 * body.k[3] + e4 * body.k[4] + e5 * body.k[5])

    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))

    for body in p.all_bodies:
        body.position = body.position + p.dt * body.velocity + p.dt * (alpha1 * body.k[1] + alpha2 * body.k[2] + alpha3 * body.k[3] + alpha4 * body.k[4] + alpha5 * body.k[5] + alpha6 * body.k[6])
        body.velocity = body.velocity + beta1 * body.k[1] + beta2 * body.k[2] + beta3 * body.k[3] + beta4 * body.k[4] + beta5 * body.k[5] + beta6 * body.k[6]

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

def gravitational_acc_temp(position, p):
    sum_acc = vector(0, 0, 0)
    # calculate the gravitational acceleration from all other bodies
    for body in p.bodies:
        # distance to the other body
        r_vec = body.temp_position - position
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
