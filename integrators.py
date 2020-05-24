from tools import conVec, vector
from types import SimpleNamespace

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

def Yoshida6(p):
    # velocity coeffs
    d = [0.784513610477,0.235573213359,-1.17767998417887,1.3151863206857402,-1.17767998417887,0.235573213359,0.784513610477]
    # position coeffs
    c = [0.3922568052385,0.510043411918,-0.471053385409935,0.068753168253435026,0.068753168253435026,-0.471053385409935,0.510043411918,0.3922568052385]
    for i in range(len(d)):
        for body in p.all_bodies:
            body.position += c[i] * p.dt * body.velocity
        for body in p.all_bodies:
            body.velocity += d[i] * p.dt * gravitational_acc(body.position, p)
    
    for body in p.all_bodies:
        body.position += c[-1] * p.dt * body.velocity

def Yoshida8(p):
    # velocity coeffs
    d = [1.48819229202922,-2.33864815101035,2.89105148970595,-2.89688250328827,0.00378039588360192,2.89195744315849,-0.00169248587770116,-3.075516961201882,-0.00169248587770116,2.89195744315849,0.00378039588360192,-2.89688250328827,2.89105148970595,-2.33864815101035,1.48819229202922]
    # position coeffs
    c = [0.74409614601461,-0.425227929490565,0.2762016693478,-0.0029155067911599275,-1.4465510537023341,1.4478689195210459,1.4451324786403945,-1.5386047235397915,-1.5386047235397915,1.4451324786403945,1.4478689195210459,-1.4465510537023341,-0.0029155067911599275,0.2762016693478,-0.425227929490565,0.74409614601461]
    for i in range(len(d)):
        for body in p.all_bodies:
            body.position += c[i] * p.dt * body.velocity
        for body in p.all_bodies:
            body.velocity += d[i] * p.dt * gravitational_acc(body.position, p)
    
    for body in p.all_bodies:
        body.position += c[-1] * p.dt * body.velocity

def McAte5(p):
    a1= 0.339839625839110000
    a2= -0.088601336903027329
    a3= 0.5858564768259621188
    a4= -0.603039356536491888
    a5= 0.3235807965546976394
    a6= 0.4423637942197494587
    b1= 0.1193900292875672758
    b2= 0.6989273703824752308
    b3= -0.1713123582716007754
    b4= 0.4012695022513534480
    b5= 0.0107050818482359840
    b6= -0.0589796254980311632
    # velocity coeffs
    d = [a1, a2, a3, a4, a5, a6]
    # position coeffs
    c = [b1, b2, b3, b4, b5, b6]
    for i in range(len(c)):
        for body in p.all_bodies:
            body.position += c[i] * p.dt * body.velocity
        for body in p.all_bodies:
            body.velocity += d[i] * p.dt * gravitational_acc(body.position, p)

def KahanLi8(p):
    a1=0.74167036435061295344822780
    a2=-0.40910082580003159399730010
    a3=0.19075471029623837995387626
    a4=-0.57386247111608226665638773
    a5=0.29906418130365592384446354
    a6=0.33462491824529818378495798
    a7=0.31529309239676659663205666
    a8=-0.79688793935291635401978884
    a9=a7
    a10=a6
    a11=a5
    a12=a4
    a13=a3
    a14=a2
    a15=a1

    b1= a1/2
    b2=(a1+a2)/2
    b3=(a2+a3)/2
    b4=(a3+a4)/2
    b5=(a4+a5)/2
    b6=(a5+a6)/2
    b7=(a6+a7)/2
    b8=(a7+a8)/2
    b9=b8
    b10=b7
    b11=b6
    b12=b5
    b13=b4
    b14=b3
    b15=b2
    b16=b1

    d = [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15]
    c = [b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16]
    for i in range(len(d)):
        for body in p.all_bodies:
            body.position += c[i] * p.dt * body.velocity
        for body in p.all_bodies:
            body.velocity += d[i] * p.dt * gravitational_acc(body.position, p)
    
    for body in p.all_bodies:
        body.position += c[-1] * p.dt * body.velocity

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
    """
    for body in p.all_bodies:
        body.temp_position = body.position + delta2 * p.dt * body.velocity + a1 * p.dt * body.k[1]
    
    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
        
    
    for body in p.all_bodies:
        body.temp_position = body.position + delta3 * p.dt * body.velocity + p.dt * (b1*body.k[1] + b2 * body.k[2])

    for body in p.all_bodies:
        body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
    """

    coeffs_list = [
    [a1],
    [b1, b2],
    ]
    delta = [delta2, delta3]
    for i, coeffs in enumerate(coeffs_list):
        for body in p.all_bodies:
            coeff_sum = vector(0,0,0)
            for j in range(len(coeffs)):
                coeff_sum += coeffs[j] * body.k[j+1]
            body.temp_position = body.position + delta[i] * p.dt * body.velocity + p.dt * coeff_sum
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
    """
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
    """

    coeffs_list = [
    [a1],
    [b1, b2],
    [c1, c2, c3]
    ]
    delta = [delta2, delta3, delta4]
    for i, coeffs in enumerate(coeffs_list):
        for body in p.all_bodies:
            coeff_sum = vector(0,0,0)
            for j in range(len(coeffs)):
                coeff_sum += coeffs[j] * body.k[j+1]
            body.temp_position = body.position + delta[i] * p.dt * body.velocity + p.dt * coeff_sum
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
    """
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
    """

    coeffs_list = [
    [a1],
    [b1, b2],
    [c1, c2, c3],
    [d1, d2, d3, d4]
    ]
    delta = [delta2, delta3, delta4, delta5]
    for i, coeffs in enumerate(coeffs_list):
        for body in p.all_bodies:
            coeff_sum = vector(0,0,0)
            for j in range(len(coeffs)):
                coeff_sum += coeffs[j] * body.k[j+1]
            body.temp_position = body.position + delta[i] * p.dt * body.velocity + p.dt * coeff_sum
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

def RKN64(p):
    dt_min = 1e-1
    dt_max = 5
    tol = 1e-3

    delta2 = 1/10
    delta3 = 3/10
    delta4 = 7/10
    delta5 = 17/25
    delta6 = 1
    a1 =  1/200
    b1 =  -1/2200
    b2 =  1/22
    c1 =  637/6600
    c2 =  -7/110
    c3 =  7/33
    d1 =  225437/1968750
    d2 =  -30073/281250
    d3 =  65569/281250
    d4 =  -9367/984375
    e1 =  151/2142
    e2 =  5/116
    e3 =  385/1368
    e4 =  55/168
    e5 =  -6250/28101
    
    # 6th order
    alpha1 =  151/2142
    alpha2 = 5/116
    alpha3 = 385/1368
    alpha4 = 55/168
    alpha5 = -6250/28101
    alpha6 = 0
    beta1 = 151/2142
    beta2 = 25/522
    beta3 = 275/684
    beta4 = 275/252
    beta5 = -78125/112404
    beta6 = 1/12
    
    # fourth order
    alpha1_hat =  1349/157500
    alpha2_hat = 7873/50000
    alpha3_hat = 192199/900000
    alpha4_hat = 521683/2100000
    alpha5_hat = -16/125
    alpha6_hat = 0
    beta1_hat = 1349/157500
    beta2_hat = 7873/45000
    beta3_hat = 27457/90000
    beta4_hat = 521683/630000
    beta5_hat = -2/5
    beta6_hat = 1/12

    coeffs_list = [
        [a1],
        [b1, b2],
        [c1, c2, c3],
        [d1, d2, d3, d4],
        [e1, e2, e3, e4, e5],
        ]
    delta = [delta2, delta3, delta4, delta5, delta6]
    alpha = [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6]
    alpha_hat = [alpha1_hat, alpha2_hat, alpha3_hat, alpha4_hat, alpha5_hat, alpha6_hat]
    beta = [beta1, beta2, beta3, beta4, beta5, beta6]


    while True:
        for body in p.all_bodies:
            body.k = [0]
            body.k.append(p.dt * gravitational_acc(body.position, p))
        """
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
        """
        
        for i, coeffs in enumerate(coeffs_list):
            for body in p.all_bodies:
                coeff_sum = vector(0,0,0)
                for j in range(len(coeffs)):
                    coeff_sum += coeffs[j] * body.k[j+1]
                body.temp_position = body.position + delta[i] * p.dt * body.velocity + p.dt * coeff_sum
            for body in p.all_bodies:
                body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
        total_error = 0
        for body in p.all_bodies:
            alpha_hat_sum = vector(0,0,0)
            alpha_sum = vector(0,0,0)
            beta_sum = vector(0,0,0)
            for i in range(len(alpha)):
                alpha_hat_sum += alpha_hat[i] * body.k[i+1]
                alpha_sum += alpha[i] * body.k[i+1]
                beta_sum += beta[i] * body.k[i+1]
            temp_pos = body.position + p.dt * body.velocity + p.dt * alpha_hat_sum #(alpha1_hat * body.k[1] + alpha2_hat * body.k[2] + alpha3_hat * body.k[3] + alpha4_hat * body.k[4] + alpha5_hat * body.k[5] + alpha6_hat * body.k[6])
            body.temp_position = body.position + p.dt * body.velocity + p.dt * alpha_sum #(alpha1 * body.k[1] + alpha2 * body.k[2] + alpha3 * body.k[3] + alpha4 * body.k[4] + alpha5 * body.k[5] + alpha6 * body.k[6])
            body.velocity = body.velocity + beta_sum #beta1 * body.k[1] + beta2 * body.k[2] + beta3 * body.k[3] + beta4 * body.k[4] + beta5 * body.k[5] + beta6 * body.k[6]

            total_error += (body.temp_position - temp_pos).mag

        p.next_dt = p.dt * 0.84 * (tol * p.dt / total_error) ** 0.25
        if p.next_dt < dt_min:
            p.next_dt = dt_min
        elif p.next_dt > dt_max:
            p.next_dt = dt_max
        if total_error < tol:
            break
        else:
            p.dt = p.next_dt
    for body in p.all_bodies:
        body.position = body.temp_position

def RKN12(p):
    dt_min = 1e-2
    dt_max = 50
    tol = 1e-6

    """
    coeffs_list = [
        [1/5000],
        [1/3750, 1/1875],
        [7/2400, -1/240, 1/160],
        [2/1215, 0, 4/729, 32/18225],
        [152/78125, 0, 1408/196875, 2048/703125, 432/546875],
        [29/51200, 0, 341/387072, -151/345600, 243/716800, -11/110592],
        [37/12000, 0, 0, 2/1125, 27/10000, 5/3168, 224/20625],
        [100467472123373/27511470744477696, 0, 101066550784375/25488568483854336, 49478218404275/15475202293768704, 21990175014231/2674726322379776, -3576386017671875/2723635603703291904, 16163228153/1654104722787, 38747524076705/10316801529179136],
        [62178936641284701329/16772293867250014666848, 0, 46108564356250/9072835168325103, 1522561724950/1296119309760729, -45978886013453735443/2174186242050927827184, 299403512366617849203125/4981371278573254356053856, 15571226634087127616/774466927638876610083, -133736375367792139885/4717207650164066625051, 7461389216/501451974639],
        [501256914705531962342417557181/14270506505142656332600844507392, 0, -1143766215625/132752960853408, -6864570325/1185294293334, 194348369382310456605879163404183/99893545535998594328205911551744, -94634958447010580589908066176109375/27549212808177898050085930321520256, -17006472665356285286219618514/155584463413110817059022733377, 33530528814694461893884349656345/14270506505142656332600844507392, -13439782155791134368/17777268379678341919, 1441341768767571/13159456712985856],
        [105854110734231079069010159870911189747853/5156624149476760916008179453333467046288864, 0, -144579793509250000/19842290513127000261, -101935644099967250/48188419817594143491, 1585474394319811696785932424388196965/1709257457318830856936350991091849456, -843499776333774172853009613469456309715703125/510505790798199330684809765880013237582597536, -15057703799298260121553794369056896088480/714327132646734138085088291809720015274157, 1749840442221344572962864758990584360232600/1450300542040339007627300471250037606768743, -11255775246405733991656178432768/27206626483067760480757659602193, 669010348769579696/7368057640845834597, 4598083098752/858563707934367],
        [-1639758773684715326849438048667467886824967397/11447568726280607813664651120965112496134881280, 0, 3942453384375/314673684985856, 11737114158175/1719466921529856, -23710715033675876683332701739887457/4940189888325748664958546898558976, 498150575499633273684774666731162498301909124515625/87415924307623977386706008889913792042985180430336, 64881557768202140428371179540010005713998551/85896810580242200654071863296887242202224768, -2336309182318568698279006266321563486172654055/18316109962048972501863441793544179993815810048, -493399374030747471036018890494175/251658285736841065236836942273664, 418285003077108927126515545155/455369916679568501838710898688, -15171723902781457/63532954684873728, 1501203688494867/9434957026426880],
        [34188549803371802849576690267872548602326398788953/42496542183406636759747616530102745233754251202880, 0, -18971246281693750/1138830954584356089, -59230464334542700/2765732318276293359, 5147939981309774383134903239728881770043/305929030949718561059100251282184099064, -362572021355026772337065830211467821556305840522907812/324512095420929759624784749347170583153994213035432256, -60305503318319653518547439098565661266182518307816/17856872599361492097414471889911176856851308259643, -1036461878759982363277481306266144563833492657780645/67994467493450618815596186448164392374006801924608, 128398681100219349205889126776607047000/7473801441221286756994805323613917077, -49156374556350058671822606102117/9039888303968618912866414995904, 12253036339964386945/8828680926314891943, -647188390508758231059/1092148506009694282240, 10915833599872/368729913707897],
        [-4939337286263213195547765488387521892799075623007291241961609516532/5408250052307451520718178852915698257207815452080611897685945761264, 0, 7588799849596321243074032368290625/3147217749590114939838670370597819616, 16870665568420512953501332587233725/955405388268427749593882076788623812, -80864251591837801485030858227147601466956843757908779606/54447992506702009927986632715967769032585338753056786562, 4610328329649866588704236006423149172472141907645890762410296050212/2135428689710103309390449198881479603148467934048051598947383737508, 4159963831215576225909381034291748993887819834160487158570788681/1040533184037697645660563795162185415624171583014576682740416336, 738139214212435127943380193414870655354213707189052136566460666444958/259596002510757672994472584939953516345975141699869371088925396540699, -333683433458405281346882867597135977469443722954786270692/132102862435303266640535426836147775872819092781208127980, 42661937996741208687503901295747546613008142604821349179/55162410119399855550108207148248549410926885937244965785, -630755628691078947314733435975762542732598947/333503232300511886435069380727586592765317456, 1522350657470125698997653827133798314909646891/1520094067152619944607524353149267399623188480, 305575414262755427083262606101825880/65839748482572312891297405431209259829, 256624643108055110568255672032710477795/22874609758516552135947898572671559986304],
        [-571597862947184314270186718640978947715678864684269066846/207705506488030390761613596901272001190776700439774478634, 0, 66981514290625/1829501741761029,  43495576635800/4443075658562499,-127865248353371207265315478623656127/10401415428935853634424440540325344, 131656514265807573955723157408023481433806699348396032656/92668695535091962564795912774190176478892159517481612467, 3881494143728609118531066904799685950051960514138645179820/2446349095978358868919950548516272963929118212742344026549, 16292266704968075585259245375842819400619822954470178684291/66288722243155885736983218667976563740242178853010092663614, -43986024977384568043684084266385512680544563954/4922783599524658241955780540171948284522386185, 285912200202585226675651763671663063668290787/65371192072964016939690070594254881767827200, -6776815256667778089672518929/3693654613173093729492918708, 398946554885847045598775476868169/344154261237450078839899047372800, -76630698033396272/4432017119727044925, 28401702316003037/1469612686944417840, 66049942462586341419969330578128801/12691068622536592094919763114637498325],
        [83940754497395557520874219603241359529066454343054832302344735/64192596456995578553872477759926464976144474354415663868673233, 0, 892543892035485503125/51401651664490002607536, -12732238157949399705325/686579204375687891972088, 5290376174838819557032232941734928484252549/357179779572898187570048915214361602000384, 2687322933801750693719999180471745666665021538793817303193221/2863980005760296740624015421425947092438943496681472214589916, -197649786681880330585741729796159873563741413724149351549277865/378029217824623393200881653405474359138017953416246216408422692, -100286075630483975704018828319990067604207336241794360144098685695/20486915674765670626893195919603679319429068544972409068469849579, 8739866119696575810411768434844068608106287881671139259/2282122412587168891929052689609009868137678763277087160, -7922242431969626895355493632206885458496418610471389/748272134517487495468365669337985635214015258726400, 2777643183645212014464950387658055285/1141545470045611737197667093465955392, -1372659703515496442825084239977218110461/1313121960368535725613950174847107891200, 6144417902699179309851023/85608793932459282773805825, 140294243355138853053241/64884622846351585391642880, 168671028523891369934964082754523881107337/24062875279623260368388427013982199424119600, 0]
        ]
    """
    """
    coeffs_list = [
        [0.0002,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0.00026666666666666668,0.00053333333333333336,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0.0029166666666666668,-0.0041666666666666666,0.00625,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0.0016460905349794238,0,0.0054869684499314125,0.0017558299039780521,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0.0019456,0,0.0071517460317460316,0.0029127111111111113,0.00078994285714285709,0,0,0,0,0,0,0,0,0,0,0,0],
        [0.00056640625,0,0.00088097304894179892,-0.00043692129629629631,0.00033900669642857143,-9.9464699074074069E-5,0,0,0,0,0,0,0,0,0,0,0],
        [0.0030833333333333333,0,0,0.0017777777777777779,0.0027,0.0015782828282828283,0.01086060606060606,0,0,0,0,0,0,0,0,0,0],
        [0.00365183937480113,0,0.0039651717140723429,0.0031972582629306284,0.0082214673068554348,-0.0013130926959572379,0.0097715869680648684,0.0037557690692328338,0,0,0,0,0,0,0,0,0],
        [0.003707241068718501,0,0.0050820458545552854,0.001174708002175412,-0.021147629915126993,0.060104636981078811,0.020105734768506186,-0.028350750122933581,0.014879568918581932,0,0,0,0,0,0,0,0],
        [0.035125376560733439,0,-0.008615749195138479,-0.0057914480510079169,1.9455548237826159,-3.4351238674565141,-0.10930701107475221,2.3496383118995166,-0.756009408687023,0.10952897222156927,0,0,0,0,0,0,0],
        [0.020527792537482496,0,-0.0072864467644801782,-0.0021153556079618403,0.92758079687235218,-1.6522824844257369,-0.021079563005686568,1.2065364326207872,-0.41371447700106612,0.090798739828096539,0.0053555526005339849,0,0,0,0,0,0],
        [-0.14324078875545515,0,0.012528703773091817,0.0068260191639698273,-4.7995553955743873,5.6986250439519415,0.75534303695236449,-0.12755487858281084,-1.9605926051117386,0.91856090566352622,-0.23880085505284429,0.15911081357234216,0,0,0,0,0],
        [0.80450192055204894,0,-0.016658527067011247,-0.021415834042629738,16.827235928962466,-1.1172835357176099,-3.3771592972263238,-15.243326655360846,17.179835738215417,-5.4377192398239949,1.3878671618364655,-0.59258277326528108,0.029603873171297354,0,0,0,0],
        [-0.91329676669735815,0,0.0024112725757805182,0.017658122693861741,-1.4851649779720382,2.1589708670045757,3.99791558311788,2.8434151800232237,-2.52593643549416,0.77338785423622369,-1.8913028948478674,1.0014845070224718,0.0046411995991090527,0.011218755022148955,0,0,0],
        [-2.75196297205594,0,0.036611888779154923,0.009789519688231562,-12.293062345886211,1.4207226453937902,1.5866476906789537,0.24577735327595945,-8.93519369440327,4.3736727316134072,-1.8347181765449492,1.1592085289061491,-0.01729025316538392,0.019325977904460768,0.0052044429375549929,0,0],
        [1.3076391847404059,0,0.017364109189745843,-0.018544456454265796,14.811522032867726,0.93831763084824715,-0.52284261999445425,-4.8951280525847647,3.8297096034337925,-10.58738133697598,2.4332304376226275,-1.0453406042575446,0.0717732095086726,0.0021622109708082783,0.0070095957596025132,0,0]
        ]
    """
    coeffs_list = [[0.0002,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.00026666666666666668,0.00053333333333333336,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.0029166666666666668,-0.0041666666666666666,0.00625,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.0016460905349794238,0,0.0054869684499314125,0.0017558299039780521,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.0019456,0,0.0071517460317460316,0.0029127111111111113,0.00078994285714285709,0,0,0,0,0,0,0,0,0,0,0,0],[0.00056640625,0,0.00088097304894179892,-0.00043692129629629631,0.00033900669642857143,-9.9464699074074069E-5,0,0,0,0,0,0,0,0,0,0,0],[0.0030833333333333333,0,0,0.0017777777777777779,0.0027,0.0015782828282828283,0.01086060606060606,0,0,0,0,0,0,0,0,0,0],[0.00365183937480113,0,0.0039651717140723429,0.0031972582629306284,0.0082214673068554348,-0.0013130926959572379,0.0097715869680648684,0.0037557690692328338,0,0,0,0,0,0,0,0,0],[0.003707241068718501,0,0.0050820458545552862,0.001174708002175412,-0.021147629915126993,0.060104636981078811,0.020105734768506189,-0.028350750122933581,0.014879568918581932,0,0,0,0,0,0,0,0],[0.035125376560733439,0,-0.008615749195138479,-0.0057914480510079169,1.9455548237826159,-3.4351238674565137,-0.10930701107475221,2.3496383118995166,-0.75600940868702293,0.10952897222156927,0,0,0,0,0,0,0],[0.020527792537482496,0,-0.00728644676448018,-0.0021153556079618403,0.92758079687235218,-1.6522824844257367,-0.021079563005686568,1.2065364326207872,-0.41371447700106612,0.090798739828096539,0.0053555526005339849,0,0,0,0,0,0],[-0.14324078875545515,0,0.012528703773091817,0.0068260191639698273,-4.7995553955743873,5.6986250439519415,0.75534303695236449,-0.12755487858281084,-1.9605926051117384,0.91856090566352622,-0.23880085505284432,0.15911081357234216,0,0,0,0,0],[0.80450192055204894,0,-0.016658527067011247,-0.021415834042629735,16.827235928962466,-11.172835357176098,-3.3771592972263238,-15.243326655360846,17.179835738215417,-5.4377192398239949,1.3878671618364655,-0.59258277326528119,0.029603873171297354,0,0,0,0],[-0.913296766697358,0,0.0024112725757805178,0.017658122693861741,-14.851649779720384,2.1589708670045757,3.99791558311788,28.434151800232232,-25.2593643549416,7.7338785423622376,-1.8913028948478674,1.0014845070224718,0.0046411995991090518,0.011218755022148957,0,0,0],[-0.27519629720559396,0,0.036611888779154923,0.009789519688231562,-12.293062345886211,14.207226453937903,1.5866476906789537,2.4577735327595946,-8.9351936944032726,4.3736727316134072,-1.8347181765449492,1.1592085289061491,-0.017290253165383924,0.019325977904460768,0.0052044429375549929,0,0],[1.3076391847404059,0,0.017364109189745843,-0.018544456454265796,14.811522032867726,9.38317630848247,-5.2284261999445425,-48.95128052584765,38.297096034337926,-10.58738133697598,2.4332304376226275,-1.0453406042575444,0.0717732095086726,0.0021622109708082783,0.0070095957596025141,0,0]]

    delta = [0.02,0.04,0.1,0.13333333333333333,0.16,0.05,0.2,0.25,0.33333333333333331,0.5,0.55555555555555558,0.75,0.8571428571428571,0.94521622227201429,1,1]
    # high order
    #alpha = [63818747/5262156900, 0, 0, 0, 0, 0, 22555300000000/261366897038247, 1696514453125/6717619827072, -45359872/229764843, 19174962087/94371046000, -19310468/929468925, 16089185487681/146694672924800, 1592709632/41841694125, 52675701958271/4527711056573100, 12540904472870916741199505796420811396/2692319557780977037279406889319526430375, 0, 0]
    #beta = [63818747/5262156900, 0, 0, 0, 0, 0, 451106000000000/4965971043726693, 8482572265625/26870479308288, -181439488/689294529, 57524886261/188742092000, -38620936/929468925, 144802669389129/586778691699200, 6370838528/41841694125, 368729913707897/4527711056573100, 111940113324845802831946788738852162520696/1316544263754897771229629968877248424453375, -113178587/12362232960, 1/40]
    alpha = [0.012127868517185414,0,0,0,0,0,0.086297462515688747,0.25254695811871469,-0.19741867993268231,0.2031869190789726,-0.020775808077714918,0.10967804874502014,0.038065132526466504,0.01163406880432423,0.0046580297040248785,0,0]
    beta = [0.012127868517185414,0,0,0,0,0,0.090839434227040786,0.31568369764839338,-0.26322490657690972,0.30478037861845891,-0.041551616155429835,0.2467756096762953,0.15226053010586602,0.081438481630269607,0.085025711938908108,-0.0091551896300779631,0.025]
    # low order
    #alpha_hat = [27121957/1594593000, 0, 0, 0, 0, 0, 4006163300000/55441463008113, 9466403125/25445529648, -163199648/406149975, 23359833/69636250, -18491714/140828625, 11052304606701/58344472186000, 1191129152/44377554375, 2033811086741/124730332137000, 3616943474975740389660406409450169802/951830146690244407118982233597812374375, 0, 0]
    alpha_hat = [0.017008701907006991,0,0,0,0,0,0.0722593359308314,0.372026177326753,-0.40182114500930355,0.33545506830135169,-0.13130650107533182,0.18943190661604864,0.026840802040029046,0.016305665605917924,0.0037999883566965944,0,0]

    # coeffs_list [
    # [a1],
    # [b1, b2],
    # [c1, c2, c3],
    # ...
    # ]
    # for i, coeffs in enumerate(coeffs_list):
    #   for body in p.all_bodies:
    #       coeff_sum = vector(0,0,0)
    #       for j in range(len(coeffs)):
    #           coeff_sum += coeffs[j] * body.k[j]
    #       body.temp_position = body.position + delta[i] * p.dt * body.velocity + p.dt * coeff_sum

    while True:
        for body in p.all_bodies:
            body.k = [0]
            body.k.append(p.dt * gravitational_acc(body.position, p))
        
        for i, coeffs in enumerate(coeffs_list):
            for body in p.all_bodies:
                coeff_sum = vector(0,0,0)
                for j in range(len(coeffs)):
                    if coeffs[j] != 0:
                        coeff_sum += coeffs[j] * body.k[j+1]
                body.temp_position = body.position + delta[i] * p.dt * body.velocity + p.dt * coeff_sum
            for body in p.all_bodies:
                body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
        total_error = 0
        for body in p.all_bodies:
            alpha_hat_sum = vector(0,0,0)
            alpha_sum = vector(0,0,0)
            beta_sum = vector(0,0,0)
            for i in range(len(alpha)):
                alpha_hat_sum += alpha_hat[i] * body.k[i+1]
                alpha_sum += alpha[i] * body.k[i+1]
                beta_sum += beta[i] * body.k[i+1]
            temp_pos = body.position + p.dt * body.velocity + p.dt * alpha_hat_sum #(alpha1_hat * body.k[1] + alpha2_hat * body.k[2] + alpha3_hat * body.k[3] + alpha4_hat * body.k[4] + alpha5_hat * body.k[5] + alpha6_hat * body.k[6])
            body.temp_position = body.position + p.dt * body.velocity + p.dt * alpha_sum #(alpha1 * body.k[1] + alpha2 * body.k[2] + alpha3 * body.k[3] + alpha4 * body.k[4] + alpha5 * body.k[5] + alpha6 * body.k[6])
            body.temp_velocity = body.velocity + beta_sum #beta1 * body.k[1] + beta2 * body.k[2] + beta3 * body.k[3] + beta4 * body.k[4] + beta5 * body.k[5] + beta6 * body.k[6]

            total_error += (body.temp_position - temp_pos).mag
        print("error:", total_error)
        print("dt", p.dt)
        if total_error != 0:
            p.next_dt = p.dt * 0.9 * (tol / total_error) ** (1/12)
        else:
            p.next_dt = p.dt * 5
        if p.next_dt < dt_min:
            p.next_dt = dt_min
        elif p.next_dt > dt_max:
            p.next_dt = dt_max
        if total_error < tol:
            break
        else:
            p.dt = p.next_dt
    for body in p.all_bodies:
        body.position = body.temp_position
        body.velocity = body.temp_velocity

def Nystrom11(p):
    
    for body in p.all_bodies:
        body.k = [0]
        body.k.append(p.dt * gravitational_acc(body.position, p))
    
    coeffs_list = [
        [0.005425347222222222],
        [0.0072337962962962963,0.014467592592592593],
        [0.103334,-0.170368,0.218284],
        [0.032575757575757577,0,0.087804878048780483,0.0046193643754619367,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0.047485073675632081,0,0.1717457745568296,0.059995206214105576,0.006253892679998,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0.0068977424777664732,0,0.006890192907847599,0.011635200867787333,-0.013446471598543607,-0.00072666465485779908,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.0013497377614387065,0,0,0.00060752947621059807,-0.00067937348691585273,-4.0998949896422495E-5,0.00046394425654738483,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.0035231835467672635,0,0,0,-2.5843922040583681E-6,2.75833083288014E-7,0.0040718374796769987,0.012407287532676508,0,0,0,0,0,0,0,0,0,0,0,0,0],[-0.013404889441780981,0,0,0,0.047481672902406062,0.000697807581125852,-0.21042727202509562,0.14798101681889583,0.263834728686989,0,0,0,0,0,0,0,0,0,0,0,0],[0.0043388866791676415,0,0,0,0,-1.2442761509373317E-6,0.0073188272587661951,0.017029036725648074,0.0025620634894002437,2.430123168781617E-6,0,0,0,0,0,0,0,0,0,0,0],[0.0012131814428591819,0,0,0,0,-8.0520837190283189E-7,-0.0012705488195759054,0.0021062886900106688,0.0010628632711979337,1.6087868732388769E-6,-0.00030008816299321539,0,0,0,0,0,0,0,0,0,0],[0.0047914575621211351,0,0,0,0,0,0.019972834957831108,0.032634666336705455,-0.001997857874597378,5.7026782348021117E-8,0.0058098823747878021,-0.01316104038363047,0,0,0,0,0,0,0,0,0],[0.0026807599371736013,0,0,0,0,0,0,0.0099317327196297486,0.0019359960271392125,-1.715138004351496E-7,-0.0011462225062375292,-0.00080652278203019088,0.00020442811812559393,0,0,0,0,0,0,0,0],[0.0074503376752169609,0,0,0,0,0,0,0,0.28260651095411948,2.4415810494897561E-5,-0.16867321434047997,0.0657061274483962,0.066234646021490887,-0.14754882356923849,0,0,0,0,0,0,0],[0.03509222722230964,0,0,0,0,1.8389474302247986E-5,-0.081752988166974072,-0.24216436100893776,-0.56640197575623064,0.00017718746282034686,0.50144230682598512,0.30640121278178167,-0.127682378805156,0.32473564200747751,0.036184737962621921,0,0,0,0,0,0],[0.0014249979875386629,0,0,0,0,-1.7028262150795965E-5,0.073803445201194792,0.2242802589974297,0,0.00688606417577786,0.019499587664316604,-0.18840133996971378,0.0441913708483126,0.038180895504889581,0.055480204019077793,0.013471543833326959,0,0,0,0,0],[0.031094852390984942,0,0,0,0,-4.5931563692237728E-5,0.1216757909683287,0.44153474669073012,-0.81750662543379737,0.11059375804343397,-0.69935179356199706,-0.71094405349471579,0.90460400099121607,1.1670873821626189,-0.22082506170188815,0.21976688615518197,0.034753564440463916,0,0,0,0],[0.01980871288613,0,0,0,0,8.215177736719382E-6,-0.035911384656382085,-0.11021769498626295,0.11740095342626759,-0.0093157112425365732,0.20592141797793864,0.275172998996976,-0.072522343522498617,-0.14348542102371206,0.098174520694169,0.04788142501027539,0.025289611447675366,0.00499469981422358,0,0,0],[0.054336854166050025,0,0,0,0,0,0,0,0,0,-0.70281961772876933,-0.075744714415168585,0.59700915641816688,0.50779989239057977,0.0078035137216304288,0.049109538683188535,0.068051376501884059,-0.022406315947171806,0.016860316209610125,0,0],[0.021624965804621487,0,0,0,0,0,0,0,0,0,0.11464551964856939,0.10794748832410697,0.026896200959149876,0.027886520719467719,0.099634628234718067,0.04982293074498375,0.040090376534865739,0.0017967007447241353,0.0095613413618966039,9.3326922896248891E-5,0]
        ]

    delta = [0.10416666666666667,0.20833333333333334,0.55,0.5,0.755618881615018,0.15,0.058323906888760706,0.2,0.68725987009651612,0.25,0.075,0.31,0.16,0.46,0.61,0.76,0.85,0.92,1,1]
    alpha = [0.021624965804621487,0,0,0,0,0,0,0,0,0,0.11464551964856939,0.10794748832410697,0.026896200959149876,0.027886520719467719,0.099634628234718067,0.04982293074498375,0.040090376534865739,0.0017967007447241353,0.0095613413618966039,9.3326922896248891E-5,0]
    beta =  [0.021478209011278296,0,0,0,0,0,0,0,0,0,0.16391106250765011,0.11779812333994047,0.028048739354934391,0.029061956157160108,0.19028575657103503,0.12211492751295235,0.17554978795153864,0.00032935677706047076,0.12873842235351626,0.022683658462933892,0]

    for i, coeffs in enumerate(coeffs_list):
        for body in p.all_bodies:
            coeff_sum = vector(0,0,0)
            for j in range(len(coeffs)):
                if coeffs[j] != 0:
                    coeff_sum += coeffs[j] * body.k[j+1]
            body.temp_position = body.position + delta[i] * p.dt * body.velocity + p.dt * coeff_sum
        for body in p.all_bodies:
            body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
    

    for body in p.all_bodies:
        alpha_sum = vector(0,0,0)
        beta_sum = vector(0,0,0)
        for i in range(len(alpha)):
            alpha_sum += alpha[i] * body.k[i+1]
            beta_sum += beta[i] * body.k[i+1]
        body.position = body.position + p.dt * body.velocity + p.dt * alpha_sum
        body.velocity = body.velocity + beta_sum

def Nystrom8(p):
    
    for body in p.all_bodies:
        body.k = [0]
        body.k.append(p.dt * gravitational_acc(body.position, p))
    
    coeffs_list = [[0.00125,0,0,0,0,0,0,0,0],[0.0016666666666666668,0.0033333333333333335,0,0,0,0,0,0,0],[0.045,-0.09,0.09,0,0,0,0,0,0],[-0.33797985325712432,0.86115164781709841,-0.48613644655235316,0.0879646519923791,0,0,0,0,0],[0.74614626410430862,-1.804347430246737,1.2761518285888953,-0.031220932716737166,0.05827027027027027,0,0,0,0],[-1.2551623461807584,3.2574631870648232,-2.0688666195750889,0.42064784556032508,-0.018382978723404254,0.06930091185410335,0,0,0],[3.7729018300959409,-9.6696131216771946,6.7991544547851674,-0.84903439078238929,0.46717348654708518,-0.049397421524663677,0.028815162556053812,0,0],[0.028092718568909044,0,0.14570932539682541,0.15294312169312169,0.091517857142857137,0.065547052154195012,0.01618992504409171,0,0]]

    delta = [0.05,0.1,0.3,0.5,0.7,0.9,1,1]
    alpha = [0.028092718568909044,0,0.14570932539682541,0.15294312169312169,0.091517857142857137,0.065547052154195012,0.01618992504409171,0,0]
    beta =  [0.028092718568909044,0,0.1618992504409171,0.2184901738473167,0.18303571428571427,0.2184901738473167,0.1618992504409171,0.028092718568909044,0]

    for i, coeffs in enumerate(coeffs_list):
        for body in p.all_bodies:
            coeff_sum = vector(0,0,0)
            for j in range(len(coeffs)):
                if coeffs[j] != 0:
                    coeff_sum += coeffs[j] * body.k[j+1]
            body.temp_position = body.position + delta[i] * p.dt * body.velocity + p.dt * coeff_sum
        for body in p.all_bodies:
            body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
    

    for body in p.all_bodies:
        alpha_sum = vector(0,0,0)
        beta_sum = vector(0,0,0)
        for i in range(len(alpha)):
            alpha_sum += alpha[i] * body.k[i+1]
            beta_sum += beta[i] * body.k[i+1]
        body.position = body.position + p.dt * body.velocity + p.dt * alpha_sum
        body.velocity = body.velocity + beta_sum

def SympNystrom10(p):
    for body in p.all_bodies:
        body.temp_position = body.position #+ 0.2261935861217336 * p.dt * body.velocity
    for body in p.all_bodies:
        body.k = [0]
        body.k.append(p.dt * gravitational_acc(body.temp_position, p))
    
    coeffs_list = [[0.14796172645774042],[0.054460937638007877,-0.041698403523680996,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-0.047502401173752387,-0.087170837500594109,0.13864085925833261,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.15765763877813588,0.004324071160999178,-0.14031787783162042,0.074527971379058824,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.19390513088072819,0.020489309673899871,-0.18960405764513794,0.087695506609692317,0.059506649942065781,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-0.02718318571992864,-0.078109108118698825,0.11101256145245104,0.0073813102597495781,-0.30344880708650473,0.28463897169388042,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-0.069650057722773731,-0.097047993956092246,0.16875531186745613,-0.0080455234546930911,-0.37316566582132149,0.339312719107585,0.037080356539410643,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.020372228292139614,-0.056900892998352756,0.046350855217126646,0.024656645823350865,-0.22537820936193079,0.22341400806397665,-0.041523465793021845,0.041243575006342438,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.20327349365815531,0.024667304067108091,-0.20234234064412962,0.091098728178124755,0.0748864693996101,-0.012061248592519044,-0.2012254596395866,0.12503951477650702,0.077111520281787763,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.29845586738302127,0.067115641434328882,-0.33176303472149088,0.12567538953196947,0.231145119732426,-0.13460328698238475,-0.28433486204033337,0.16864717620995351,0.11724058952177271,0.13000162478889862,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.27702042318680953,0.057556109230528382,-0.30261698638118217,0.11788858966622896,0.19595505829241167,-0.10700633890469903,-0.26561829837812345,0.15882655967557369,0.10820336560439005,0.10072474857193897,0.0093381021979772917,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.104169186910473,-0.019530091760930006,-0.067588943944694446,0.055097356085563154,-0.0878107255592564,0.11553006150844798,-0.11469158380929274,0.079635025624537809,0.03532895664291219,-0.13535826431496414,0.084638734861349318,-0.03909180949430624,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-0.093641862491913433,-0.10774758236762692,0.20137727645898876,-0.016760965992068916,-0.41255244391440443,0.37020083994349162,0.058029030950938559,-0.0109918092878647,-0.048068556859408557,-0.40553179347952634,0.17081282009607793,-0.083828497675954586,0.37887976649323052,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.1171175689343147,-0.013755521669870178,-0.085195025019625961,0.059801082235666742,-0.066553673230558039,0.098859736269041459,-0.1259975877462057,0.085567307383330529,0.040788019055901535,-0.11767315481819363,0.078997922684550359,-0.036163420353216742,-0.024800838843328545,-0.0037437597826608872,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.33577929388734945,0.083760712763775017,-0.38251217633775636,0.13923377674487553,0.29241810571847265,-0.18265513297600761,-0.3169241698179408,0.1857468479542364,0.13297621690584155,0.180978564093619,-0.016259517082352137,0.013288829336770783,-0.44361719688552642,-0.0076278895050973768,0.45425388780223386,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.12386650931651338,-0.010745707068661226,-0.094371645798969866,0.062252752822986916,-0.055474097527315233,0.090170849040856635,-0.13189048976212031,0.088659324273686224,0.0436333851178642,-0.10845534293817716,0.076057825225400666,-0.034637088824482055,-0.037727502767357625,-0.003863642491621863,0.014020434567676213,-0.012600893847034948,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.24773301863302677,0.0444948502972006,-0.26279452685698973,0.10724942792272808,0.14787461935994897,-0.069300426426162243,-0.24004572023072732,0.14540857930730214,0.095855740311965087,0.060723537960747004,0.022096819690991146,-0.0066236011038348043,-0.27497670733886914,-0.0060639068065313981,0.27134413142227254,-0.0052354640629795526,-0.26388499160342194,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.099771292175286325,-0.021491414137322452,-0.06160907007905312,0.05349974393346487,-0.095030645264269353,0.12119210758256307,-0.11085151984619716,0.07762013676301055,0.03347479588360381,-0.1413649800423683,0.086554626602019541,-0.040086431614850863,0.0084235604399418972,-0.0034356345762589771,-0.036035632931984765,-0.014033659095402336,0.051332407832273143,-0.0087981950324227832,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.19327208099501889,0.020206989386358547,-0.18874329159793973,0.087465540170931985,0.05846738692662324,0.00081501669228728414,-0.19249263026091162,0.12045738260011331,0.072894905739091034,-0.013660091102060113,0.045822094371851883,-0.018940415521191422,-0.17066429709874689,-0.00509650683455651,0.1582054801071707,-0.0084738552500725275,-0.1478615067194069,-0.0032383911870929745,0.19424111303915548,-0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.29523541980677914,0.065679423363271663,-0.32738415085627232,0.12450550564597383,0.22585818716247558,-0.13045713896068103,-0.2815228973555895,0.16717173299422874,0.11588284271726114,0.12560308542625215,0.0014029505670562306,0.0041194798090360245,-0.36596100128479631,-0.0069077009745583014,0.37002692749805155,-0.0024108458320316367,-0.36508402275434571,0.0028246182309479163,0.4060625604300363,-0.0018111941400017922,-0,0,0,0,0,0,0,0,0,0,0,0,0],[0.090075379854890877,-0.02581548529832163,-0.048425413766319307,0.049977534266915008,-0.11094819505463047,0.13367505808372579,-0.10238544102148314,0.073177969109381957,0.029386980864673212,-0.15460781383747477,0.090778538144122373,-0.042279246462275162,0.026994742418848819,-0.0032634042454231854,-0.056178188142665619,-0.014610203649554147,0.071988562230853315,-0.0093747395865745944,-0.020142555210680854,0.001833102589133324,0.39295574370364511,0,0,0,0,0,0,0,0,0,0,0,0],[0.053827887752298555,-0.041980723811222323,0.00086076604719820564,0.036809999036281522,-0.17045484499669625,0.18034169808171327,-0.0707355956847805,0.056571235093791535,0.014104969464481235,-0.20411522454232803,0.10656933690602539,-0.05047693186184956,0.096421811256944323,-0.0026195331634142035,-0.13147972623407092,-0.016765575262130572,0.14921015507235344,-0.011530111199151019,-0.095444093302086142,0.0024769736711423057,0.46238281254174063,-0.008197685399574401,-0,0,0,0,0,0,0,0,0,0,0],[0.27491620435295538,0.056617693981376373,-0.29975585305039082,0.11712419538622426,0.19250061203187424,-0.10429727361216716,-0.26378097952424062,0.15786251498013634,0.10731622295417735,0.097850772471036435,0.010254780683273674,-0.00047588737899364461,-0.32704234820788941,-0.0065467668938556135,0.32781522986796496,-0.0036190800473820477,-0.32179600221507987,0.0016163840155975053,0.36385086279994971,-0.0014502600592991041,0.038918653076906916,0.041803359083281519,-0.096314556222751729,0,0,0,0,0,0,0,0,0,0],[0.31738307635580049,0.075556579818769787,-0.35749860346539586,0.13255102910066693,0.262217470766691,-0.15897102102587174,-0.30086133606365129,0.17731864854737242,0.12522033671303126,0.15585271278858259,-0.0082454186669721867,0.0091283648434220578,-0.40838178094809979,-0.0073011140014739434,0.41603707900153225,-0.0010938876801202394,-0.4122673509576118,0.0041415763828593138,0.452072711933517,-0.0022046071669174342,-0.042420779663303473,0.051407611305697216,-0.11481475557299758,0.058001940317546151,0,0,0,0,0,0,0,0,0],[0.22736079034088713,0.035409478861030304,-0.23509414681506638,0.099848859822622982,0.11443001430730031,-0.0430723099822634,-0.2222575137312188,0.13607507354103002,0.0872667702357284,0.032898819213933145,0.030971744234245757,-0.01123095760603336,-0.23595651548613064,-0.005702031048719843,0.2290223034622729,-0.0064468505617068326,-0.22048403260749774,-0.0012113864987272802,0.26505793639425768,-0.00060552421416333417,0.13000448579866566,0.0310482888562418,-0.075597592671779637,-0.06495195325710329,-0.037953566477302865,0,0,0,0,0,0,0,0],[0.044459524974871448,-0.046158718204430543,0.013599049046189879,0.033406777467849084,-0.18583466445424057,0.1924029466742323,-0.062555519884654048,0.052279133770865439,0.010155249953940635,-0.21691068587635431,0.1106505552956495,-0.052595668533094692,0.11436561728931166,-0.0024531211572792117,-0.15094182100581041,-0.017322642855376283,0.16916849850594518,-0.01208717879239673,-0.11490618807382563,0.0026433856772772974,0.48032661857410797,-0.010316422070819528,0.0040812183896241075,-0.31476145834739078,-0.11506508675909062,-0.083795939770164585,-0,0,0,0,0,0,0],[-0.050722848749994526,-0.088607055571651341,0.14301974312355112,-0.00116988388599563,-0.34209331478705646,0.314944985064098,0.020553882516092747,0.00867147233741894,-0.029973819286044315,-0.34691231066525297,0.15211565045080133,-0.074121989669264129,0.29667431864901206,-0.00076237858714620825,-0.3486763048257483,-0.022982438180348414,0.37194485855496046,-0.01774697411736886,-0.31264067189376354,0.004334128247410301,0.66263531993380842,-0.031842743206988967,0.045546313544775938,-0.44476308313628937,-0.15519415599907557,-0.12740360120361108,0.0831094024007468,-0,0,0,0,0,0],[-0.029287404553782783,-0.079047523367850842,0.11387369478324247,0.0066169159797448817,-0.30690325334704216,0.28734803698641231,0.0018373188538827916,0.018492088871798752,-0.020936595368661645,-0.31763543444829329,0.14277754825282404,-0.069274177648964111,0.25561768345964636,-0.0011431404508606619,-0.30414572287324426,-0.021707830055937603,0.32627882628404981,-0.016472365992958052,-0.26811008994125951,0.0039533663836958475,0.62157868474444267,-0.026994931186688952,0.036208211346798649,-0.41548620691932975,-0.14615693208169292,-0.11758298466923127,0.064392838738536837,-0.027596948077685715,0,0,0,0,0],
        [0.14356383172255377,-0.0019613223763924356,-0.12115434765324531,0.069408149560410692,-0.023137469495374031,0.064811636573265247,-0.14908939571494789,0.097683622922834654,0.051937813592816219,-0.081552421561390129,0.067476915589452013,-0.030182368154657871,-0.0754550055347153,-0.0042135297680811895,0.054940194457797388,-0.011429639034768674,-0.041963140694035124,-0.0061941749717891207,0.090975827389782146,0.00088297706647531924,0.290505995750081,0.012096878307617291,-0.039092421316573388,-0.17940319403242658,-0.073282523120215048,-0.038391450618195364,-0.086533875830293844,-0.25013334849083274,0.28376578385166812,0,0,0,0],[0.34137488112494019,0.086256168230304464,-0.39012056805692852,0.14126647163804276,0.30160424885977394,-0.18985914186177832,-0.32181001047517921,0.18831045783523714,0.13533532709513696,0.188621107603172,-0.018697169645276592,0.014554320026990473,-0.45433477202794564,-0.0077272849832437261,0.46587831479217962,0.00033272839548354437,-0.46337944925593155,0.0055681924584630969,0.50191394772416431,-0.0026307781486872173,-0.0883737707431494,0.056833566489265637,-0.125266506551302,0.090770335132135546,0.010114990382105686,0.052235384294207127,-0.25925449059052513,-0.5048041269258764,0.60850750220681615,0.071858322077632056,-0,0,0],[0.13061544969871208,-0.007735892467452263,-0.1035482665783138,0.0647044234103071,-0.044394521824072386,0.081481961812671769,-0.13778339177803495,0.091751341164041947,0.046478751179826874,-0.099237531058160661,0.073117727766250973,-0.033110757295747369,-0.050654166691386754,-0.00398352520058284,0.028040869135352479,-0.012199584027486868,-0.014377930611879025,-0.0069641199645073155,0.06407650206733724,0.0011129816339736692,0.31530683459340958,0.0091684891665277917,-0.033451609139774428,-0.1970883035291971,-0.078741585533204386,-0.044323732376988077,-0.075227871893380871,-0.23346302325142623,0.26250873152296977,-0.0047037261501035911,0.28657230147861468,0,0],[-0.0880462752543227,-0.10525212690109746,0.19376888473981665,-0.01472827109890168,-0.40336630077310309,0.36299683105772085,0.05314319029370014,-0.0084281994068639669,-0.045709446670113146,-0.39788924996997332,0.16837516753315349,-0.082563006985734891,0.36816219135081124,-9.9395478146349669E-5,-0.42621301866688144,-0.025201787694069895,0.45145862013629645,-0.019966323631090341,-0.39017738573489663,0.0049971113564101595,0.73412319263560755,-0.040283760523459736,0.061805830627128078,-0.49574002244100979,-0.17092978338314441,-0.14450327294789397,0.11569871017835419,0.048051845993622853,-0.096463047426060974,-0.084136420659312372,0.58388945279674509,-0.0975162344336452,0]]

    delta = [0.55326237926184918,0.34657926674995865,0.12118972197991972,0.57469515391881332,0.654820087441651,0.16610526490536828,0.072232417467532131,0.27122631371016226,0.67552881516920293,0.88592906425858342,0.83854610224860948,0.45645914912140878,0.0191986308326014,0.48508149479851631,0.96843234639020259,0.5,0.77380641387826643,0.44673762073815082,0.65342073325004135,0.87881027802008027,0.42530484608118668,0.345179912558349,0.83389473509463174,0.92776758253246783,0.72877368628983774,0.32447118483079707,0.11407093574141658,0.16145389775139052,0.54354085087859128,0.98080136916739857,0.51491850520148374,0.031567653609797408]
    alpha = [0.35006009543824695,0.09012949994972079,-0.40192996532924241,0.14442152758322055,0.31586255965134469,-0.2010408746149657,-0.32939358864258067,0.19228957571631033,0.13899702991787025,0.20048351372085754,-0.022480782344253335,0.016518556735197878,-0.47097010145423834,-0.0078815620935486566,0.48392121813551792,0.00084917417680989253,-0.48188241441404744,0.0060846382397894451,0.51995685106750267,-0.0027850552589921479,-0.10500910016944201,0.058797803197473041,-0.12905011925027873,0.10263274124982109,0.01377669320483898,0.056214502175280331,-0.26683806875792659,-0.51598585967906374,0.62276581299838685,0.075013378022809871,-0.011809397272313945,0.097865392417173055,0.43810637069256969]
    beta =  [0.45238717224346719,0.20175041403676403,-0.61511663906054515,0.16433754952046728,0.74267331435731987,-0.58242344731164464,-0.39500619776092066,0.20726050288524825,0.19072728960001209,0.6178777133180694,-0.197077215139308,0.10231129111935987,-0.8664851973737614,-0.008035839203853588,0.9398015671356833,0.02690013604768968,-0.96376482882809489,0.02690013604768968,0.9398015671356833,-0.008035839203853588,-0.8664851973737614,0.10231129111935987,-0.197077215139308,0.6178777133180694,0.19072728960001209,0.20726050288524825,-0.39500619776092066,-0.58242344731164464,0.74267331435731987,0.16433754952046728,-0.61511663906054515,0.20175041403676403,0.45238717224346719]
    for i, coeffs in enumerate(coeffs_list):
        for body in p.all_bodies:
            coeff_sum = vector(0,0,0)
            for j in range(len(coeffs)):
                if coeffs[j] != 0:
                    coeff_sum += coeffs[j] * body.k[j+1]
            body.temp_position = body.position + delta[i] * p.dt * body.velocity + p.dt * coeff_sum
        for body in p.all_bodies:
            body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
    

    for body in p.all_bodies:
        alpha_sum = vector(0,0,0)
        beta_sum = vector(0,0,0)
        for i in range(len(alpha)):
            alpha_sum += alpha[i] * body.k[i+1]
            beta_sum += beta[i] * body.k[i+1]
        body.position = body.position + p.dt * body.velocity + p.dt * alpha_sum
        body.velocity = body.velocity + beta_sum
""" 
def RKN86(p):
    dt_min = 0.1
    dt_max = 2
    tol = 1e-3

    coeffs_list = [[0.00125,0,0,0,0,0,0,0,0],[0.0016666666666666668,0.0033333333333333335,0,0,0,0,0,0,0],[0.045,-0.09,0.09,0,0,0,0,0,0],[-0.33797985325712432,0.86115164781709841,-0.48613644655235316,0.0879646519923791,0,0,0,0,0],[0.74614626410430862,-1.804347430246737,1.2761518285888953,-0.031220932716737166,0.05827027027027027,0,0,0,0],[-1.2551623461807584,3.2574631870648232,-2.0688666195750889,0.42064784556032508,-0.018382978723404254,0.06930091185410335,0,0,0],[3.7729018300959409,-9.6696131216771946,6.7991544547851674,-0.84903439078238929,0.46717348654708518,-0.049397421524663677,0.028815162556053812,0,0],[0.028092718568909044,0,0.14570932539682541,0.15294312169312169,0.091517857142857137,0.065547052154195012,0.01618992504409171,0,0]]

    delta = [0.05,0.1,0.3,0.5,0.7,0.9,1,1]
    alpha = [0.028092718568909044,0,0.14570932539682541,0.15294312169312169,0.091517857142857137,0.065547052154195012,0.01618992504409171,0,0]
    beta =  [0.028092718568909044,0,0.1618992504409171,0.2184901738473167,0.18303571428571427,0.2184901738473167,0.1618992504409171,0.028092718568909044,0]
    alpha_hat = [0.021624965804621487,0,0,0,0,0,0,0,0,0,0,0.11464551964856939]


    while True:
        for body in p.all_bodies:
            body.k = [0]
            body.k.append(p.dt * gravitational_acc(body.position, p))
        
        for i, coeffs in enumerate(coeffs_list):
            for body in p.all_bodies:
                coeff_sum = vector(0,0,0)
                for j in range(len(coeffs)):
                    if coeffs[j] != 0:
                        coeff_sum += coeffs[j] * body.k[j+1]
                body.temp_position = body.position + delta[i] * p.dt * body.velocity + p.dt * coeff_sum
            for body in p.all_bodies:
                body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
        total_error = 0
        for body in p.all_bodies:
            alpha_hat_sum = vector(0,0,0)
            alpha_sum = vector(0,0,0)
            beta_sum = vector(0,0,0)
            for i in range(len(alpha)):
                alpha_hat_sum += alpha_hat[i] * body.k[i+1]
                alpha_sum += alpha[i] * body.k[i+1]
                beta_sum += beta[i] * body.k[i+1]
            temp_pos = body.position + p.dt * body.velocity + p.dt * alpha_hat_sum #(alpha1_hat * body.k[1] + alpha2_hat * body.k[2] + alpha3_hat * body.k[3] + alpha4_hat * body.k[4] + alpha5_hat * body.k[5] + alpha6_hat * body.k[6])
            body.temp_position = body.position + p.dt * body.velocity + p.dt * alpha_sum #(alpha1 * body.k[1] + alpha2 * body.k[2] + alpha3 * body.k[3] + alpha4 * body.k[4] + alpha5 * body.k[5] + alpha6 * body.k[6])
            body.temp_velocity = body.velocity + beta_sum #beta1 * body.k[1] + beta2 * body.k[2] + beta3 * body.k[3] + beta4 * body.k[4] + beta5 * body.k[5] + beta6 * body.k[6]

            total_error += (body.temp_position - temp_pos).mag
        print("error:", total_error)
        if p.dt != 0.5:
            print("dt", p.dt)
        p.next_dt = p.dt * 0.84 * (tol * p.dt / total_error) ** (1/8)
        if p.next_dt < dt_min:
            p.next_dt = dt_min
        elif p.next_dt > dt_max:
            p.next_dt = dt_max
        if total_error < tol:
            break
        else:
            p.dt = p.next_dt
    for body in p.all_bodies:
        body.position = body.temp_position
        body.velocity = body.temp_velocity
"""
def Nystrom12(p):
    
    for body in p.all_bodies:
        body.k = [0]
        body.k.append(p.dt * gravitational_acc(body.position, p))
    
    coeffs_list = [[0.0002,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.00026666666666666668,0.00053333333333333336,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.0029166666666666668,-0.0041666666666666666,0.00625,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.0016460905349794238,0,0.0054869684499314125,0.0017558299039780521,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.0019456,0,0.0071517460317460316,0.0029127111111111113,0.00078994285714285709,0,0,0,0,0,0,0,0,0,0,0,0],[0.00056640625,0,0.00088097304894179892,-0.00043692129629629631,0.00033900669642857143,-9.9464699074074069E-5,0,0,0,0,0,0,0,0,0,0,0],[0.0030833333333333333,0,0,0.0017777777777777779,0.0027,0.0015782828282828283,0.01086060606060606,0,0,0,0,0,0,0,0,0,0],[0.00365183937480113,0,0.0039651717140723429,0.0031972582629306284,0.0082214673068554348,-0.0013130926959572379,0.0097715869680648684,0.0037557690692328338,0,0,0,0,0,0,0,0,0],[0.003707241068718501,0,0.0050820458545552862,0.001174708002175412,-0.021147629915126993,0.060104636981078811,0.020105734768506189,-0.028350750122933581,0.014879568918581932,0,0,0,0,0,0,0,0],[0.035125376560733439,0,-0.008615749195138479,-0.0057914480510079169,1.9455548237826159,-3.4351238674565137,-0.10930701107475221,2.3496383118995166,-0.75600940868702293,0.10952897222156927,0,0,0,0,0,0,0],[0.020527792537482496,0,-0.00728644676448018,-0.0021153556079618403,0.92758079687235218,-1.6522824844257367,-0.021079563005686568,1.2065364326207872,-0.41371447700106612,0.090798739828096539,0.0053555526005339849,0,0,0,0,0,0],[-0.14324078875545515,0,0.012528703773091817,0.0068260191639698273,-4.7995553955743873,5.6986250439519415,0.75534303695236449,-0.12755487858281084,-1.9605926051117384,0.91856090566352622,-0.23880085505284432,0.15911081357234216,0,0,0,0,0],[0.80450192055204894,0,-0.016658527067011247,-0.021415834042629735,16.827235928962466,-11.172835357176098,-3.3771592972263238,-15.243326655360846,17.179835738215417,-5.4377192398239949,1.3878671618364655,-0.59258277326528119,0.029603873171297354,0,0,0,0],[-0.913296766697358,0,0.0024112725757805178,0.017658122693861741,-14.851649779720384,2.1589708670045757,3.99791558311788,28.434151800232232,-25.2593643549416,7.7338785423622376,-1.8913028948478674,1.0014845070224718,0.0046411995991090518,0.011218755022148957,0,0,0],[-0.27519629720559396,0,0.036611888779154923,0.009789519688231562,-12.293062345886211,14.207226453937903,1.5866476906789537,2.4577735327595946,-8.9351936944032726,4.3736727316134072,-1.8347181765449492,1.1592085289061491,-0.017290253165383924,0.019325977904460768,0.0052044429375549929,0,0],[1.3076391847404059,0,0.017364109189745843,-0.018544456454265796,14.811522032867726,9.38317630848247,-5.2284261999445425,-48.95128052584765,38.297096034337926,-10.58738133697598,2.4332304376226275,-1.0453406042575444,0.0717732095086726,0.0021622109708082783,0.0070095957596025141,0,0]]

    delta = [0.02,0.04,0.1,0.13333333333333333,0.16,0.05,0.2,0.25,0.33333333333333331,0.5,0.55555555555555558,0.75,0.8571428571428571,0.94521622227201429,1,1]
    alpha = [0.012127868517185414,0,0,0,0,0,0.086297462515688747,0.25254695811871469,-0.19741867993268231,0.2031869190789726,-0.020775808077714918,0.10967804874502014,0.038065132526466504,0.01163406880432423,0.0046580297040248785,0,0]
    beta = [0.012127868517185414,0,0,0,0,0,0.090839434227040786,0.31568369764839338,-0.26322490657690972,0.30478037861845891,-0.041551616155429835,0.2467756096762953,0.15226053010586602,0.081438481630269607,0.085025711938908108,-0.0091551896300779631,0.025]
    
    for i, coeffs in enumerate(coeffs_list):
        for body in p.all_bodies:
            coeff_sum = vector(0,0,0)
            for j in range(len(coeffs)):
                if coeffs[j] != 0:
                    coeff_sum += coeffs[j] * body.k[j+1]
            body.temp_position = body.position + delta[i] * p.dt * body.velocity + p.dt * coeff_sum
        for body in p.all_bodies:
            body.k.append(p.dt * gravitational_acc_temp(body.temp_position, p))
    

    for body in p.all_bodies:
        alpha_sum = vector(0,0,0)
        beta_sum = vector(0,0,0)
        for i in range(len(alpha)):
            alpha_sum += alpha[i] * body.k[i+1]
            beta_sum += beta[i] * body.k[i+1]
        body.position = body.position + p.dt * body.velocity + p.dt * alpha_sum
        body.velocity = body.velocity + beta_sum

def Kuraev(p):
    for body in p.all_bodies:
        body.position += body.velocity * p.dt + 0.125 * (5 * body.acc - body.old_acc) * p.dt ** 2
    for body in p.all_bodies:
        temp_acc = gravitational_acc(body.position, p)
        body.velocity += 0.125 * p.dt * (3 * temp_acc + 6 * body.acc - body.old_acc)
        body.old_acc = body.acc
        body.acc = temp_acc

def Beeman(p):
    for body in p.all_bodies:
        body.position += body.velocity * p.dt + 1/6 * (4 * body.acc - body.old_acc) * p.dt ** 2
    for body in p.all_bodies:
        temp_acc = gravitational_acc(body.position, p)
        body.velocity += 1/6 * p.dt * (2 * temp_acc + 5 * body.acc - body.old_acc)
        body.old_acc = body.acc
        body.acc = temp_acc

def Optimized_Verlet(p):
    xi = 0.1931833275037836
    for body in p.all_bodies:
        body.position += body.velocity * xi * p.dt
    for body in p.all_bodies:
        body.velocity += p.dt/2 * gravitational_acc(body.position, p)
    for body in p.all_bodies:
        body.position += p.dt * (1-2*xi) * body.velocity
    for body in p.all_bodies:
        body.velocity += p.dt/2 * gravitational_acc(body.position, p)
    for body in p.all_bodies:
        body.position += xi * p.dt * body.velocity


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

def generate_children(tree):
    a = tree.length/2
    # negative z
    tree1 = SimpleNamespace(**{
        'length': a,
        'position': tree.position,
        'CoM': vector(0,0,0),
        'mass': 0,
        'children': [],
    })
    tree2 = SimpleNamespace(**{
        'length': a,
        'position': tree.position + vector(0, a, 0),
        'CoM': vector(0,0,0),
        'mass': 0,
        'children': [],
    })
    tree3 = SimpleNamespace(**{
        'length': a,
        'position': tree.position + vector(a, 0, 0),
        'CoM': vector(0,0,0),
        'mass': 0,
        'children': [],
    })
    tree4 = SimpleNamespace(**{
        'length': a,
        'position': tree.position + vector(a, a, 0),
        'CoM': vector(0,0,0),
        'mass': 0,
        'children': [],
    })
    # positive z
    tree5 = SimpleNamespace(**{
        'length': a,
        'position': tree.position + vector(0, 0, a),
        'CoM': vector(0,0,0),
        'mass': 0,
        'children': [],
    })
    tree6 = SimpleNamespace(**{
        'length': a,
        'position': tree.position + vector(0, a, a),
        'CoM': vector(0,0,0),
        'mass': 0,
        'children': [],
    })
    tree7 = SimpleNamespace(**{
        'length': a,
        'position': tree.position + vector(a, 0, a),
        'CoM': vector(0,0,0),
        'mass': 0,
        'children': [],
    })
    tree8 = SimpleNamespace(**{
        'length': a,
        'position': tree.position + vector(a, a, a),
        'CoM': vector(0,0,0),
        'mass': 0,
        'children': [],
    })
    return [tree1, tree2, tree3, tree4, tree5, tree6, tree7, tree8]

def is_in_node(body, node):
    a = node.length
    if node.position.x <= body.position.x < node.position.x + a and node.position.y <= body.position.y < node.position.y + a and node.position.z <= body.position.z < node.position.z + a:
        return True
    else:
        return False

def insert_body_in_tree(body, tree):
    # if empty
    if len(tree.children) == 0:
        tree.children = [body]
        tree.mass = body.mass
        tree.CoM = body.mass * body.position
        return
    # there is one body in this node
    elif len(tree.children) == 1:
        other_body = tree.children[0]
        tree.children = generate_children(tree)
        tree.CoM = vector(0,0,0)
        tree.mass = 0
        insert_body_in_tree(other_body, tree)
        insert_body_in_tree(body, tree)
        return
    # There are multiple bodies in this tree already
    elif len(tree.children) == 8:
        for node in tree.children:
            if is_in_node(body, node):
                tree.CoM += body.mass * body.position
                tree.mass += body.mass
                insert_body_in_tree(body, node)
                return
    else:
        raise ValueError(f"Something went wrong. The tree had {len(tree.children)} children.")

def fix_CoM(tree):
    n_children = len(tree.children)
    if n_children == 0:
        return
    elif n_children == 1:
        tree.CoM /= tree.mass
        return
    elif n_children == 8:
        #print(tree)
        tree.CoM /= tree.mass
        for node in tree.children:
            fix_CoM(node)
        return
    else:
        raise ValueError(f"Something went wrong. The tree had {len(tree.children)} children.")


# add vectors together for CoM first and in the end divide by the total mass

def construct_BH_tree(p):
    a = 500
    tree = SimpleNamespace(**{
        'length': a,
        'position': vector(-a/2, -a/2, -a/2),
        'CoM': vector(0,0,0),
        'mass': 0,
        'children': [],
    })
    for body in p.all_bodies:
        insert_body_in_tree(body, tree)
    fix_CoM(tree)
    return tree

def BH_accelerate(position, tree, p):
    n_children = len(tree.children)
    if n_children == 0:
        return vector(0,0,0)
    elif n_children == 1:
        other = tree.children[0]
        r_vec = other.position - position
        r = r_vec.mag
        # skip if body is itself
        if r == 0:
            return vector(0,0,0)
        acc = r_vec * other.GM / r ** 3
        return acc
    elif n_children == 8:
        theta = 1
        r_vec = tree.CoM - position
        d = r_vec.mag
        if d == 0:
            return vector(0,0,0)
        s = tree.length
        if s/d < theta:
            acc = r_vec * p.G * tree.mass / d **3
            return acc
        else:
            total_acc = vector(0,0,0)
            for node in tree.children:
                total_acc += BH_accelerate(position, node, p)
            return total_acc
    else:
        raise ValueError(f"Something went wrong. The tree had {len(tree.children)} children.")

def PEFRL_BH(p):
    for body in p.all_bodies:
        body.position += p.Epsilon * p.dt * body.velocity
    tree = construct_BH_tree(p)
    for body in p.all_bodies:
        body.velocity += (
            (1 - 2 * p.Lambda) * p.dt / 2 * BH_accelerate(body.position, tree, p)
        )
    for body in p.all_bodies:
        body.position += p.Chi * p.dt * body.velocity
    tree = construct_BH_tree(p)
    for body in p.all_bodies:
        body.velocity += p.Lambda * p.dt * BH_accelerate(body.position, tree, p)
    for body in p.all_bodies:
        body.position += (1 - 2 * (p.Chi + p.Epsilon)) * p.dt * body.velocity
    tree = construct_BH_tree(p)
    for body in p.all_bodies:
        body.velocity += p.Lambda * p.dt * BH_accelerate(body.position, tree, p)
    for body in p.all_bodies:
        body.position += p.Chi * p.dt * body.velocity
    tree = construct_BH_tree(p)
    for body in p.all_bodies:
        body.velocity += (
            (1 - 2 * p.Lambda) * p.dt / 2 * BH_accelerate(body.position, tree, p)
        )
    for body in p.all_bodies:
        body.position += p.Epsilon * p.dt * body.velocity

if __name__ == "__main__":
    p = {
    'all_bodies': [
        SimpleNamespace(**{
            'mass': 1,
            'GM': 1,
            'position':vector(-2.5, -2.5, -2.5)
        }),
        SimpleNamespace(**{
            'mass': 1,
            'GM': 1,
            'position':vector(-2.6, -2.5, -2.6)
        }),
        SimpleNamespace(**{
            'mass': 2,
            'GM': 2,
            'position':vector(2.5, 2.5, 2.5)
        }),
        SimpleNamespace(**{
            'mass': 10,
            'GM': 10,
            'position':vector(0,0,0)
        }),
        SimpleNamespace(**{
            'mass': 10,
            'GM': 10,
            'position':vector(3,3,3)
        }),
        SimpleNamespace(**{
            'mass': 10,
            'GM': 10,
            'position':vector(-1,1,1)
        }),
    ],
    'G': 1,
    }
    p = SimpleNamespace(**p)

    tree = construct_BH_tree(p)
    acc_BH = BH_accelerate(p.all_bodies[0].position, tree, p)
    print(f"BH: {acc_BH}")
    position = p.all_bodies[0].position
    sum_acc = vector(0, 0, 0)
    # calculate the gravitational acceleration from all other bodies
    for body in p.all_bodies:
        # distance to the other body
        r_vec = body.position - position
        r = r_vec.mag
        # skip if body is itself
        if r == 0:
            continue
        acc = r_vec * body.GM / r ** 3
        # add force vector to the sum of forces
        sum_acc += acc
    print(f"Brute: {sum_acc}")
    print(f"Diff: {sum_acc - acc_BH}, {(sum_acc - acc_BH).mag}")
        
