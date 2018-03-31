"""Point mass simulator for straight line acceleration."""

from math import pi, sin, acos, sqrt
import matplotlib.pyplot as plt
from track_gen import gen_track
from functools import reduce

inf = 1e300


def dist(p1, p2):
    """Return the euclidean distance beetween points."""
    return ((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)**.5


def ic_engine(vel):
    """Calculate engine parameters."""
    for g, r in enumerate(gear_ratios):
        gear = g
        rpm = vel / (2 * pi * tire_radius) * 60 * final_drive * r
        if rpm < upshift_RPM:
            v = vel
            break
        elif rpm > torque_curve[-1][0]:
            rpm = torque_curve[-1][0]
            v = rpm * 2 * pi * tire_radius / (60 * final_drive * r)
            continue
        else:
            v = rpm * 2 * pi * tire_radius / (60 * final_drive * r)
            continue
    if rpm <= torque_curve[0][0]:
        torque = torque_curve[0][1]
        rpm = torque_curve[0][0]
    elif rpm >= torque_curve[-1][0]:
        torque = torque_curve[-1][1]
        rpm = torque_curve[-1][0]
    else:
        for r in torque_curve:
            if rpm > r[0]:
                e1 = r
                break
        for r in torque_curve:
            if rpm < r[0]:
                e2 = r
        torque = e1[1] + (rpm - e1[0]) * (e2[1] - e1[1]) / (e2[0] - e1[0])
    return {'torque': torque,
            'rpm': rpm,
            'vel': min(v, vel),
            'gear': gear}


def correct_frame(i):
    """Do reverse calculations to correct a frame for an engine limited velocity."""
    eng = ic_engine(d[i]['vel'])
    if d[i]['vel'] <= eng['vel']:
        d[i]['A_long_corr'] = d[i]['A_long']
        return

    d[i]['vel'] = eng['vel']
    d[i]['A_long_corr'] = (d[i]['vel']**2 - d[i - 1]['vel']**2) / (2 * d[i]['dist'])


# simulation parameters
track_file = 'acceltrack.csv'
track_angle_offset = 0
dd = .001  # meters
closed_track = False
num_laps = 1
plot_mode = "time"  # track or time

# constants
G = 9.8  # meters per second
rho = 1.2041  # air density in kg/m^3

# vehicle parameters
VEHICLE_MASS = 190 + 69  # mass with driver in kg
MU_LAT = 1.5  # lateral coefficient of friction
MU_LONG = 1.5  # longitudinal coefficient of friction
Cd = .436  # coefficient of drag
Cl = 1.07  # coefficient of lift
A = 3.84  # frontal area in m^2
tire_radius = 0.22098  # in meters

# transmission parameters
upshift_RPM = 9500
final_drive = 61 / 23 * 37 / 13
gear_ratios = [35 / 14,
               30 / 15,
               31 / 19,
               28 / 21,
               23 / 21]
torque_curve = [(2500, 33.72),  # (rpm, Nm)
                (3000, 30.26),
                (4000, 33.60),
                (5000, 32.59),
                (6000, 30.07),
                (7000, 30.59),
                (8000, 33.93),
                (9000, 32.75),
                (10000, 29.83)]


# import track
x, y = gen_track(track_file, dd=dd, angle=track_angle_offset,
                 closed=closed_track)

if closed_track:
    x *= num_laps
    y *= num_laps

# calculate the radius at each point on the track
radii = []
# constant that determines how gradual each corner is
h = 1
print("Generating track")
if closed_track:
    for i in range(len(x)):
        s = str(int(i * 100 / len(x))) + ('%\t[') + ((i * 50 // len(x)) * '#') + ((49 - i * 50 // len(x)) * ' ' + ']')
        print('\r' + s, end='')

        a = ((x[(i + h) % len(x)] - x[i - h])**2 +
             (y[(i + h) % len(x)] - y[i - h])**2)**.5
        b = ((x[(i + h) % len(x)] - x[i])**2 +
             (y[(i + h) % len(x)] - y[i])**2)**.5
        c = ((x[(i) % len(x)] - x[i - h])**2 +
             (y[(i) % len(x)] - y[i - h])**2)**.5
        A = acos(max(-1, min(1, (b**2 + c**2 - a**2) / (2 * b * c))))
        try:
            radii.append(a / (2 * sin(pi - A)))
        except ZeroDivisionError:
            radii.append(inf)
else:
    for i in range(len(x) - 2 * h):
        s = str(int(i * 100 / len(x))) + ('%\t[') + ((i * 50 // len(x)) * '#') + ((49 - i * 50 // len(x)) * ' ' + ']')
        print('\r' + s, end='')

        a = ((x[(i + 2 * h)] - x[i])**2 +
             (y[(i + 2 * h)] - y[i])**2)**.5
        b = ((x[(i + 2 * h)] - x[i + h])**2 +
             (y[(i + 2 * h)] - y[i + h])**2)**.5
        c = ((x[(i + h)] - x[i])**2 +
             (y[(i + h)] - y[i])**2)**.5
        A = acos(max(-1, min(1, (b**2 + c**2 - a**2) / (2 * b * c))))
        try:
            radii.append(a / (2 * sin(pi - A)))
        except ZeroDivisionError:
            radii.append(inf)
    for i in range(h):
        radii.append(radii[-1])
print()

# initial data required to bootstrap the simulation
d = [{'t': 0,
      'vel': 0,
      'dist': 0}]

print("Simulating")
i = 0
while i < len(x) - 1:
    s = str(int(i * 100 / len(x))) + ('%\t[') + ((i * 50 // len(x)) * '#') + ((49 - i * 50 // len(x)) * ' ' + ']')
    print('\r' + s, end='')

    d[i]['len'] = dist((x[i], y[i]), (x[i + 1], y[i + 1]))
    d[i]['dist'] = d[i - 1]['dist'] + d[i]['len']

    eng = ic_engine(d[i - 1]['vel'])
    d[i]['rpm'] = eng['rpm']
    d[i]['gear'] = eng['gear']

    d[i]['T_eng_max'] = eng['torque']
    d[i]['P_eng_max'] = d[i]['T_eng_max'] * d[i]['rpm'] * 2 * pi / 60
    d[i]['F_eng_max'] = d[i]['T_eng_max'] * final_drive * gear_ratios[eng['gear']] / tire_radius

    d[i]['F_drag'] = 0.5 * rho * A * Cd * d[i - 1]['vel']
    d[i]['F_df'] = 0.5 * rho * A * Cl * d[i - 1]['vel']

    d[i]['F_normal'] = VEHICLE_MASS * G + d[i]['F_df']

    d[i]['F_long_fric_lim'] = d[i]['F_normal'] * MU_LONG

    d[i]['F_cp'] = min(d[i]['F_long_fric_lim'], d[i]['F_eng_max'])

    d[i]['F_net'] = d[i]['F_cp'] - d[i]['F_drag']

    d[i]['A_long'] = d[i]['F_net'] / VEHICLE_MASS

    d[i]['vel'] = sqrt(d[i - 1]['vel']**2 + 2 * d[i]['A_long'] * d[i]['dist'])

    # correct the previous frame for engine limitations
    correct_frame(i)

    d[i]['t'] = d[i - 1]['t'] + d[i]['len'] / ((d[i - 1]['vel'] + d[i]['vel']) / 2)

    i += 1
    d.append({})
print()
d = d[:-1]

l = {key: [item[key] for item in d]
     for key in list(reduce(
         lambda x, y: x.union(y),
         (set(dicts.keys()) for dicts in d)
     ))
     }

print("Lap time = " + str(l['t'][-1]))

print("Plotting")
if plot_mode == "track":
    plt.scatter(x, y, c=[i for i, d in enumerate(x)], s=1)
    plt.axis('equal')
    plt.colorbar()
    plt.show()
elif plot_mode == "time":
    plt.scatter(l['t'], [x for x in l['A_long']])
    plt.scatter(l['t'], [x for x in l['A_long_corr']])
    plt.show()
