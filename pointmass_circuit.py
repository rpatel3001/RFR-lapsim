"""Point mass simulator for straight line acceleration."""

from math import pi, sqrt
import matplotlib.pyplot as plt
from functools import reduce
import csv

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
track_angle_offset = -pi / 2
dt = .0001  # seconds
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


# read the list of track sections
tracklist = []
# read in the list of track elements
with open(track_file) as csvfile:
    for sec in csv.reader(csvfile):
        if sec[0] == "straight":
            tracklist.append({'type': sec[0],
                              'length': float(sec[1])})
        elif sec[0] == "turn":
            tracklist.append({'type': sec[0],
                              'radius': float(sec[1]),
                              'angle': float(sec[2])})

totaldist = 0
for sec in tracklist:
    if 'radius' in sec.keys():
        totaldist += 2 * pi * sec['radius'] * sec['angle'] / 360
    else:
        totaldist += sec['length']

# initial data required to bootstrap the simulation
d = [{'t': 0,
      'vel': 0,
      'dist': 0,
      'x': 0,
      'y': 0,
      'angle': track_angle_offset}, {}]

print("Simulating")
i = 1
while True:
    s = str(int(round(d[i - 1]['dist'] * 100 / totaldist))) + ('%\t[') + (int(d[i - 1]['dist'] * 50 // totaldist) * '#') + (int(49 - d[i - 1]['dist'] * 50 // totaldist) * ' ' + ']')
    print('\r' + s, end='')

    d[i]['t'] = d[i - 1]['t'] + dt

    d[i]['F_normal'] = VEHICLE_MASS * G

    d[i]['F_long_fric_lim'] = d[i]['F_normal'] * MU_LONG

    d[i]['F_cp'] = d[i]['F_long_fric_lim']

    d[i]['F_net'] = d[i]['F_cp']

    d[i]['A_long'] = d[i]['F_net'] / VEHICLE_MASS

    d[i]['vel'] = d[i - 1]['vel'] + d[i]['A_long'] * dt
    d[i]['len'] = (d[i - 1]['vel'] + d[i]['vel']) * dt / 2
    d[i]['dist'] = d[i - 1]['dist'] + d[i]['len']

    if d[i]['dist'] >= totaldist:
        break
    else:
        i += 1
        d.append({})
print()
d = d[1:]

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
    plt.scatter([x for i, x in enumerate(l['t'])], [x for x in l['vel']])
    plt.show()
