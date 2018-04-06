"""Point mass simulator for straight line acceleration."""

from math import pi, inf, sin, acos
import matplotlib.pyplot as plt
from functools import reduce
import argparse
import numpy as np
from scipy import interpolate


def mu_lat(fn):
    """Allow adding load sensitivity."""
    return 1.8


def mu_long(fn):
    """Allow adding load sensitivity."""
    return 1.3


def dist(x1, y1, x2, y2):
    """Calculate Euclidean distance."""
    return ((x1 - x2)**2 + (y1 - y2)**2)**0.5


def get_point(dist):
    """Return the x, y coords of the track at a given distance."""
    x, y = interpolate.splev(dist, tck)
    return float(x), float(y)


def get_radius(i):
    """Get the instantaneous radius at a point."""
    if i == 0:
        bx, by = get_point(d[i - 1]['dist'] - 2 * d[i - 1]['len'])
    else:
        bx = d[i - 1]['x']
        by = d[i - 1]['y']
    ax, ay = get_point(d[i - 1]['dist'] + d[i - 1]['len'])
    cx, cy = get_point(d[i - 1]['dist'] + 2 * d[i - 1]['len'])

    a = dist(bx, by, cx, cy)
    b = dist(ax, ay, cx, cy)
    c = dist(ax, ay, bx, by)

    pa = acos(max(-1, min(1, (b**2 + c**2 - a**2) / (2 * b * c))))
    try:
        return a / (2 * sin(pi - pa))
    except ZeroDivisionError:
        return inf


def ic_engine(vel):
    """Calculate engine parameters."""
    for g, gr in enumerate(gear_ratios):
        gear = g
        rpm = vel / (2 * pi * tire_radius) * 60 * final_drive * gr
        if rpm < upshift_RPM:
            break
        elif rpm > torque_curve[-1][0]:
            rpm = torque_curve[-1][0]
            continue
        else:
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
            'vel': min(rpm * 2 * pi * tire_radius / (60 * final_drive * gr), vel),
            'gear': gear}


def correct_frame(i):
    """Do reverse calculations to correct a frame for an engine limited velocity."""
    eng = ic_engine(d[i]['vel'])
    if d[i]['vel'] <= eng['vel']:
        return
    d[i]['vel'] = eng['vel']
    d[i]['len'] = (d[i - 1]['vel'] + d[i]['vel']) * dt / 2
    d[i]['dist'] = d[i - 1]['dist'] + d[i]['len']


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filename', help='A numpy file output by dxf_to_tck.py. ')
parser.add_argument('-d', '--delta', type=float, default=0.01, help='The timestep to use for simulation. ')
args = parser.parse_args()

npload = np.load(args.filename)
tck = npload[:3]
totaldist = npload[3]
radii = npload[4]

# simulation parameters
dt = args.delta  # seconds
plot_mode = "track"  # track or time

# constants
G = 9.8  # meters per second
rho = 1.2041  # air density in kg/m^3

# vehicle parameters
VEHICLE_MASS = 190 + 69  # mass with driver in kg
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
torque_curve = [(7000, 36.18),  # (rpm, Nm)
                (8000, 33.59),
                (9000, 31.29),
                (10000, 29.83)]


# initial data required to bootstrap the simulation
d = [{'t': 0,
      'vel': 0,
      'len': 1,
      'dist': 0,
      'x': get_point(0)[0],
      'y': get_point(0)[1]}]
d[0]['r'] = get_radius(0)
d.append({})

print("Simulating")
i = 1
while True:
    s = str(int(round(d[i - 1]['dist'] * 100 / totaldist))) + ('%\t[') + (int(d[i - 1]['dist'] * 50 // totaldist) * '#') + (int(49 - d[i - 1]['dist'] * 50 // totaldist) * ' ' + ']')
    print('\r' + s, end='')

    d[i]['t'] = d[i - 1]['t'] + dt

    eng = ic_engine(d[i - 1]['vel'])
    correct_frame(i - 1)
    if d[i - 1]['dist'] >= totaldist:
        break
    d[i]['gear'] = eng['gear']
    d[i]['T_eng_max'] = eng['torque'] * final_drive * gear_ratios[d[i]['gear']]
    d[i]['F_eng_max'] = d[i]['T_eng_max'] / tire_radius

    d[i]['F_drag'] = 0.5 * rho * A * Cd * d[i - 1]['vel']
    d[i]['F_df'] = 0.5 * rho * A * Cl * d[i - 1]['vel']

    d[i]['F_normal'] = (VEHICLE_MASS * G + d[i]['F_df'])

    d[i]['r'] = [r[1] for r in radii if d[i - 1]['dist'] < r[0]][0]

    d[i]['F_lat_tire_max'] = d[i]['F_normal'] * mu_lat(d[i]['F_normal'])
    d[i]['V_corner_max'] = np.sqrt(d[i]['F_lat_tire_max'] / VEHICLE_MASS * d[i]['r'])
    d[i]['F_lat_vel_max'] = VEHICLE_MASS * d[i - 1]['vel']**2 / d[i]['r']
    d[i]['F_lat'] = min(d[i]['F_lat_tire_max'], d[i]['F_lat_vel_max'])

    d[i]['A_lat'] = d[i]['F_lat'] / VEHICLE_MASS

    d[i]['F_long_fric_lim'] = ((1 - (d[i]['F_lat'] / d[i]['F_lat_tire_max'])**2) * (d[i]['F_normal'] * mu_long(d[i]['F_normal']))**2)**.5

    d[i]['F_long_cp'] = min(d[i]['F_long_fric_lim'], d[i]['F_eng_max'])

    d[i]['F_long_net'] = d[i]['F_long_cp'] - d[i]['F_drag']

    d[i]['A_long'] = d[i]['F_long_net'] / VEHICLE_MASS

    d[i]['vel'] = min(d[i]['V_corner_max'], d[i - 1]['vel'] + d[i]['A_long'] * dt)
    d[i]['len'] = (d[i - 1]['vel'] + d[i]['vel']) * dt / 2
    d[i]['dist'] = d[i - 1]['dist'] + d[i]['len']
    d[i]['x'], d[i]['y'] = get_point(d[i]['dist'])

    d.append({})
    if d[i]['dist'] >= totaldist:
        break
    else:
        i += 1
print()
print("Lap length = %f" % d[i]['dist'])
d = d[1:-1]

l = {key: [item[key] for item in d]
     for key in list(reduce(
         lambda x, y: x.union(y),
         (set(dicts.keys()) for dicts in d)
     ))
     }

print("Lap time = " + str(round(l['t'][-1], 4)))
print("Max velocity = " + str(round(max(l['vel']), 3)))
print("Max lateral Gs = " + str(round(max(l['A_lat']) / G, 3)))
print("Max longitudinal Gs = " + str(round(max(l['A_long']) / G, 3)))

print("Plotting")
if plot_mode == "track":
    plt.set_cmap('jet')
    plt.scatter(l['x'], l['y'], c=[d for i, d in enumerate(l['A_lat'])], s=1)
    plt.axis('equal')
    plt.colorbar()
    plt.show()
elif plot_mode == "time":
    plt.scatter([x for i, x in enumerate(l['t'])], [d / G for i, d in enumerate(l['vel'])], label="vel m/s")
    plt.legend(loc='best')
    plt.show()
