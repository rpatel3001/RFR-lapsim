"""Point mass simulator for straight line acceleration."""

from math import pi, inf, sin, acos
import matplotlib.pyplot as plt
from functools import reduce
import argparse
import numpy as np
from scipy import interpolate
import warnings


def mu_lat(fn):
    """Allow adding load sensitivity."""
    return 1.5


def mu_long(fn):
    """Allow adding load sensitivity."""
    return 1.5


def dist(x1, y1, x2, y2):
    """Calculate Euclidean distance."""
    return ((x1 - x2)**2 + (y1 - y2)**2)**0.5


def get_point(dist):
    """Return the x, y coords of the track at a given distance."""
    x, y = interpolate.splev(dist, tck)
    return float(x), float(y)


def calc_radius(i):
    """Get the instantaneous radius at a point."""
    bx, by = track_points[i - 1]
    ax, ay = track_points[i]
    cx, cy = track_points[i + 1]

    a = dist(bx, by, cx, cy)
    b = dist(ax, ay, cx, cy)
    c = dist(ax, ay, bx, by)

    pa = acos(max(-1, min(1, (b**2 + c**2 - a**2) / (2 * b * c))))
    warnings.simplefilter("ignore")
    try:
        r = a / (2 * sin(pi - pa))
    except ZeroDivisionError:
        r = inf
    warnings.simplefilter("default")
    return r


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
    if d[i]['vel'] <= eng['vel'] + .00001:
        return
    d[i]['vel'] = eng['vel']
    d[i]['A_long'] = 0
    d[i]['dt'] = min([x for x in np.roots([0.5 * d[i]['A_long'], d[i - 1]['vel'], -dd]) if x > 0])


parser = argparse.ArgumentParser()
req = parser.add_argument_group('required named arguments')
req.add_argument('-f', '--filename', required=True, help='A numpy file output by dxf_to_tck.py. ')
parser.add_argument('-d', '--delta', type=float, default=0.01, help='The timestep to use for simulation. ')
args = parser.parse_args()

npload = np.load(args.filename)
tck = npload[:3]
totaldist = npload[3]
radii = npload[4]

# simulation parameters
dd = args.delta  # meters
plot_mode = "time"  # track or time

# constants
G = 9.8  # meters per second
rho = 1.2041  # air density in kg/m^3

# vehicle parameters
VEHICLE_MASS = 190 + 69  # mass with driver in kg
Cd = 0.436  # coefficient of drag
Cl = 1.07  # coefficient of lift
A = 3.84  # frontal area in m^2
tire_radius = 0.2286  # in meters
CG_long = 0.49  # percent rear
CP_long = 0.54  # percent rear
CG_vert = 0.3124
wheelbase = 1.575
trackwidth_front = 1.270
trackwidth_rear = 1.219

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

radii = [inf]
numdiv = int(round(totaldist / dd))
dd = totaldist / numdiv
print("Generating track with step size = %f" % dd)

td = np.linspace(0, totaldist, numdiv)
tx, ty = interpolate.splev(td, tck)
track_points = list(zip(tx, ty))

for i in range(len(track_points)):
    radii.append(calc_radius(i - 2))
radii.append(calc_radius(i - 1))

# initial data required to bootstrap the simulation
d = [{'t': 0,
      'vel': 0,
      'dist': 0,
      'A_long': 0}]

print("Simulating")
for i in range(len(track_points)):
    s = str(int(round(d[i - 1]['dist'] * 100 / totaldist))) + ('%\t[') + (int(d[i - 1]['dist'] * 50 // totaldist) * '#') + (int(49 - d[i - 1]['dist'] * 50 // totaldist) * ' ' + ']')
    print('\r' + s, end='')

    d[i]['len'] = dd
    d[i]['dist'] = d[i - 1]['dist'] + d[i]['len']
    d[i]['x'], d[i]['y'] = track_points[i]

    correct_frame(i - 1)

    eng = ic_engine(d[i - 1]['vel'])
    d[i]['gear'] = eng['gear']
    d[i]['T_eng_max'] = eng['torque'] * final_drive * gear_ratios[d[i]['gear']]
    d[i]['F_eng_max'] = d[i]['T_eng_max'] / tire_radius

    d[i]['F_drag'] = 0.5 * rho * A * Cd * d[i - 1]['vel']
    d[i]['F_df'] = 0.5 * rho * A * Cl * d[i - 1]['vel']

    d[i]['F_normal_front'] = VEHICLE_MASS * G * (1 - CG_long) + d[i]['F_df'] * (1 - CP_long)
    d[i]['F_normal_rear'] = VEHICLE_MASS * G * CG_long + d[i]['F_df'] * CP_long

    d[i]['F_normal_total'] = d[i]['F_normal_front'] + d[i]['F_normal_rear']

    d[i]['radius'] = radii[i]

    d[i]['F_lat_tire_max'] = d[i]['F_normal_total'] * mu_lat(d[i]['F_normal_total'])
    d[i]['V_corner_max'] = np.sqrt(d[i]['F_lat_tire_max'] / VEHICLE_MASS * d[i]['radius'])
    d[i]['F_lat_vel_max'] = VEHICLE_MASS * d[i - 1]['vel']**2 / d[i]['radius']
    d[i]['F_lat'] = min(d[i]['F_lat_tire_max'], d[i]['F_lat_vel_max'])

    d[i]['A_lat'] = d[i]['F_lat'] / VEHICLE_MASS

    d[i]['F_long_fric_lim'] = ((1 - (d[i]['F_lat'] / d[i]['F_lat_tire_max'])**2) * (d[i]['F_normal_rear'] * mu_long(d[i]['F_normal_rear']))**2)**.5

    d[i]['F_long_cp'] = min(d[i]['F_long_fric_lim'], d[i]['F_eng_max'])

    d[i]['F_long_net'] = d[i]['F_long_cp'] - d[i]['F_drag']

    d[i]['A_long'] = d[i]['F_long_net'] / VEHICLE_MASS

    if d[i]['V_corner_max'] < ((d[i - 1]['vel'])**2 + 2 * d[i]['A_long'] * dd)**.5:
        d[i]['vel'] = d[i]['V_corner_max']
        d[i]['A_long'] = 0
    else:
        d[i]['vel'] = ((d[i - 1]['vel'])**2 + 2 * d[i]['A_long'] * dd)**.5

    d[i]['dt'] = min([x for x in np.roots([0.5 * d[i]['A_long'], d[i - 1]['vel'], -dd]) if x > 0])

    d.append({})

d = d[1:-1]
correct_frame(-1)

cum = 0
for f in d:
    cum += f['dt']
    f['t'] = cum

l = {key: [item[key] for item in d]
     for key in list(reduce(
         lambda x, y: x.union(y),
         (set(dicts.keys()) for dicts in d)
     ))
     }

print()
print("Lap length = " + str(round(l['dist'][-1], 2)))
print("Lap time = " + str(round(l['t'][-1], 4)))
print("Max velocity = " + str(round(max(l['vel']), 3)))
print("Max lateral Gs = " + str(round(max(l['A_lat']) / G, 3)))
print("Max longitudinal Gs = " + str(round(max(l['A_long']) / G, 3)))

if plot_mode == "track":
    plt.set_cmap('jet')
    plt.scatter(l['x'], l['y'], c=[min(d, 50) for i, d in enumerate(l['vel'])], s=1)
    plt.axis('equal')
    plt.colorbar()
    plt.show()
elif plot_mode == "time":
    plt.scatter([x for i, x in enumerate(l['t'])], [d for i, d in enumerate(l['vel'])], label="dist m")
    plt.legend(loc='best')
    plt.show()
