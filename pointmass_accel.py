"""Point mass simulator for straight line acceleration."""

from math import pi
import csv
import matplotlib.pyplot as plot


def data_frame(index):
    """Return data in a different format."""
    return {'time': data['time'][index],
            'normal_force': data['normal_force'][index],
            'friction_force': data['friction_force'][index],
            'downforce': data['downforce'][index],
            'drag_force': data['drag_force'][index],
            'net_force': data['net_force'][index],
            'acceleration': data['acceleration'][index],
            'velocity': data['velocity'][index],
            'distance': data['distance'][index],
            'gear': data['gear'][index],
            'engine_speed': data['engine_speed'][index],
            'wheel_torque': data['wheel_torque'][index],
            'engine_torque': data['engine_torque'][index]}


def add_plot(series):
    """Add a specified data series to the plot."""
    labels.append(series)
    plot.plot(data['time'], data[series])


def engine_stuff(vel):
    """Calculate engine parameters."""
    for g, r in enumerate(gear_ratios):
        gear = g
        rpm = vel / (2 * pi * tire_radius) * 60 * final_drive * r
        if rpm < upshift_RPM:
            v = vel
            break
        elif rpm > torque_curve[-1][0]:
            rpm = torque_curve[-1][0]
            v = rpm * 2 * pi * tire_radius / 60 / final_drive / gear_ratios[gear]
        elif rpm > upshift_RPM:
            continue
        else:
            v = rpm * 2 * pi * tire_radius / 60 / final_drive / gear_ratios[gear]
            break
    if rpm <= torque_curve[0][0]:
        torque = torque_curve[0][1]
    elif rpm >= torque_curve[-1][0]:
        torque = torque_curve[-1][1]
    else:
        for r in torque_curve:
            if rpm > r[0]:
                e1 = r
                break
        for r in torque_curve:
            if rpm < r[0]:
                e2 = r
        torque = e1[1] + (rpm - e1[0]) * (e2[1] - e1[1]) / (e2[0] - e1[0])
    return {'wheel_torque': torque * final_drive * gear_ratios[gear],
            'engine_torque': torque,
            'rpm': rpm,
            'velocity': v,
            'gear': gear}


# simulation parameters
DT = .00001  # seconds
finish_time = 10  # seconds
finish_distance = 75  # meters
finish_mode = "distance"  # distance or time

# constants
G = 9.8  # meters per second
rho = 1.2041  # air density in kg/m^3

# vehicle parameters
VEHICLE_MASS = 190 + 69  # mass with driver in kg
MU = 1.5  # coefficient of friction
Cd = .436  # coefficient of drag
Cl = 1.07  # coefficient of lift
A = 3.84  # frontal area in m^2
tire_radius = 0.22098  # in meters

# transmission parameters
upshift_RPM = 10000
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

# stored data
data = {'time': [],
        'normal_force': [],
        'friction_force': [],
        'downforce': [],
        'drag_force': [],
        'net_force': [],
        'acceleration': [],
        'velocity': [],
        'distance': [],
        'gear': [],
        'engine_speed': [],
        'wheel_torque': [],
        'engine_torque': []}


i = 0
while True:
    time = i * DT
    i += 1
    try:
        engine = engine_stuff(data['velocity'][-1])
        data['velocity'][-1] = engine['velocity']
    except IndexError:
        engine = engine_stuff(0)
    try:
        downforce = 0.5 * rho * A * Cl * data['velocity'][-1]
        drag_force = 0.5 * rho * A * Cd * data['velocity'][-1]
    except IndexError:
        downforce = 0
        drag_force = 0
    engine_force = engine['wheel_torque'] / tire_radius
    normal_force = VEHICLE_MASS * G + downforce
    friction_force = MU * normal_force
    net_force = min(friction_force, engine_force) - drag_force
    acceleration = net_force / VEHICLE_MASS
    try:
        velocity = data['velocity'][-1] + data['acceleration'][-1] * DT
        distance = data['distance'][-1] + velocity * DT
    except IndexError:
        velocity = acceleration * DT
        distance = velocity * DT

    # store new values
    data['time'].append(time)
    data['normal_force'].append(normal_force)
    data['friction_force'].append(friction_force)
    data['downforce'].append(downforce)
    data['drag_force'].append(drag_force)
    data['net_force'].append(net_force)
    data['acceleration'].append(acceleration)
    data['velocity'].append(velocity)
    data['distance'].append(distance)
    data['gear'].append(engine['gear'])
    data['engine_speed'].append(engine['rpm'])
    data['wheel_torque'].append(engine['wheel_torque'])
    data['engine_torque'].append(engine['engine_torque'])

    if finish_mode == "distance" and distance > finish_distance \
            or finish_mode == "time" and time > finish_time:
        break

print(data['time'][-1])

labels = []
add_plot('distance')
add_plot('velocity')
add_plot('engine_torque')
add_plot('acceleration')
plot.legend(labels)
plot.show()
