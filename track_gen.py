"""Create a list of X,Y points based on straights and constant radius turns."""

import csv
from math import sin, cos, tan, atan, radians, inf, pi


def gen_track(filename):
    """Turn a list of turns and straights into a list of X,Y points."""
    tracklist = []
    with open(filename) as csvfile:
        for sec in csv.reader(csvfile):
            if sec[0] == "straight":
                tracklist.append({'type': sec[0],
                                  'length': float(sec[1])})
            elif sec[0] == "turn":
                tracklist.append({'type': sec[0],
                                  'radius': float(sec[1]),
                                  'angle': float(sec[2])})

    points = [(0, 0)]
    angle = 0 * pi
    dd = 5
    for sec in tracklist:
        if sec['type'] == "straight":
            numdiv = round(sec['length'] / dd)
            dlen = sec['length'] / numdiv
            for i in range(numdiv):
                p1 = points[-1]
                dx = cos(angle) * dlen
                dy = sin(angle) * dlen
                points.append((p1[0] + dx, p1[1] + dy))
        elif sec['type'] == "turn":
            dangle = dd / sec['radius']
            numdiv = abs(round(radians(sec['angle']) / dangle))
            dangle = radians(sec['angle']) / numdiv
            try:
                cangle = atan(-1 / tan(angle))
            except ZeroDivisionError:
                cangle = atan(inf)
            p1 = points[-1]
            sign_corr = (1 if sin(cangle) * cos(angle) >=
                         0 else -1) * (1 if dangle >= 0 else -1)
            cx = p1[0] + sign_corr * cos(cangle) * sec['radius']
            cy = p1[1] + sign_corr * sin(cangle) * sec['radius']
            for i in range(numdiv):
                dx = sign_corr * cos(cangle) * sec['radius']
                dy = sign_corr * sin(cangle) * sec['radius']
                points.append((cx - dx, cy - dy))
                cangle += dangle
                angle = (angle + dangle) % (2 * pi)

    dist = ((points[-1][0] - points[0][0])**2 +
            (points[-1][1] - points[0][1])**2)**.5
    numdiv = round(dist / dd)
    dlen = dist / numdiv
    for i in range(numdiv):
        p1 = points[-1]
        dx = cos(angle) * dlen
        dy = sin(angle) * dlen
        points.append((p1[0] + dx, p1[1] + dy))

    x, y = zip(*points)
    return (x, y)
