"""Convert a DXF to a spline."""

import argparse
import ezdxf
from math import sin, cos, tan, atan, radians, inf, pi, atan2
import numpy as np
from matplotlib import pyplot as plt
import warnings
from scipy import interpolate


def dist(p1, p2):
    """Return the euclidean distance between points."""
    return ((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)**.5


def endpoints(e):
    """Return endpoints of a section."""
    if e.dxftype() == 'LINE':
        #        out.write("straight,%f\n" %
        #                  dist((e.dxf.start[0], e.dxf.start[1]), (e.dxf.end[0], e.dxf.end[1])))
        p1 = (e.dxf.start[0], e.dxf.start[1])
        p2 = (e.dxf.end[0], e.dxf.end[1])
        return ((p1, p2, 0),)
    elif e.dxftype() == 'ARC':
        #        out.write("turn,%f,%f\n" %
        #                  (e.dxf.radius, (e.dxf.end_angle - e.dxf.start_angle)))
        c = (e.dxf.center[0], e.dxf.center[1])
        p1 = (c[0] + e.dxf.radius * cos(radians(e.dxf.start_angle)),
              c[1] + e.dxf.radius * sin(radians(e.dxf.start_angle)))
        p2 = (c[0] + e.dxf.radius * cos(radians(e.dxf.end_angle)),
              c[1] + e.dxf.radius * sin(radians(e.dxf.end_angle)))
        q1 = (c[0] + e.dxf.radius * cos(radians(e.dxf.start_angle)),
              c[1] + e.dxf.radius * sin(radians(e.dxf.start_angle)))
        q2 = (c[0] + e.dxf.radius * cos(radians(e.dxf.end_angle)),
              c[1] + e.dxf.radius * sin(radians(e.dxf.end_angle)))
        return ((p1, p2, -1), (q2, q1, 1))


if __name__ == "__main__":
    """Store the spline representation of a DXF in a numpy file."""
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='A DXF file consisting of only LINEs and ARCs. ')
    parser.add_argument('-d', '--delta', type=float, default=.001, help='The distance between points to be fed into the interpolation. ')
    parser.add_argument('-a', '--angle', type=float, default=-pi / 2, help='The angle offset of the track. ')
    parser.add_argument('-c', '--closed', action='store_true', help='Specify that the track is closed. ')
    args = parser.parse_args()

    if args.filename[-4:] != ".dxf":
        print("File must end in '.dxf'")
        exit()

    dwg = ezdxf.readfile(args.filename)
    modelspace = [x for x in dwg.modelspace()]
    e = [x for x in modelspace if x.dxftype() == 'LINE'][0]
    modelspace.remove(e)
    sections = [{"type": e.dxftype(), "start": e.dxf.start, "end": e.dxf.end}]
    del e
    while modelspace:
        td = inf
        tf = None
        tq = None
        for f in modelspace:
            for q1, q2, q3 in endpoints(f):
                d = dist(q1, sections[-1]['end'])
                if d < td:
                    td = d
                    tf = f
                    tq = (q1, q2, q3)
        if tf.dxftype() == 'LINE':
            sections.append({"type": tf.dxftype(), "start": tq[0], "end": tq[1]})
        elif tf.dxftype() == 'ARC':
            ang = tq[2] * (tf.dxf.start_angle - tf.dxf.end_angle)
            if ang < -180:
                ang += 360
            if ang > 180:
                ang -= 360
            sections.append({"type": tf.dxftype(), "start": tq[0], "end": tq[1], "radius": tf.dxf.radius, "angle": ang})
        modelspace.remove(tf)

    tracklist = []
    for sec in sections:
        if sec['type'] == 'LINE':
            tracklist.append({'type': 'straight',
                              'length': dist(sec['start'], sec['end'])})
        elif sec['type'] == 'ARC':
            tracklist.append({'type': 'turn',
                              'radius': sec['radius'],
                              'angle': sec['angle']})

    points = [(0, 0)]
    lens = [0]
    angle = args.angle
    # generate points for each track element
    for sec in tracklist:
        # for straights, use linear interpolation between endpoints
        if sec['type'] == "straight":
            numdiv = round(sec['length'] / args.delta)
            dlen = sec['length'] / numdiv
            for i in range(numdiv):
                lens.append(dlen)
                p1 = points[-1]
                dx = cos(angle) * dlen
                dy = sin(angle) * dlen
                points.append((p1[0] + dx, p1[1] + dy))
        # for turns, calculate points along an arc starting at the current point
        elif sec['type'] == "turn":
            dangle = args.delta / sec['radius']
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
                lens.append(abs(sec['radius'] * dangle))
                cangle += dangle
                angle = (angle + dangle) % (2 * pi)
                dx = sign_corr * cos(cangle) * sec['radius']
                dy = sign_corr * sin(cangle) * sec['radius']
                points.append((cx - dx, cy - dy))

    if args.closed:
        dist = ((points[-1][0] - points[0][0])**2 +
                (points[-1][1] - points[0][1])**2)**.5
        numdiv = round(dist / args.delta)
        if numdiv != 0:
            dlen = dist / numdiv
            angle = atan2(points[0][1] - points[-1][1],
                          points[0][0] - points[-1][0])
            for i in range(numdiv - 1):
                p1 = points[-1]
                dx = cos(angle) * dlen
                dy = sin(angle) * dlen
                points.append((p1[0] + dx, p1[1] + dy))

    x, y = zip(*points)
    lens = np.cumsum(lens)

    x = np.array(x)
    y = np.array(y)
    plt.plot(x, y, label='poly')

    # fit splines to x=f(u) and y=g(u), treating both as periodic. also note that s=0
    # is needed in order to force the spline fit to pass through all the input points.
    warnings.simplefilter("ignore")
    tck, u = interpolate.splprep([x, y], u=lens, s=0, per=True)
    warnings.simplefilter("default")

    tck.append(lens[-1])
    np.save(args.filename[:-4], tck)

    tck2 = np.load(args.filename[:-3] + 'npy')

    # evaluate the spline fits for 1000 evenly spaced distance values
    for i in np.arange(0, tck2[-1], 100):
        xi, yi = interpolate.splev(i, tck2[:-1])
        # plot the result
        plt.scatter(xi, yi, s=10, color='red', label='interp')
    plt.axis('equal')
    plt.legend(loc='best')
    plt.show()
