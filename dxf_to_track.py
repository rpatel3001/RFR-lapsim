"""Convert a DXF file into a text file."""

import ezdxf
from math import sin, cos, radians, inf


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


dwg = ezdxf.readfile("endurancemichigan2018.dxf")
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


out = open("endurancemichigan2018.csv", 'w')
for sec in sections:
    if sec['type'] == 'LINE':
        out.write("straight,%f\n" % dist(sec['start'], sec['end']))
    elif sec['type'] == 'ARC':
        out.write("turn,%f,%f\n" % (sec['radius'], sec['angle']))
