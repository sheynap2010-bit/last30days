"""Polygonise the two solids and write the STL / OBJ deliverables."""

import sys, time, os
from geom import Grid, polygonise, orient, write_stl, volume, tri_area
import model as M

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out')

BODY_BOX = ((-11.0, -7.0, -1.0), (11.0, 7.0, 83.0))
CAP_BOX = ((-11.0, -7.0, 74.0), (11.0, 7.0, 96.0))


def build(fn, box, step, name):
    t0 = time.time()
    g = Grid(fn, box[0], box[1], step)
    t1 = time.time()
    verts, tris = polygonise(g)
    t2 = time.time()
    tris = orient(verts, tris, fn, step * 0.35)
    print('  %-5s sample %5.1fs  march %5.1fs  %6d verts %7d tris'
          % (name, t1 - t0, t2 - t1, len(verts), len(tris)))
    return verts, tris


def write_obj(path, meshes, names):
    with open(path, 'w') as f:
        off = 0
        for (verts, tris), nm in zip(meshes, names):
            f.write('o %s\n' % nm)
            for v in verts:
                f.write('v %.4f %.4f %.4f\n' % v)
            for a, b, c in tris:
                f.write('f %d %d %d\n' % (a + 1 + off, b + 1 + off, c + 1 + off))
            off += len(verts)


def main():
    step = float(sys.argv[1]) if len(sys.argv) > 1 else 0.25
    print('grid step %.2f mm' % step)
    os.makedirs(OUT, exist_ok=True)
    bm = build(M.body, BODY_BOX, step, 'body')
    cm = build(M.cap, CAP_BOX, step, 'cap')

    n = write_stl(os.path.join(OUT, 'strip_case.stl'), [bm, cm])
    write_obj(os.path.join(OUT, 'strip_case.obj'), [bm, cm], ['body', 'cap'])
    write_stl(os.path.join(OUT, 'strip_case_body.stl'), [bm])
    write_stl(os.path.join(OUT, 'strip_case_cap.stl'), [cm])

    vb, vc = volume(*bm), volume(*cm)
    print('  volume body %8.1f mm^3  cap %8.1f mm^3  total %8.1f'
          % (vb, vc, vb + vc))
    print('  mass at 100%% infill: %.2f g' % ((vb + vc) * M.PLA_DENSITY))
    print('  wrote %d triangles to out/strip_case.stl' % n)
    return bm, cm


if __name__ == '__main__':
    main()
