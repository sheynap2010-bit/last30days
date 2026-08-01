"""Locate, rather than merely count, the problems the checks report."""

import math
from geom import surface_points, tri_area, rot_about_pin
import model as M
import build as B


def clusters(pts, tol=1.5):
    """Crude spatial grouping so a report reads as regions, not point dumps."""
    out = []
    for p, v in pts:
        for c in out:
            if (abs(p[0] - c['c'][0]) < tol and abs(p[1] - c['c'][1]) < tol
                    and abs(p[2] - c['c'][2]) < tol):
                c['n'] += 1
                c['worst'] = min(c['worst'], v)
                if v <= c['worst']:
                    c['c'] = p
                break
        else:
            out.append({'c': p, 'n': 1, 'worst': v})
    out.sort(key=lambda c: c['worst'])
    return out


def closed(bm, cm):
    print('\nclosed-position collisions')
    bad = []
    for p in surface_points(cm[0], cm[1], 20000):
        d = M.body(p)
        if d < -0.02:
            bad.append((p, d))
    for p in surface_points(bm[0], bm[1], 20000):
        d = M.cap(p)
        if d < -0.02:
            bad.append((p, d))
    for c in clusters(bad)[:12]:
        print('   n=%4d worst %7.3f at (%6.2f,%6.2f,%6.2f)'
              % (c['n'], c['worst'], *c['c']))
    if not bad:
        print('   none')


def swing_worst(cm):
    print('\nswing collisions')
    pts = surface_points(cm[0], cm[1], 8000)
    for a in range(0, 121, 5):
        ang = math.radians(a)
        bad = []
        for p in pts:
            q = rot_about_pin(p, M.Y_PIN, M.Z_PIN, ang)
            dd = M.body(q)
            if dd < -0.02:
                bad.append((p, dd))
        if not bad:
            continue
        print('  %3d deg:' % a)
        for c in clusters(bad)[:3]:
            p = c['c']
            q = rot_about_pin(p, M.Y_PIN, M.Z_PIN, ang)
            print('     n=%4d worst %7.3f  cap-local (%6.2f,%6.2f,%6.2f)'
                  ' -> world (%6.2f,%6.2f,%6.2f)'
                  % (c['n'], c['worst'], p[0], p[1], p[2], q[0], q[1], q[2]))


def overhang_where(mesh, name, z_ignore=0.06):
    print('\nsteep overhangs on the %s' % name)
    verts, tris = mesh
    bad = []
    for a, b, c in tris:
        pa, pb, pc = verts[a], verts[b], verts[c]
        if max(pa[2], pb[2], pc[2]) < z_ignore:
            continue
        ux, uy, uz = pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]
        vx, vy, vz = pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]
        nx = uy * vz - uz * vy
        ny = uz * vx - ux * vz
        nz = ux * vy - uy * vx
        L = math.sqrt(nx * nx + ny * ny + nz * nz)
        if L == 0.0:
            continue
        s = -nz / L
        if s > 0.7071:
            ctr = ((pa[0] + pb[0] + pc[0]) / 3.0, (pa[1] + pb[1] + pc[1]) / 3.0,
                   (pa[2] + pb[2] + pc[2]) / 3.0)
            bad.append((ctr, -math.degrees(math.asin(min(s, 1.0)))))
    for c in clusters(bad, 2.0)[:12]:
        print('   n=%4d worst %5.0f deg at (%6.2f,%6.2f,%6.2f)'
              % (c['n'], -c['worst'], *c['c']))
    if not bad:
        print('   none')


if __name__ == '__main__':
    bm, cm = B.main()
    closed(bm, cm)
    swing_worst(cm)
    overhang_where(bm, 'body')
    overhang_where(cm, 'cap')
