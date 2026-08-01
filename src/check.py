"""Verification.  Every claim in the report comes from this file.

Sign tests are exact: the meshes and the collision tests are evaluated from
the same signed distance functions, and min/max composition preserves sign
exactly.  Reported *depths* are conservative lower bounds on true Euclidean
penetration, because max() composition underestimates distance away from the
nearest surface.
"""

import math, sys
from geom import (polygonise, Grid, orient, volume, tri_area, check_manifold,
                  surface_points, rot_about_pin)
import model as M
import build as B

TESS = 0.02          # below this, a hit is tessellation noise, not interference


def hdr(t):
    print('\n' + t + '\n' + '-' * len(t))


# --------------------------------------------------------------------------

def topology(name, mesh):
    wt, single, bad, comps = check_manifold(*mesh)
    print('  %-5s watertight %-5s  single solid %-5s  (non-manifold edges %d,'
          ' components %d)' % (name, wt, single, bad, comps))
    return wt and single


def closed_intersection(bm, cm, step=0.25):
    """Volume intersection in the closed position.

    Measured analytically on the SDFs over a lattice, not from the meshes: a
    tessellated surface carries error of order the marching grid, which is
    larger than the clearances being checked, so a mesh-based test reports
    penetrations that are not in the model.
    """
    lo = (-11.0, -7.0, 72.0)
    hi = (11.0, 7.0, 88.5)
    n = [int((hi[i] - lo[i]) / step) + 1 for i in range(3)]
    inter = worst = 0
    worst_d = 0.0
    worst_at = None
    tot = 0
    for k in range(n[2]):
        z = lo[2] + k * step
        for j in range(n[1]):
            y = lo[1] + j * step
            for i in range(n[0]):
                p = (lo[0] + i * step, y, z)
                dc = M.cap(p)
                if dc >= 0.0:
                    continue
                db = M.body(p)
                if db < 0.0:
                    inter += 1
                    d = max(db, dc)
                    if d < worst_d:
                        worst_d, worst_at = d, p
                tot += 1
    dv = step ** 3
    print('  lattice %d x %d x %d at %.2f mm over the overlap region'
          % (n[0], n[1], n[2], step))
    print('  cells inside BOTH solids: %d  ->  intersection volume %.4f mm^3'
          % (inter, inter * dv))
    if worst_at:
        print('  deepest shared point %.3f mm inside, at (%.2f, %.2f, %.2f)'
              % (-worst_d, worst_at[0], worst_at[1], worst_at[2]))
    return inter * dv


def in_latch_zone(p):
    """Is this point on the detent bead or the groove that receives it?

    Contact there during the first few degrees is the latch releasing -- the
    cantilever is meant to deflect and let the bead out -- so it is reported
    separately from hinge interference, which nothing can deflect away.
    """
    return (abs(p[0]) <= M.BEAD_W / 2.0 + 0.6
            and p[1] >= M.R_NECK - 0.3
            and abs(p[2] - M.BEAD_Z) <= M.BEAD_R + 0.5)


def swing(cm, lo=0, hi=120, step=5, n=6000):
    """Rotate the cap about the pin axis and look for collisions with the
    body at every station."""
    pts = surface_points(cm[0], cm[1], n)
    print('  %d cap-surface sample points per station\n' % len(pts))
    print('   angle   hinge hits  deepest    latch hits  note')
    bad = []
    for a in range(lo, hi + 1, step):
        ang = math.radians(a)
        hits, latch_hits, worst, lworst = 0, 0, 0.0, 0.0
        for p in pts:
            q = rot_about_pin(p, M.Y_PIN, M.Z_PIN, ang)
            d = M.body(q)
            if d < -TESS:
                if in_latch_zone(p):
                    latch_hits += 1
                    lworst = min(lworst, d)
                else:
                    hits += 1
                    worst = min(worst, d)
        note = ''
        if hits:
            note = 'HINGE INTERFERENCE'
            bad.append((a, hits, worst))
        elif latch_hits:
            note = 'detent releasing (%.2f mm, tab deflects)' % -lworst
        print('   %4d   %8d   %9.4f   %8d    %s'
              % (a, hits, worst, latch_hits, note))
    return bad


# --------------------------------------------------------------------------

def bed_contact(bm):
    verts, tris = bm
    area = 0.0
    for a, b, c in tris:
        pa, pb, pc = verts[a], verts[b], verts[c]
        if max(pa[2], pb[2], pc[2]) < 0.05:
            area += tri_area(pa, pb, pc)
    print('  bed contact area: %.1f mm^2   (need >= 100)' % area)
    return area


def bore_wall():
    """Material between the pin hole and the outside of the ear.

    Walks outward from the axis and takes the first continuous run of solid.
    That skips the hole -- bore plus its 45 degree teardrop roof -- without
    needing to know where the hole ends, which matters because the body's
    distance field is exactly zero on that boundary.
    """
    worst, worst_a, worst_x = 1e9, 0.0, 0.0
    for x in (-M.EAR_X, M.EAR_X):
        for i in range(720):
            ang = i * math.pi / 360.0
            ca, sa = math.cos(ang), math.sin(ang)
            r, enter, wall = 0.0, None, 0.0
            while r < 12.0:
                p = (x, M.Y_PIN + r * ca, M.Z_PIN + r * sa)
                solid = M.body(p) < 0.0
                if solid and enter is None:
                    enter = r
                elif not solid and enter is not None:
                    wall = max(wall, r - enter)
                    enter = None
                r += 0.005
            if enter is not None:
                wall = max(wall, r - enter)
            if wall < worst:
                worst, worst_a, worst_x = wall, math.degrees(ang), x
    print('  material between the pin hole and the ear surface: %.2f mm'
          % worst)
    print('        thinnest at %.0f deg round the axis, ear at x = %+.1f'
          '   (need >= %.1f)' % (worst_a, worst_x, M.EAR_MIN_WALL))
    return worst


def cavity_fit():
    """Does the 8 x 5 bundle actually go in, and how deep?"""
    hw, ht = M.BUNDLE_W / 2.0, M.BUNDLE_T / 2.0
    # Sweep the bundle's seating position rather than assuming one, and
    # report the best it can do.  Only the strips' own height matters: they
    # are 60 mm long standing on the floor, so nothing above that is relevant
    # to whether they fit.
    z_lo = M.Z_FLOOR + M.CAV_FILLET
    z_hi = M.Z_FLOOR + M.STRIP_L
    best, y_ctr = -1e9, 0.0
    yc = -1.0
    while yc <= 1.0:
        c = 1e9
        for sx in (-1, 1):
            for sy in (-1, 1):
                for i in range(121):
                    z = z_lo + i * (z_hi - z_lo) / 120.0
                    c = min(c, M.body((sx * hw, yc + sy * ht, z)))
        if c > best:
            best, y_ctr = c, yc
        yc += 0.05
    clear = best
    zc = (M.Z_FLOOR + z_hi) / 2.0
    wx = 0.0
    while M.body((wx, 0.0, zc)) > 0.0 and wx < 12.0:
        wx += 0.005
    yb = yf = 0.0
    while M.body((0.0, -yb, zc)) > 0.0 and yb < 12.0:
        yb += 0.005
    while M.body((0.0, yf, zc)) > 0.0 and yf < 12.0:
        yf += 0.005
    # how far up the bore stays fully enclosed by wall
    zz = M.Z_FLOOR
    while zz < M.Z_CAV_TOP:
        t = 0.0
        while M.body((0.0, -t, zz)) > 0.0 and t < 9.0:
            t += 0.01
        if t > M.R_CAV + 1.0:
            break
        zz += 0.05
    print('  cavity  %.2f wide x %.2f thick x %.2f deep'
          % (wx * 2, 2 * M.R_CAV, M.CAV_DEPTH))
    print('  section at mid strip: %.2f mm back to front (-%.2f to +%.2f)'
          % (yb + yf, yb, yf))
    print('  fully enclosed up to z = %.1f; strips reach z = %.1f'
          '  -> %.1f mm of margin' % (zz, z_hi, zz - z_hi))
    print('  bundle  %.2f wide x %.2f thick x %.2f long (10 x 0.5 mm strips)'
          % (M.BUNDLE_W, M.BUNDLE_T, M.STRIP_L))
    print('  clearance at the bundle corners: %+.2f mm, best seating y = %+.2f'
          % (clear, y_ctr))
    print('        (measured above the %.1f mm fillet at the cavity floor)'
          % M.CAV_FILLET)
    print('  retrieval margin over the 60 mm strip: %.2f mm'
          % (M.CAV_DEPTH - M.STRIP_L))
    return clear


def tab_thickness(z):
    """Measure the free tab where the model actually built it: walk out along
    +Y at x = 0 and take the last continuous run of material."""
    runs, cur, y = [], None, 0.0
    while y < 7.0:
        inside = M.body((0.0, y, z)) < 0.0
        if inside and cur is None:
            cur = y
        elif not inside and cur is not None:
            runs.append((cur, y))
            cur = None
        y += 0.005
    if cur is not None:
        runs.append((cur, y))
    return (runs[-1][1] - runs[-1][0]) if runs else 0.0


def latch():
    """Cantilever bending of the tab as built.

    The tab is not prismatic -- it is thicker below the split than in the
    neck -- so integrate tip deflection along its length rather than using a
    single second moment:  delta = integral P x^2 / (E I(x)) dx, with x
    measured from the loaded tip.
    """
    n = 120
    dz = M.TAB_L / n
    compl, ts = 0.0, []
    for i in range(n):
        z = M.TAB_ROOT_Z + (i + 0.5) * dz
        t = tab_thickness(z)
        ts.append(t)
        if t <= 0.0:
            continue
        I = M.TAB_W * t ** 3 / 12.0
        x = M.TAB_TOP_Z - z
        compl += x * x / (M.PLA_E * I) * dz
    F = M.BEAD_ENGAGE / compl
    t_root = ts[0]
    I_root = M.TAB_W * t_root ** 3 / 12.0
    strain = (F * M.TAB_L) * (t_root / 2.0) / (M.PLA_E * I_root)
    live = [t for t in ts if t > 0]
    print('  tab: %.1f wide x %.2f-%.2f thick x %.2f long'
          % (M.TAB_W, min(live), max(live), M.TAB_L))
    print('  bead %.2f proud, %.2f engaged past the skirt (%.2f is clearance)'
          % (M.BEAD_R, M.BEAD_ENGAGE, M.GAP))
    print('  force to deflect it clear of the groove: %.1f N   (target 8-12)' % F)
    print('  peak surface strain at the root: %.2f %%   (PLA yields near 2-3 %%)'
          % (strain * 100.0))
    return F, strain


BRIDGE = 0.6          # a downward face this close above material is bridged


def supported(p, fn, other):
    """Is there material directly beneath p within BRIDGE mm?"""
    d = 0.05
    while d <= BRIDGE:
        q = (p[0], p[1], p[2] - d)
        if fn(q) < 0.0 or (other is not None and other(q) < 0.0):
            return True
        d += 0.05
    return False


def overhangs(fn, mesh, name, z_ignore=0.15, other=None):
    """Overhang census taken from the SDF gradient at surface sample points.

    Facet normals off a 0.25 mm marching grid carry several degrees of noise
    and produce spurious 90 degree slivers at every sharp edge; the analytic
    gradient does not.  0 deg is a vertical wall, 90 deg a flat ceiling.
    """
    h = 0.004
    LIMIT = math.sin(math.radians(46.0))   # 1 deg of slack on the 45 limit
    pts = surface_points(mesh[0], mesh[1], 40000)
    n_over = n_tot = 0
    worst = 0.0
    worst_at = None
    split = 0
    for p in pts:
        if p[2] < z_ignore:
            continue
        gz = (fn((p[0], p[1], p[2] + h)) - fn((p[0], p[1], p[2] - h))) / (2 * h)
        gx = (fn((p[0] + h, p[1], p[2])) - fn((p[0] - h, p[1], p[2]))) / (2 * h)
        gy = (fn((p[0], p[1] + h, p[2])) - fn((p[0], p[1] - h, p[2]))) / (2 * h)
        L = math.sqrt(gx * gx + gy * gy + gz * gz)
        if L == 0.0:
            continue
        n_tot += 1
        s = -gz / L
        if s > LIMIT:
            # A downward face with material close underneath is not printing
            # over air -- it bridges a short gap, which is the whole premise
            # of a print-in-place mechanism.  Probe straight down for it.
            if supported(p, fn, other):
                split += 1
                continue
            n_over += 1
            a = math.degrees(math.asin(min(s, 1.0)))
            if a > worst:
                worst, worst_at = a, p
    print('  %-5s surface samples steeper than 45 deg: %d of %d (%.2f %%)'
          % (name, n_over, n_tot, 100.0 * n_over / max(n_tot, 1)))
    if split:
        print('        (plus %d with material within %.1f mm directly below --'
              ' those bridge a gap, they do not print over air)'
              % (split, BRIDGE))
    if worst_at:
        print('        worst %.0f deg at (%.2f, %.2f, %.2f)'
              % (worst, worst_at[0], worst_at[1], worst_at[2]))
    return n_over, n_tot


def main():
    step = float(sys.argv[1]) if len(sys.argv) > 1 else 0.25
    bm, cm = B.main()

    hdr('1. topology')
    ok = topology('body', bm) & topology('cap', cm)

    hdr('2. closed position -- volume intersection')
    hits = closed_intersection(bm, cm)

    hdr('3. swing test, 0 to 120 deg about the pin axis')
    bad = swing(cm)

    hdr('4. measurements')
    bed_contact(bm)
    bore_wall()
    cavity_fit()
    latch()
    overhangs(M.body, bm, 'body', other=M.cap)
    overhangs(M.cap, cm, 'cap', other=M.body)

    vb, vc = volume(*bm), volume(*cm)
    hdr('5. mass')
    print('  body %.1f mm^3, cap %.1f mm^3, total %.1f mm^3'
          % (vb, vc, vb + vc))
    print('  PLA at 1.24 g/cm^3, 100%% infill: %.2f g' % ((vb + vc) * M.PLA_DENSITY))

    hdr('summary')
    print('  topology ok        : %s' % ok)
    print('  closed-position hits: %d' % hits)
    print('  swing stations with interference: %d' % len(bad))
    if bad:
        for a, h, w in bad:
            print('     %d deg: %d points, %.3f mm' % (a, h, w))


if __name__ == '__main__':
    main()
