"""Analytic SDF kernel + marching-tetrahedra polygoniser + STL writer.

Everything the case is made of is defined as a signed distance function.  That
matters for this design: rounding an edge is `abs(d) - r` or a smooth-min, i.e.
a genuine constant-radius blend evaluated on the surface itself, not a
post-hoc mesh operation.  It also makes the collision checks exact -- testing
whether a point is inside the body is one function evaluation, no ray casting.

Pure standard library: no numpy, no OpenCascade.
"""

from array import array
from math import sqrt, hypot, sin, cos, inf


# --------------------------------------------------------------------------
# scalar helpers
# --------------------------------------------------------------------------

def clamp(x, lo, hi):
    return lo if x < lo else (hi if x > hi else x)


def smin(a, b, k):
    """Polynomial smooth minimum -- a fillet of radius ~k on a union."""
    if k <= 0.0:
        return a if a < b else b
    h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0)
    return b * (1.0 - h) + a * h - k * h * (1.0 - h)


def smax(a, b, k):
    """Smooth maximum -- a fillet on an intersection / subtraction."""
    return -smin(-a, -b, k)


# --------------------------------------------------------------------------
# 2D primitives (used on cross sections)
# --------------------------------------------------------------------------

def sd_stadium_2d(x, y, half_flat, r):
    """Stadium (obround): two semicircles of radius r joined by flats.

    Total width  = 2*half_flat + 2*r,  total thickness = 2*r.
    """
    dx = abs(x) - half_flat
    if dx < 0.0:
        dx = 0.0
    return hypot(dx, y) - r


def sd_segment_2d(px, py, ax, ay, bx, by):
    vx, vy = bx - ax, by - ay
    wx, wy = px - ax, py - ay
    den = vx * vx + vy * vy
    t = 0.0 if den == 0.0 else clamp((wx * vx + wy * vy) / den, 0.0, 1.0)
    return hypot(wx - vx * t, wy - vy * t)


# --------------------------------------------------------------------------
# 3D primitives
# --------------------------------------------------------------------------

def sd_capsule_shell(p, half_flat, r, z0, z1, cap_r):
    """A stadium prism between z0 and z1 whose ends are rounded by cap_r.

    The result is the 'flat oval capsule' silhouette: every cross section is a
    stadium, and the top and bottom are domed with the same profile so the
    outline is symmetric end to end.
    """
    x, y, z = p
    # radial distance from the stadium spine (a line segment along x at y=0)
    dx = abs(x) - half_flat
    if dx < 0.0:
        dx = 0.0
    rad = hypot(dx, y)

    # distance in the (rad, z) half-plane to a rounded-end rod
    if z < z0:
        dz = z0 - z
    elif z > z1:
        dz = z - z1
    else:
        dz = 0.0

    if dz == 0.0:
        return rad - r
    # inside the end-dome region: blend radius r -> cap_r
    # treat the end as a torus-like round of radius cap_r on the rim
    rr = rad - (r - cap_r)
    if rr < 0.0:
        rr = 0.0
    return hypot(rr, dz) - cap_r


def sd_stadium_prism(p, half_flat, r, z0, z1):
    """Stadium prism with square ends (used for cavities and cutters)."""
    x, y, z = p
    d2 = sd_stadium_2d(x, y, half_flat, r)
    dz = max(z0 - z, z - z1)
    ox, oy = max(d2, 0.0), max(dz, 0.0)
    return min(max(d2, dz), 0.0) + hypot(ox, oy)


def sd_box(p, cx, cy, cz, hx, hy, hz):
    qx = abs(p[0] - cx) - hx
    qy = abs(p[1] - cy) - hy
    qz = abs(p[2] - cz) - hz
    ox, oy, oz = max(qx, 0.0), max(qy, 0.0), max(qz, 0.0)
    return min(max(qx, max(qy, qz)), 0.0) + sqrt(ox * ox + oy * oy + oz * oz)


def sd_cyl_x(p, cy, cz, r, x0, x1):
    """Cylinder with its axis parallel to X -- the hinge-family primitive."""
    rad = hypot(p[1] - cy, p[2] - cz) - r
    dx = max(x0 - p[0], p[0] - x1)
    ox, oz = max(rad, 0.0), max(dx, 0.0)
    return min(max(rad, dx), 0.0) + hypot(ox, oz)


def sd_cyl_y(p, cx, cz, r, y0, y1):
    rad = hypot(p[0] - cx, p[2] - cz) - r
    dy = max(y0 - p[1], p[1] - y1)
    orad, ody = max(rad, 0.0), max(dy, 0.0)
    return min(max(rad, dy), 0.0) + hypot(orad, ody)


def sd_teardrop(a, b, r):
    """Cross section of a self-supporting horizontal hole.

    A circle of radius r with a 45 degree gable roof tangent to it, so the
    hole's roof never exceeds 45 degrees from vertical and needs no support.
    The apex sits r*sqrt(2) above the centre -- not 2r.  `b` is the vertical
    coordinate relative to the axis, `a` the horizontal one.
    """
    circ = hypot(a, b) - r
    # The 45 degree faces are tangent to the circle at b = r/sqrt(2), so the
    # roof region starts there and the teardrop is never wider than the
    # circle.  |a| is softened to sqrt(a^2 + k^2) - k so the ridge where the
    # two faces meet carries a small radius instead of a knife edge, which
    # otherwise shows up as a line of spurious horizontal facets.
    k = 0.3
    sa = sqrt(a * a + k * k) - k
    gable = max(r * 0.7071067811865476 - b,
                (sa + b - r * 1.4142135623730951) * 0.7071067811865476)
    return min(circ, gable)


def sd_plane_z(p, z, up=True):
    return (z - p[2]) if up else (p[2] - z)


# --------------------------------------------------------------------------
# rotation about the pin axis (a line parallel to X at (y0, z0))
# --------------------------------------------------------------------------

def rot_about_pin(p, y0, z0, ang):
    """Rotate a point about the hinge axis by `ang` radians."""
    x, y, z = p
    dy, dz = y - y0, z - z0
    c, s = cos(ang), sin(ang)
    return (x, y0 + dy * c - dz * s, z0 + dy * s + dz * c)


# --------------------------------------------------------------------------
# marching tetrahedra
# --------------------------------------------------------------------------

# Each cube is split into 6 tetrahedra sharing the main diagonal 0-6.
_TETS = ((0, 5, 1, 6), (0, 1, 2, 6), (0, 2, 3, 6),
         (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6))

_CORNER = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
           (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))


class Grid:
    """Sampled scalar field on a regular lattice."""

    def __init__(self, fn, lo, hi, step):
        self.lo, self.step = lo, step
        self.n = tuple(int((hi[i] - lo[i]) / step) + 3 for i in range(3))
        nx, ny, nz = self.n
        vals = array('f', [0.0]) * (nx * ny * nz)
        i = 0
        for k in range(nz):
            z = lo[2] + (k - 1) * step
            for j in range(ny):
                y = lo[1] + (j - 1) * step
                base = (k * ny + j) * nx
                for ii in range(nx):
                    vals[base + ii] = fn((lo[0] + (ii - 1) * step, y, z))
                i += 1
        self.v = vals

    def at(self, i, j, k):
        nx, ny, _ = self.n
        return self.v[(k * ny + j) * nx + i]

    def pos(self, i, j, k):
        s = self.step
        return (self.lo[0] + (i - 1) * s,
                self.lo[1] + (j - 1) * s,
                self.lo[2] + (k - 1) * s)


def polygonise(grid):
    """Marching tetrahedra -> (verts, tris).

    Vertices are keyed by the lattice edge they sit on, so shared edges yield
    one shared vertex and the output is watertight by construction.
    """
    nx, ny, nz = grid.n
    at, pos = grid.at, grid.pos
    verts, vindex, tris = [], {}, []

    def vid(a, b):
        key = (a, b) if a < b else (b, a)
        got = vindex.get(key)
        if got is not None:
            return got
        (i0, j0, k0), (i1, j1, k1) = a, b
        d0, d1 = at(i0, j0, k0), at(i1, j1, k1)
        t = 0.5 if d1 == d0 else d0 / (d0 - d1)
        p0, p1 = pos(i0, j0, k0), pos(i1, j1, k1)
        v = (p0[0] + (p1[0] - p0[0]) * t,
             p0[1] + (p1[1] - p0[1]) * t,
             p0[2] + (p1[2] - p0[2]) * t)
        n = len(verts)
        verts.append(v)
        vindex[key] = n
        return n

    for k in range(nz - 1):
        for j in range(ny - 1):
            for i in range(nx - 1):
                c = [(i + dx, j + dy, k + dz) for dx, dy, dz in _CORNER]
                d = [at(*q) for q in c]
                if min(d) > 0.0 or max(d) <= 0.0:
                    continue
                for t in _TETS:
                    p = [c[m] for m in t]
                    dv = [d[m] for m in t]
                    inside = [m for m in range(4) if dv[m] <= 0.0]
                    ni = len(inside)
                    if ni == 0 or ni == 4:
                        continue
                    out = [m for m in range(4) if dv[m] > 0.0]
                    if ni == 1:
                        a = inside[0]
                        tri = [vid(p[a], p[o]) for o in out]
                        tris.append((tri[0], tri[1], tri[2]))
                    elif ni == 3:
                        o = out[0]
                        tri = [vid(p[a], p[o]) for a in inside]
                        tris.append((tri[0], tri[2], tri[1]))
                    else:
                        a, b = inside
                        c1, c2 = out
                        v00, v01 = vid(p[a], p[c1]), vid(p[a], p[c2])
                        v10, v11 = vid(p[b], p[c1]), vid(p[b], p[c2])
                        tris.append((v00, v01, v11))
                        tris.append((v00, v11, v10))
    return verts, tris


# --------------------------------------------------------------------------
# mesh utilities
# --------------------------------------------------------------------------

def orient(verts, tris, fn, eps):
    """Flip any triangle whose normal points into the solid.

    Degenerate triangles are kept, not dropped.  Marching tetrahedra emits a
    fair number of zero-area triangles where an interpolated vertex lands on
    a lattice point, and discarding them punches holes in an otherwise
    watertight mesh -- they carry real connectivity even with no area.
    """
    out = []
    for a, b, c in tris:
        pa, pb, pc = verts[a], verts[b], verts[c]
        ux, uy, uz = pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]
        vx, vy, vz = pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]
        nx_, ny_, nz_ = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
        L = sqrt(nx_ * nx_ + ny_ * ny_ + nz_ * nz_)
        if L == 0.0:
            out.append((a, b, c))
            continue
        cx = (pa[0] + pb[0] + pc[0]) / 3.0
        cy = (pa[1] + pb[1] + pc[1]) / 3.0
        cz = (pa[2] + pb[2] + pc[2]) / 3.0
        s = eps / L
        if fn((cx + nx_ * s, cy + ny_ * s, cz + nz_ * s)) < 0.0:
            out.append((a, c, b))
        else:
            out.append((a, b, c))
    return out


def project_to_surface(p, fn, h=0.01, steps=2):
    """Newton-step a point onto the fn = 0 isosurface along the gradient."""
    for _ in range(steps):
        d = fn(p)
        gx = (fn((p[0] + h, p[1], p[2])) - d) / h
        gy = (fn((p[0], p[1] + h, p[2])) - d) / h
        gz = (fn((p[0], p[1], p[2] + h)) - d) / h
        g2 = gx * gx + gy * gy + gz * gz
        if g2 < 1e-12:
            break
        k = d / g2
        p = (p[0] - gx * k, p[1] - gy * k, p[2] - gz * k)
    return p


def smooth_project(verts, tris, fn, passes=3, weight=0.65, crease_deg=42.0):
    """Tidy the triangulation without moving the surface.

    Marching tetrahedra puts vertices wherever the lattice edges happen to
    cross, which leaves slivers and an uneven distribution -- that, not the
    grid pitch, is what makes a curved surface look faceted under shading.
    This relaxes each vertex toward its neighbours and then pushes it back
    onto the exact isosurface, so the triangles get better shaped while every
    vertex still lies on the true surface.

    Vertices on a crease are held still, so the split line, the latch slot
    and the flats stay sharp instead of being rounded off.
    """
    n = len(verts)
    nbr = [set() for _ in range(n)]
    for a, b, c in tris:
        nbr[a].update((b, c))
        nbr[b].update((a, c))
        nbr[c].update((a, b))

    # Face normals gathered per vertex, to find creases.  Weighted by area:
    # marching tetrahedra emits a lot of slivers whose normals are numerical
    # noise, and an unweighted average lets them masquerade as creases -- on
    # the first attempt that flagged 85% of the surface and the smoothing did
    # essentially nothing.
    acc = [[0.0, 0.0, 0.0] for _ in range(n)]
    warea = [0.0] * n
    faces = [[] for _ in range(n)]
    fnorm, farea = [], []
    for a, b, c in tris:
        pa, pb, pc = verts[a], verts[b], verts[c]
        ux, uy, uz = pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]
        vx, vy, vz = pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]
        nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
        L = sqrt(nx * nx + ny * ny + nz * nz)
        i = len(fnorm)
        if L == 0.0:
            fnorm.append(None)
            farea.append(0.0)
            continue
        fnorm.append((nx / L, ny / L, nz / L))
        farea.append(0.5 * L)
        for v in (a, b, c):
            faces[v].append(i)
            acc[v][0] += nx / L * L
            acc[v][1] += ny / L * L
            acc[v][2] += nz / L * L
            warea[v] += 0.5 * L

    lim = cos(crease_deg * 3.141592653589793 / 180.0)
    free = [True] * n
    for v in range(n):
        ax, ay, az = acc[v]
        L = sqrt(ax * ax + ay * ay + az * az)
        if L == 0.0 or warea[v] <= 0.0:
            free[v] = False
            continue
        ax, ay, az = ax / L, ay / L, az / L
        big = 0.10 * warea[v]        # ignore slivers when judging a crease
        for i in faces[v]:
            f = fnorm[i]
            if f is None or farea[i] < big:
                continue
            if ax * f[0] + ay * f[1] + az * f[2] < lim:
                free[v] = False
                break

    out = list(verts)
    for _ in range(passes):
        moved = list(out)
        for v in range(n):
            if not free[v] or not nbr[v]:
                continue
            sx = sy = sz = 0.0
            for u in nbr[v]:
                q = out[u]
                sx += q[0]; sy += q[1]; sz += q[2]
            k = len(nbr[v])
            p = out[v]
            moved[v] = (p[0] + weight * (sx / k - p[0]),
                        p[1] + weight * (sy / k - p[1]),
                        p[2] + weight * (sz / k - p[2]))
        out = [project_to_surface(moved[v], fn, steps=1) if free[v] else moved[v]
               for v in range(n)]
    return [project_to_surface(out[v], fn, steps=1) if free[v] else out[v]
            for v in range(n)], sum(1 for f in free if f)


def tri_area(pa, pb, pc):
    ux, uy, uz = pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]
    vx, vy, vz = pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]
    nx_, ny_, nz_ = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
    return 0.5 * sqrt(nx_ * nx_ + ny_ * ny_ + nz_ * nz_)


def volume(verts, tris):
    """Signed volume via the divergence theorem (mm^3)."""
    tot = 0.0
    for a, b, c in tris:
        pa, pb, pc = verts[a], verts[b], verts[c]
        tot += (pa[0] * (pb[1] * pc[2] - pb[2] * pc[1])
                - pa[1] * (pb[0] * pc[2] - pb[2] * pc[0])
                + pa[2] * (pb[0] * pc[1] - pb[1] * pc[0]))
    return tot / 6.0


def surface_points(verts, tris, want):
    """Roughly `want` points spread over the mesh, area-weighted."""
    areas = [tri_area(verts[a], verts[b], verts[c]) for a, b, c in tris]
    total = sum(areas)
    if total <= 0.0:
        return []
    pts, acc = [], 0.0
    per = total / want
    # deterministic low-discrepancy barycentric samples
    golden = 0.6180339887498949
    seq = 0.0
    for (a, b, c), ar in zip(tris, areas):
        acc += ar
        n = int(acc / per)
        acc -= n * per
        pa, pb, pc = verts[a], verts[b], verts[c]
        for _ in range(n):
            seq = (seq + golden) % 1.0
            s = seq
            seq = (seq + golden) % 1.0
            t = seq
            if s + t > 1.0:
                s, t = 1.0 - s, 1.0 - t
            u = 1.0 - s - t
            pts.append((pa[0] * u + pb[0] * s + pc[0] * t,
                        pa[1] * u + pb[1] * s + pc[1] * t,
                        pa[2] * u + pb[2] * s + pc[2] * t))
    return pts


def check_manifold(verts, tris):
    """Return (is_watertight, is_single_shell, n_bad_edges, n_components)."""
    edge = {}
    for a, b, c in tris:
        for u, v in ((a, b), (b, c), (c, a)):
            key = (u, v) if u < v else (v, u)
            edge[key] = edge.get(key, 0) + 1
    bad = sum(1 for n in edge.values() if n != 2)

    adj = {}
    for a, b, c in tris:
        for u, v in ((a, b), (b, c), (c, a)):
            adj.setdefault(u, set()).add(v)
            adj.setdefault(v, set()).add(u)
    seen, comps = set(), 0
    for s in adj:
        if s in seen:
            continue
        comps += 1
        stack = [s]
        seen.add(s)
        while stack:
            q = stack.pop()
            for r in adj[q]:
                if r not in seen:
                    seen.add(r)
                    stack.append(r)
    return bad == 0, comps == 1, bad, comps


# --------------------------------------------------------------------------
# STL
# --------------------------------------------------------------------------

def write_stl(path, meshes):
    """Binary STL.  `meshes` is a list of (verts, tris)."""
    import struct
    n = sum(len(t) for _, t in meshes)
    with open(path, 'wb') as f:
        f.write(b'strip case - print in place flip top'.ljust(80, b' '))
        f.write(struct.pack('<I', n))
        for verts, tris in meshes:
            for a, b, c in tris:
                pa, pb, pc = verts[a], verts[b], verts[c]
                ux, uy, uz = pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]
                vx, vy, vz = pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]
                nx_ = uy * vz - uz * vy
                ny_ = uz * vx - ux * vz
                nz_ = ux * vy - uy * vx
                L = sqrt(nx_ * nx_ + ny_ * ny_ + nz_ * nz_) or 1.0
                f.write(struct.pack('<12fH',
                                    nx_ / L, ny_ / L, nz_ / L,
                                    pa[0], pa[1], pa[2],
                                    pb[0], pb[1], pb[2],
                                    pc[0], pc[1], pc[2], 0))
    return n
