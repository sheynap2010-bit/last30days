"""Strip case -- print-in-place flip top, cylindrical, external hinge.

Coordinate system
-----------------
    Z   long axis, Z = 0 is the print bed; the case stands upright, closed
    X   across the hinge pin
    Y   -Y is the BACK, where the hinge lives; +Y is the FRONT, the latch

The body is a cylinder with hemispherical ends -- a capsule.  The cross
section machinery is the same one that drew the earlier flat oval: a stadium
scaled to zero over each end along a circular profile.  Setting the width
equal to the thickness collapses the stadium to a circle, so the same code
gives a true cylinder with domed ends and no special casing.

Why the hinge is OUTSIDE the skin
---------------------------------
Opening rotates the cap about the pin, and every point of the cap BEHIND the
pin sinks as it rotates.  With the pin buried in the back wall, the cap's
whole rear skin is behind it, so the body had to be hollowed out of its way
-- and that hollow was a 2 mm slot from the cavity to open air, wide enough
for a strip to escape.  Four fixes were built and measured against that and
all failed (see the README).

Putting the pin 0.15 mm behind the skin makes every point of the cap sit in
FRONT of it.  Everything rises on opening, nothing sweeps into the body, and
the split line closes to its 0.3 mm print clearance the whole way round.

That leaves only the cap's own lug, which hangs down behind the pin.  It is
bounded by a cylinder on the pin axis, so its swept envelope is that same
cylinder at every angle -- the body just carries a matching cylindrical
recess.  The recess is 3.0 mm deep into a 5.5 mm back wall, so it is a dish
in the outer surface, not a hole into the cavity.

Because of that, this model needs no swept-envelope subtraction, no dish and
no knuckle cut.  Every hinge clearance is an explicit cylinder on the axis.
"""

from math import hypot, sqrt
from geom import (smin, smax, sd_stadium_2d, sd_teardrop, sd_segment_2d,
                  sd_box, sd_cyl_x)

# ==========================================================================
# PARAMETERS
# ==========================================================================

# ---- contents -----------------------------------------------------------
STRIP_L, STRIP_W, STRIP_T = 60.0, 8.0, 0.5
STRIP_N = 10
BUNDLE_W, BUNDLE_T = STRIP_W, STRIP_N * STRIP_T          # 8.0 x 5.0 mm

# ---- outer envelope -----------------------------------------------------
# Width == thickness, so the stadium collapses to a circle: a cylinder.
OUT_T = 18.0
OUT_W = 18.0
H_TOTAL = 94.0
R_OUT = OUT_T / 2.0                # 9.00
HF_OUT = OUT_W / 2.0 - R_OUT       # 0.00 -- no flat, the section is a circle
HALF_W = OUT_W / 2.0               # 9.00

BOT_R = HALF_W                     # hemispherical bottom
TOP_R = HALF_W                     # hemispherical top

# How much is cut off the bottom tip so the part can stand on the bed.
#   0.00 -> a true hemisphere, touches at a point, NEEDS SUPPORT + BRIM
#   2.64 -> the roundest support-free bottom (45 deg at the bed)
BOT_FLAT = 0.0

Z_BOT_C = BOT_R - BOT_FLAT         # 9.00 centre plane of the bottom dome
Z_TOP_C = H_TOTAL - TOP_R          # 85.00 centre plane of the top dome

# ---- cavity -------------------------------------------------------------
CAV_W, CAV_T = 10.8, 7.0
R_CAV = CAV_T / 2.0                # 3.50
HF_CAV = CAV_W / 2.0 - R_CAV       # 1.90
Z_FLOOR = 11.5
CAV_FILLET = 0.6

# ---- split, neck, skirt -------------------------------------------------
Z_SPLIT = 76.35
GAP = 0.3                          # clearance across VERTICAL faces

# Clearance across HORIZONTAL faces has to be bigger, and this is the single
# thing most likely to weld the cap to the body.  A vertical gap is safe at
# 0.3 mm because the nozzle simply never enters it.  A horizontal gap is not:
# the layer above it is laid down over air, and at 0.3 mm with 0.2 mm layers
# that is 1.5 layers -- the cap's first layer lands almost on the body's top
# face and fuses to it.  0.6 mm is three layers at 0.20 and four at 0.15, so
# at least two are fully empty.
SPLIT_GAP = 0.6

Z_BODY_TOP = Z_SPLIT - GAP / 2.0   # 76.20
Z_CAP_BOT = Z_BODY_TOP + SPLIT_GAP # 76.80
SKIRT_T = 2.0
NECK_INSET = SKIRT_T + GAP         # 2.30
R_NECK = R_OUT - NECK_INSET        # 6.70  -> 1.30 mm of neck wall at the sides
Z_NECK_TOP = 78.00
CAP_WALL = 2.2
Z_APEX = 91.0                      # the cap's interior tents shut here

Z_CAV_TOP = Z_NECK_TOP
CAV_DEPTH = Z_CAV_TOP - Z_FLOOR    # 66.50

# The Ø18 barrel leaves 5.5 mm of wall front and back, so unlike the flat
# oval the cavity needs no taper to make room for the neck.
CAV_FRONT_TOP = R_CAV

# ---- hinge (external) ---------------------------------------------------
BORE_D, PIN_D = 3.5, 2.7
R_BORE, R_PIN = BORE_D / 2.0, PIN_D / 2.0
EAR_MIN_WALL = 1.2
EAR_R = R_BORE + EAR_MIN_WALL      # 2.95

Y_PIN = -(R_OUT + 0.15)            # -9.15, just behind the skin
Z_PIN = 72.50                      # low enough that the ears clear the cap
EAR_EXT = R_BORE * 1.41421 + EAR_MIN_WALL - EAR_R          # 0.725

# The ear leans forward and down into the body, at 45 degrees.  That roots it
# properly AND means its protrusion fades out downward at 45 degrees instead
# of hanging off the back as an overhang.
# 2.6, not 4.0: leaning 4 mm forward put the ear's root inside the cavity at
# x = +-4, right where the strip bundle's corners sit.  The limit is set by
# the cavity's back wall at the ear's inner edge (y = -3.21 at x = 3.3).
EAR_LEAN, EAR_DROP = 2.6, 2.6
EAR_X, EAR_T = 4.5, 2.4            # ear centre |x| and thickness
EAR_HALF = EAR_T / 2.0

# The cap's lug, and the recess in the body that receives it.  The lug is
# bounded by a cylinder on the pin axis, so it sweeps that same cylinder at
# every opening angle -- the recess is simply that cylinder plus clearance.
R_LUG = R_PIN + EAR_MIN_WALL       # 2.55
LUG_HALF_X = 2.8
RECESS_R = R_LUG + 0.5             # 3.05 -- the lug's underside is a
                                   # horizontal interface, so it gets the
                                   # larger clearance too
RECESS_HALF_X = LUG_HALF_X + GAP   # 3.10
LUG_BACK = Y_PIN - R_LUG           # -11.70, the boss's outermost face
# The strap must clear the body's skin where it passes it, but ABOVE the body
# it has to run forward far enough to actually join the cap's wall.  Held at
# the clearance value all the way up, it floats free and the cap comes out as
# two disconnected solids -- which a coarse marching grid hides, because the
# grid bridges the 0.3 mm gap.
STRAP_FRONT_LO = -(R_OUT + GAP)      # -9.30, alongside the body
STRAP_FRONT_HI = -(R_OUT - SKIRT_T)  # -7.00, merged into the cap's wall
STRAP_RAMP_TOP = 78.8                # 41 deg ramp, inside the print limit
PIN_X0, PIN_X1 = 2.4, 6.2          # the captive pin spans this, mirrored

# ---- latch --------------------------------------------------------------
TAB_W = 8.0
TAB_T = 1.8
TAB_ROOT_Z = 66.35
TAB_TOP_Z = Z_NECK_TOP
TAB_L = TAB_TOP_Z - TAB_ROOT_Z     # 11.65
BEAD_R = 0.75
BEAD_ENGAGE = BEAD_R - GAP         # 0.45
BEAD_Z = 77.25
BEAD_W = 6.0
# The slot has to stay UNDER the strip thickness, because it runs from the
# cap's interior down past the tab and out through the skin -- at 0.6 mm a
# 0.5 mm strip walks straight out of it.  0.45 blocks the strip and is a
# vertical gap, which is far more forgiving to print than a horizontal one:
# the nozzle never passes over it.  Disable gap-fill in the slicer so it does
# not try to bridge it with a thin extrusion.
SLOT_W = 0.45
# The slot and the pocket behind the tab both stop 1.0 mm short of the cavity
# wall, leaving a membrane between them.  Without it the slot is an open
# channel from the cavity to the outside -- and at 0.6 mm it is wide enough
# for a 0.5 mm strip to walk out through.
POCKET_Y_IN = R_CAV + 1.0          # 4.50
GROOVE_GAP = 0.15

# ---- keyring ------------------------------------------------------------
RING_D = 7.0
R_RING = RING_D / 2.0
Z_RING = 5.8                       # 2.30 below, 4.91 either side

# ---- material -----------------------------------------------------------
PLA_DENSITY = 1.24e-3
PLA_E = 3000.0


# ==========================================================================
# shared surfaces
# ==========================================================================

def profile_scale(z):
    """How much of the full cross section survives at this height: 1.0 down
    the barrel, falling to 0 at each tip along a circular profile."""
    if z < Z_BOT_C:
        t = (Z_BOT_C - z) / BOT_R
    elif z > Z_TOP_C:
        t = (z - Z_TOP_C) / TOP_R
    else:
        return 1.0
    if t >= 1.0:
        return 0.0
    return sqrt(1.0 - t * t)


def outer(p, grow=0.0):
    """The one continuous surface, shared by body and cap.  `grow` offsets it
    inward, which is how the neck and the cap's wall are derived."""
    x, y, z = p
    s = profile_scale(z)
    if s <= 1e-9:
        zc = (Z_BOT_C - BOT_R) if z < Z_BOT_C else (Z_TOP_C + TOP_R)
        return hypot(hypot(x, y), z - zc) + grow
    return sd_stadium_2d(x, y, HF_OUT * s, R_OUT * s) + grow


def spine_dist(x, y, s=1.0):
    return sd_stadium_2d(x, y, HF_OUT * s, 0.0)


def cav_front(z):
    return CAV_FRONT_TOP


def tab_face(z):
    """Outer face of the latch cantilever: flush with the skin below the
    split, flush with the neck above it."""
    return R_NECK if z >= Z_BODY_TOP else R_OUT


def tab_inner(z):
    """Inner face of the cantilever.  The step at the split is ramped at 45
    degrees so its underside is not a flat ceiling over the pocket."""
    hi, lo = R_NECK - TAB_T, R_OUT - TAB_T
    z0 = Z_BODY_TOP - (lo - hi)
    if z >= Z_BODY_TOP:
        return hi
    if z <= z0:
        return lo
    return lo + (hi - lo) * (z - z0) / (Z_BODY_TOP - z0)


# ==========================================================================
# features
# ==========================================================================

def cavity(p):
    """Strip bore: a stadium prism, filleted into its floor, open at the top
    through the neck."""
    walls = sd_stadium_2d(p[0], p[1], HF_CAV, R_CAV)
    walls = max(walls, p[1] - cav_front(p[2]))
    d = smax(walls, Z_FLOOR - p[2], CAV_FILLET)
    return max(d, p[2] - (Z_CAV_TOP + 5.0))


def keyring(p):
    """Ø7 through the bottom dome, axis along Y, teardrop roofed."""
    return max(sd_teardrop(p[0], p[2] - Z_RING, R_RING), abs(p[1]) - 20.0)


def latch_slot(p):
    """U-slot freeing the cantilever on three sides, plus the pocket that
    thins it.  No slot under the root -- that would set the tab loose."""
    x, y, z = p
    hw = TAB_W / 2.0
    zc = (TAB_ROOT_Z + TAB_TOP_Z + 6.0) / 2.0
    zh = (TAB_TOP_Z + 6.0 - TAB_ROOT_Z) / 2.0
    yc = (POCKET_Y_IN + 12.0) / 2.0
    yh = (12.0 - POCKET_Y_IN) / 2.0
    side = min(sd_box(p, hw + SLOT_W / 2.0, yc, zc, SLOT_W / 2.0, yh, zh),
               sd_box(p, -hw - SLOT_W / 2.0, yc, zc, SLOT_W / 2.0, yh, zh))
    y_out = tab_inner(z)
    pocket = max(max(abs(x) - hw, max(POCKET_Y_IN - y, y - y_out)),
                 max(TAB_ROOT_Z - z, z - (TAB_TOP_Z + 6.0)))
    return min(side, pocket)


def ear_profile(y, z):
    """Ear outline in the plane normal to the pin.

    Two stadia unioned: a short vertical one centred ON the pin, which
    guarantees EAR_R of material in every direction round the bore, and a
    leaning one that carries the ear forward and down into the body at 45
    degrees.  The lean alone is not enough -- its segment starts above the
    pin, so below and behind the pin the wall thinned to 0.73 mm.
    """
    up = sd_segment_2d(y, z, Y_PIN, Z_PIN, Y_PIN, Z_PIN + EAR_EXT) - EAR_R
    lean = sd_segment_2d(y, z, Y_PIN, Z_PIN,
                         Y_PIN + EAR_LEAN, Z_PIN - EAR_DROP) - EAR_R
    return min(up, lean)


def lug_recess(p):
    """The dish in the body's back that receives the cap's lug.

    A cylinder on the pin axis: the lug lives entirely within R_LUG of that
    axis, and radius is invariant under rotation, so this one cylinder clears
    the lug at every opening angle without any swept-envelope machinery.
    """
    x, y, z = p
    return max(hypot(y - Y_PIN, z - Z_PIN) - RECESS_R,
               abs(x) - RECESS_HALF_X)


# ==========================================================================
# BODY
# ==========================================================================

def body(p):
    x, y, z = p

    d = max(outer(p), z - Z_BODY_TOP)
    neck = max(outer(p, NECK_INSET),
               max(Z_BODY_TOP - 0.5 - z, z - Z_NECK_TOP))
    d = smin(d, neck, 0.4)

    d = max(d, -z)                          # flat bottom, if BOT_FLAT > 0

    d = max(d, -cavity(p))
    d = max(d, -latch_slot(p))
    d = max(d, -keyring(p))
    d = max(d, -lug_recess(p))

    # ears, added after the recess so it cannot eat them
    prof = ear_profile(y, z)
    for sx in (-1.0, 1.0):
        d = smin(d, max(prof, abs(x - sx * EAR_X) - EAR_HALF), 0.8)

    # detent bead: a 45 degree triangular ridge on the cantilever
    yf = tab_face(BEAD_Z)
    bead = max(max(abs(x) - BEAD_W / 2.0, yf - y),
               (y - yf) + abs(z - BEAD_Z) - BEAD_R)
    d = smin(d, bead, 0.2)

    # pin bores, teardrop roofed so they need no support
    tp = sd_teardrop(y - Y_PIN, z - Z_PIN, R_BORE)
    for sx in (-1.0, 1.0):
        lo, hi = sorted((sx * (PIN_X0 - 0.4), sx * (PIN_X1 + 0.4)))
        d = max(d, -max(tp, max(lo - x, x - hi)))

    return d


# ==========================================================================
# CAP
# ==========================================================================

def strap_front(z):
    """Front face of the cap's hinge strap: held clear of the body's skin
    while it passes the body, then ramped forward at 41 degrees to merge into
    the cap's wall."""
    if z <= Z_BODY_TOP:
        return STRAP_FRONT_LO
    if z >= STRAP_RAMP_TOP:
        return STRAP_FRONT_HI
    t = (z - Z_BODY_TOP) / (STRAP_RAMP_TOP - Z_BODY_TOP)
    return STRAP_FRONT_LO + (STRAP_FRONT_HI - STRAP_FRONT_LO) * t


def cap_wall(z):
    """Skirt is thinner so it clears the neck; above the neck it thickens to
    CAP_WALL over a 45 degree taper."""
    if z <= Z_NECK_TOP:
        return SKIRT_T
    if z >= Z_NECK_TOP + 1.0:
        return CAP_WALL
    t = z - Z_NECK_TOP
    return SKIRT_T * (1.0 - t) + CAP_WALL * t


def cap(p):
    x, y, z = p

    d = max(outer(p), Z_CAP_BOT - z)

    # hollow interior, offset from the skin so the wall follows the dome,
    # with a 45 degree tent so nothing bridges mid air
    void = max(outer(p, cap_wall(z)),
               (spine_dist(x, y, profile_scale(z)) - (Z_APEX - z)) * 0.7071)
    void = max(void, Z_CAP_BOT - 4.0 - z)
    d = max(d, -void)

    # Lug: the disc around the pin, plus a strap up to the cap.  The disc is
    # inside R_LUG of the axis so the recess always clears it; the strap sits
    # behind the body's skin, and being above the pin it only ever swings
    # further back as the cap opens.
    disc = max(hypot(y - Y_PIN, z - Z_PIN) - R_LUG, abs(x) - LUG_HALF_X)
    strap = max(max(abs(x) - LUG_HALF_X, y - strap_front(z)),
                max(max(LUG_BACK - y, Z_PIN - z), z - STRAP_RAMP_TOP))
    d = smin(d, min(disc, strap), 0.5)

    # captive pins: Ø2.7 with a 45 degree V underside so they self support
    for sx in (-1.0, 1.0):
        dy, dz = y - Y_PIN, z - Z_PIN
        prof = max(hypot(dy, dz) - R_PIN, (abs(dy) - R_PIN - dz) * 0.7071)
        lo, hi = sorted((sx * PIN_X0, sx * PIN_X1))
        d = min(d, max(prof, max(lo - x, x - hi)))

    # latch groove: the bead's profile grown by GROOVE_GAP, bounded in -Y so
    # it stays a groove instead of a channel through the whole cap
    notch = max(abs(x) - (BEAD_W / 2.0 + 0.2),
                max((R_NECK - 0.6) - y,
                    (y - R_NECK) + abs(z - BEAD_Z) - (BEAD_R + GROOVE_GAP)))
    d = max(d, -notch)

    return d
