"""Flat-oval strip case -- print-in-place flip top.

Coordinate system
-----------------
    Z   long axis, Z = 0 is the print bed; the case stands upright, closed
    X   width axis     (the wide direction of the oval, 21.1 mm)
    Y   thickness axis (the narrow direction, 12.4 mm).  -Y is the BACK,
        where the hinge lives; +Y is the FRONT, where the latch lives.

Every cross section normal to Z is a stadium -- two semicircles joined by
flats -- and that whole cross section is SCALED to zero over each end, along
a circular profile whose radius is the body's half width.

That last part is the difference between a lip balm stick and a vape.  If the
ends are rounded by sweeping the stadium (offsetting it in 3D), the bottom
comes out flat across the middle with rounded corners, and it reads as a
tube.  Scaling the section instead makes the front view of each end a true
semicircle spanning the full width -- a U, not a flat with corners.

The hinge axis is parallel to X.  A positive rotation about it lifts the
front of the cap.  That sign is the reason only the *rear* of the cap needs
relief: every point in front of the pin rises on opening and clears
immediately, while every point behind it sinks and must be given a
cylindrical underside concentric with the pin.
"""

from math import hypot, sqrt, radians
from geom import (smin, smax, sd_stadium_2d, sd_teardrop,
                  sd_box, sd_cyl_x, rot_about_pin)

# ==========================================================================
# PARAMETERS
# ==========================================================================

# ---- contents -----------------------------------------------------------
STRIP_L, STRIP_W, STRIP_T = 60.0, 8.0, 0.5     # 0.5 mm confirmed by the user
STRIP_N = 10
BUNDLE_W, BUNDLE_T = STRIP_W, STRIP_N * STRIP_T          # 8.0 x 5.0 mm

# ---- outer envelope -----------------------------------------------------
OUT_T = 12.4                       # thickness (Y)
OUT_W = 21.1                       # width (X)      -> 1.702 : 1
H_TOTAL = 96.0
R_OUT = OUT_T / 2.0                # 6.20  stadium end radius
HF_OUT = OUT_W / 2.0 - R_OUT       # 4.35  stadium half flat
HALF_W = OUT_W / 2.0               # 10.55

# Each end is a half ellipsoid: semi-axis HALF_W across the width, R_OUT
# across the thickness, and BOT_R / TOP_R up the long axis.  Setting those
# equal to HALF_W makes the front view of each end an exact semicircle.
BOT_R = HALF_W                     # 10.55
TOP_R = HALF_W                     # 10.55

# How much is cut off the bottom tip so the part can stand on the bed.
#   0.00 -> perfectly round, touches the bed at a point, NEEDS SUPPORT
#   3.09 -> the roundest support-free bottom (45 deg at the bed), 103 mm^2
BOT_FLAT = 0.0

Z_BOT_C = BOT_R - BOT_FLAT         # centre plane of the bottom dome
Z_TOP_C = H_TOTAL - TOP_R          # 84.45 centre plane of the top dome

# ---- cavity -------------------------------------------------------------
CAV_W, CAV_T = 10.8, 7.0
R_CAV = CAV_T / 2.0                # 3.50
HF_CAV = CAV_W / 2.0 - R_CAV       # 1.90
Z_FLOOR = 11.5
CAV_FILLET = 0.6                   # blends the floor into the walls

# ---- split, neck, skirt -------------------------------------------------
Z_SPLIT = 76.35
GAP = 0.3                          # print-in-place clearance everywhere
Z_BODY_TOP = Z_SPLIT - GAP / 2.0   # 75.20
Z_CAP_BOT = Z_SPLIT + GAP / 2.0    # 75.50
SKIRT_T = 2.0                      # cap skirt over the neck
NECK_INSET = SKIRT_T + GAP         # 2.30
R_NECK = R_OUT - NECK_INSET        # 3.90
Z_NECK_TOP = 78.00                 # 1.80 tall, 1.50 of skirt engagement
CAP_WALL = 2.2                     # cap wall above the neck

Z_CAV_TOP = Z_NECK_TOP             # strips come out through the neck
CAV_DEPTH = Z_CAV_TOP - Z_FLOOR    # 66.00

# The cavity's FRONT wall tapers inward near the top.  That is what buys the
# neck enough wall to carry a 2.0 mm skirt, which in turn gives the latch
# groove 1.45 mm of cap wall behind it.  The opening still measures 6.25 mm
# back to front against a 5.0 mm bundle, so the strips clear it.
CAV_FRONT_LO, CAV_FRONT_HI = 69.0, Z_BODY_TOP
CAV_FRONT_TOP = 2.75

# ---- hinge --------------------------------------------------------------
BORE_D, PIN_D = 3.5, 2.7
R_BORE, R_PIN = BORE_D / 2.0, PIN_D / 2.0
EAR_MIN_WALL = 1.2                 # required material all round the bore
EAR_R = R_BORE + EAR_MIN_WALL      # 2.95  ear outline radius about the pin
EAR_X, EAR_T = 6.0, 2.4            # ear centre |x| and thickness
EAR_HALF = EAR_T / 2.0

# The pin sits as far back as the oval allows at the ears' outermost x:
# at x = EAR_X + EAR_HALF the skin is only +-5.51 mm from the centre plane,
# so pushing the pin further back would break the bore out of the ear -- the
# exact failure the brief calls out.
Y_PIN = -2.56
Z_PIN = 78.95
# ear outline extended upward so the bore's 45 degree teardrop roof still
# leaves EAR_MIN_WALL above it
EAR_EXT = R_BORE * 1.41421 + EAR_MIN_WALL - EAR_R          # 0.725
EAR_BOT_Z = 72.0                   # the ear roots this far down into the body

# The cap's relief where the ears pass through: a cylinder on the pin axis,
# sized to the ear's corner at the cap's own bottom edge, plus 0.4.
R_EAR_RELIEF = hypot(EAR_R, Z_PIN - Z_CAP_BOT) + 0.4      # 4.24

# The knuckle underside is a cylinder of exactly this radius about the axis.
# Flushness at the split would want 4.39 -- the distance from the pin to the
# cap's rear-bottom corner -- but the ears root through that same annulus,
# and a knuckle that large sweeps into them around 95 deg.  The ears win:
# the knuckle is cut back inside the ear relief, which costs a 1.4 mm gap
# across the back of the split line and reads as an ordinary hinge gap.
R_KNUCK = min(hypot(R_OUT + Y_PIN, Z_CAP_BOT - Z_PIN),
              R_EAR_RELIEF - 0.4)                         # 3.84
R_DISH = R_KNUCK + GAP                                    # 4.14

# Cylindrical reliefs get the hinge most of the way there, but the body's
# rear shoulder and the ears' lower corners still sweep into the cap.  Rather
# than hand-deriving another radius, the body subtracts the cap's actual
# swept envelope over the opening range -- correct by construction.
SWEEP_MAX_DEG = 105.0
SWEEP_STEP_DEG = 1.5
SWEEP_CLEAR = 0.28                 # just under GAP, so it cannot nibble the
                                   # neck, which is already built at GAP
# The gate reaches forward of the pin, because cap material that starts
# behind the hinge swings forward as well as down.  It cannot reach much
# further than this: the ears root through the band just ahead of the pin,
# and sweeping there cuts them off the body entirely.  The last bit of
# forward reach is handed to dish_neck instead, which runs on the same axis
# but is applied before the ears rather than after them.
SWEEP_Y_MAX = Y_PIN + 2.5
CAP_IN_HALF_W = OUT_W / 2.0 - CAP_WALL                    # 8.35
PIN_X_IN = EAR_X - EAR_HALF - 0.6  # pin reaches past the inner ear face

# ---- latch --------------------------------------------------------------
TAB_W = 8.0                        # cantilever width (X)
TAB_T = 1.8                        # cantilever thickness (Y)
TAB_ROOT_Z = 66.35
TAB_TOP_Z = Z_NECK_TOP
TAB_L = TAB_TOP_Z - TAB_ROOT_Z     # 11.50
# The bead is a triangular ridge with 45 degree faces, and the groove that
# receives it is the same ridge grown by GROOVE_GAP.  Both print clean, and
# unlike a round bead in a shallow V they actually fit each other.
BEAD_R = 0.75                      # bead proud of the tab face
BEAD_ENGAGE = BEAD_R - GAP         # 0.45 -- the deflection needed to open
BEAD_Z = 77.25
BEAD_W = 6.0
SLOT_W = 0.4                       # slot freeing the cantilever
GROOVE_GAP = 0.15

# ---- keyring ------------------------------------------------------------
RING_D = 7.0
R_RING = RING_D / 2.0
Z_RING = 5.8                       # 2.30 mm of material below the hole

# The cap interior tents shut at 45 degrees below this height.
Z_APEX = 93.5

# ---- material -----------------------------------------------------------
PLA_DENSITY = 1.24e-3              # g/mm^3
PLA_E = 3000.0                     # MPa, typical printed PLA


# ==========================================================================
# shared surfaces
# ==========================================================================

def profile_scale(z):
    """How much of the full cross section survives at this height.

    1.0 through the straight middle, falling to 0 at each tip along a
    circular profile.  Because the WHOLE section scales -- flats and radius
    together -- the front view of each end is a semicircle of radius BOT_R /
    TOP_R spanning the full width.  That is the U shape; sweeping the section
    instead would leave a flat with rounded corners.
    """
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
    """The one continuous surface, shared by body and cap.

    `grow` offsets it inward, which is how the cap's wall and the neck are
    derived -- they stay parallel to the skin through the domes, where a
    fixed cross section would not.
    """
    x, y, z = p
    s = profile_scale(z)
    if s <= 1e-9:
        zc = (Z_BOT_C - BOT_R) if z < Z_BOT_C else (Z_TOP_C + TOP_R)
        return hypot(hypot(x, y), z - zc) + grow
    return sd_stadium_2d(x, y, HF_OUT * s, R_OUT * s) + grow


def spine_dist(x, y, s=1.0):
    """Distance from the stadium spine -- the 'radius' of the oval."""
    return sd_stadium_2d(x, y, HF_OUT * s, 0.0)


def cav_front(z):
    """Front wall of the cavity: tapers inward near the top at about 15 deg
    from vertical, well inside the 45 degree printing limit."""
    if z <= CAV_FRONT_LO:
        return R_CAV
    if z >= CAV_FRONT_HI:
        return CAV_FRONT_TOP
    t = (z - CAV_FRONT_LO) / (CAV_FRONT_HI - CAV_FRONT_LO)
    return R_CAV + (CAV_FRONT_TOP - R_CAV) * t


def tab_face(z):
    """Outer face of the latch cantilever: flush with the skin below the
    split, flush with the neck above it.  The step between the two is the
    body's top rim, which faces up -- not an overhang."""
    return R_NECK if z >= Z_BODY_TOP else R_OUT


def tab_inner(z):
    """Inner face of the cantilever.  The tab has to step inward at the
    split because the neck is inset; ramping that step at 45 degrees is what
    keeps its underside from being a flat ceiling over the pocket."""
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
    """Strip bore: stadium prism, filleted into its floor, front wall tapered
    near the top, open through the neck."""
    walls = sd_stadium_2d(p[0], p[1], HF_CAV, R_CAV)
    walls = max(walls, p[1] - cav_front(p[2]))
    d = smax(walls, Z_FLOOR - p[2], CAV_FILLET)
    return max(d, p[2] - (Z_CAV_TOP + 5.0))


def keyring(p):
    """Ø4.5 through the bottom dome, axis along Y, teardrop roofed so the
    hole prints without support."""
    return max(sd_teardrop(p[0], p[2] - Z_RING, R_RING), abs(p[1]) - 20.0)


def latch_slot(p):
    """The U-slot that frees the cantilever on three sides, plus the pocket
    that thins it to TAB_T.

    There is deliberately no slot under the root: the tab is a cantilever
    built in at TAB_ROOT_Z, and cutting across there would set it loose as a
    separate floating part.
    """
    x, y, z = p
    hw = TAB_W / 2.0
    zc = (TAB_ROOT_Z + TAB_TOP_Z + 6.0) / 2.0
    zh = (TAB_TOP_Z + 6.0 - TAB_ROOT_Z) / 2.0
    side = min(sd_box(p, hw + SLOT_W / 2.0, 5.0, zc, SLOT_W / 2.0, 3.0, zh),
               sd_box(p, -hw - SLOT_W / 2.0, 5.0, zc, SLOT_W / 2.0, 3.0, zh))

    # Pocket behind the tab: from the cavity's front wall out to the tab's
    # inner face.  Above the taper the cavity itself has moved in far enough
    # that the pocket closes to nothing on its own -- so the void behind the
    # tab has no flat ceiling anywhere.
    y_in = cav_front(z) - 0.3
    y_out = tab_inner(z)
    pocket = max(max(abs(x) - hw, max(y_in - y, y - y_out)),
                 max(TAB_ROOT_Z - z, z - (TAB_TOP_Z + 6.0)))
    return min(side, pocket)


def ear_profile(y, z):
    """Ear outline in the plane normal to the pin: a stadium about a vertical
    segment through the pin axis.  Rounded, so nothing square swings near the
    cap; rooted below the split, so it never floats."""
    if z > Z_PIN + EAR_EXT:
        dz = z - (Z_PIN + EAR_EXT)
    elif z < EAR_BOT_Z:
        dz = z - EAR_BOT_Z
    else:
        dz = 0.0
    return hypot(y - Y_PIN, dz) - EAR_R


# ==========================================================================
# BODY
# ==========================================================================

def body(p):
    x, y, z = p
    r_pin = hypot(y - Y_PIN, z - Z_PIN)

    # capsule below the split, plus the neck collar above it
    d = max(outer(p), z - Z_BODY_TOP)
    neck = max(outer(p, NECK_INSET),
               max(Z_BODY_TOP - 0.5 - z, z - Z_NECK_TOP))
    d = smin(d, neck, 0.4)

    d = max(d, -z)                          # flat bottom, first layer

    # Two reliefs on the same axis, both R_DISH, but over different regions.
    #
    # Below the split it is the dish that receives the knuckle, and it is
    # confined to behind the pin: unbounded, it scallops away the body's top
    # rim all the way round and leaves the cap's rim overhanging nothing.
    #
    # At and above the split it applies all the way round, because that is
    # the neck, and the skirt sweeps over the whole of it.  Without this the
    # skirt's bottom edge digs into the neck's sides from about 30 deg -- the
    # 'deep skirt cannot rotate off a neck' failure.  What survives is a
    # front crescent, which is exactly where the latch needs it.
    dish_rear = max(r_pin - R_DISH, y - (Y_PIN + 0.5))
    # reaches 1.2 mm below the split so it also takes the top rim corner
    # that the cap's rear-bottom corner swings onto around 95 deg
    dish_neck = max(r_pin - R_DISH, (Z_BODY_TOP - 1.2) - z)
    d = max(d, -min(dish_rear, dish_neck))

    d = max(d, -cavity(p))
    d = max(d, -latch_slot(p))
    d = max(d, -keyring(p))

    # ears, added after the dish so the dish cannot eat them
    prof = ear_profile(y, z)
    for sx in (-1.0, 1.0):
        d = smin(d, max(prof, abs(x - sx * EAR_X) - EAR_HALF), 0.8)

    # Remove whatever the cap sweeps through.  This is what actually clears
    # the body's rear shoulder and rounds off the ears' lower corners -- the
    # square corners that jam a print-in-place hinge.
    d = max(d, SWEEP_CLEAR - cap_swept(p))

    # detent bead: a 45 degree triangular ridge on the cantilever
    yf = tab_face(BEAD_Z)
    bead = max(max(abs(x) - BEAD_W / 2.0, yf - y),
               (y - yf) + abs(z - BEAD_Z) - BEAD_R)
    d = smin(d, bead, 0.2)

    # pin bores, teardrop roofed so they need no support
    tp = sd_teardrop(y - Y_PIN, z - Z_PIN, R_BORE)
    for sx in (-1.0, 1.0):
        lo = sx * EAR_X - EAR_HALF - 1.0
        hi = sx * EAR_X + EAR_HALF + 1.0
        lo, hi = min(lo, hi), max(lo, hi)
        d = max(d, -max(tp, max(lo - x, x - hi)))

    return d


# ==========================================================================
# CAP
# ==========================================================================

def cap_wall(z):
    """Cap wall thickness: the skirt is thinner, so it clears the neck with
    GAP to spare; above the neck it thickens to CAP_WALL over a 45 degree
    taper so the step is self supporting."""
    if z <= Z_NECK_TOP:
        return SKIRT_T                         # inner face lands at R_NECK+GAP
    if z >= Z_NECK_TOP + 1.0:
        return CAP_WALL
    t = z - Z_NECK_TOP
    return SKIRT_T * (1.0 - t) + CAP_WALL * t


def cap(p):
    x, y, z = p
    r_pin = hypot(y - Y_PIN, z - Z_PIN)

    d = max(outer(p), Z_CAP_BOT - z)

    # Hollow interior.  It is the outer surface offset inward, so the wall
    # stays parallel to the skin right through the dome; a 45 degree tent
    # closes it off so nothing bridges mid air.
    void = max(outer(p, cap_wall(z)),
               (spine_dist(x, y, profile_scale(z)) - (Z_APEX - z)) * 0.7071)
    void = max(void, Z_CAP_BOT - 4.0 - z)
    d = max(d, -void)

    # Knuckle.  Behind the pin and below it the cap keeps material only
    # within R_KNUCK of the axis, so its underside is a cylinder concentric
    # with the pin rather than a square corner.
    d = max(d, -max(z - Z_PIN, max(y - Y_PIN, R_KNUCK - r_pin)))

    # relief where the body ears pass through the cap: a cylinder on the pin
    # axis, teardrop roofed so its own ceiling is not an overhang
    rel = sd_teardrop(y - Y_PIN, z - Z_PIN, R_EAR_RELIEF)
    for sx in (-1.0, 1.0):
        d = max(d, -max(rel, abs(x - sx * EAR_X) - (EAR_HALF + GAP)))

    # captive pins: Ø2.7 with a 45 degree V underside so they self support
    for sx in (-1.0, 1.0):
        dy, dz = y - Y_PIN, z - Z_PIN
        prof = max(hypot(dy, dz) - R_PIN, (abs(dy) - R_PIN - dz) * 0.7071)
        lo = min(sx * PIN_X_IN, sx * CAP_IN_HALF_W)
        hi = max(sx * PIN_X_IN, sx * CAP_IN_HALF_W)
        d = min(d, max(prof, max(lo - x, x - hi)))

    # latch groove: the bead's own profile, grown by GROOVE_GAP and clipped
    # at the skirt's inner face.  Its roof is a 45 degree face.
    # There is no clip at the skirt's inner face: everything inboard of it is
    # already void, and clipping there left the bead touching the mouth of
    # its own groove with zero clearance.
    notch = max(abs(x) - (BEAD_W / 2.0 + 0.2),
                (y - R_NECK) + abs(z - BEAD_Z) - (BEAD_R + GROOVE_GAP))
    d = max(d, -notch)

    return d


def cap_nominal(p):
    """The cap without the knurl and latch groove.  Both only remove
    material, so leaving them out makes the swept envelope conservative."""
    x, y, z = p
    r_pin = hypot(y - Y_PIN, z - Z_PIN)
    d = max(outer(p), Z_CAP_BOT - z)
    void = max(outer(p, cap_wall(z)),
               (spine_dist(x, y, profile_scale(z)) - (Z_APEX - z)) * 0.7071)
    void = max(void, Z_CAP_BOT - 4.0 - z)
    d = max(d, -void)
    d = max(d, -max(z - Z_PIN, max(y - Y_PIN, R_KNUCK - r_pin)))
    rel = sd_teardrop(y - Y_PIN, z - Z_PIN, R_EAR_RELIEF)
    for sx in (-1.0, 1.0):
        d = max(d, -max(rel, abs(x - sx * EAR_X) - (EAR_HALF + GAP)))
    for sx in (-1.0, 1.0):
        dy, dz = y - Y_PIN, z - Z_PIN
        prof = max(hypot(dy, dz) - R_PIN, (abs(dy) - R_PIN - dz) * 0.7071)
        lo = min(sx * PIN_X_IN, sx * CAP_IN_HALF_W)
        hi = max(sx * PIN_X_IN, sx * CAP_IN_HALF_W)
        d = min(d, max(prof, max(lo - x, x - hi)))
    return d


_SWEEP_ANGLES = [radians(-a) for a in
                 [i * SWEEP_STEP_DEG for i in
                  range(int(SWEEP_MAX_DEG / SWEEP_STEP_DEG) + 1)]]


def cap_swept(p):
    """Closest approach of the cap to p over the whole opening range.

    Evaluated by rotating p backwards through each station, which is the same
    thing as rotating the cap forwards.  Gated to the rear of the hinge: the
    front of the cap rises away on opening and never needs body relief, and
    the gate keeps the latch bead out of its own groove's way.
    """
    if p[1] > SWEEP_Y_MAX or p[2] < 69.0 or p[2] > 87.0:
        return 1e3
    # Exact short circuit inside an ear band.  The cap's relief there is a
    # cylinder concentric with the pin, so it is rotation invariant: at any
    # angle the cap has material only outside R_EAR_RELIEF, or on the pin
    # stub inside R_PIN.  Both bounds ignore separation along x, so this
    # under-states the true distance and can only carve more, never less.
    # Without it every probe near the bore costs 71 cap evaluations.
    r_pin = hypot(p[1] - Y_PIN, p[2] - Z_PIN)
    if R_PIN < r_pin < R_EAR_RELIEF:
        for sx in (-1.0, 1.0):
            if abs(p[0] - sx * EAR_X) <= EAR_HALF + GAP:
                return min(r_pin - R_PIN, R_EAR_RELIEF - r_pin)
    best = 1e3
    for a in _SWEEP_ANGLES:
        d = cap_nominal(rot_about_pin(p, Y_PIN, Z_PIN, a))
        if d < best:
            best = d
    return best
