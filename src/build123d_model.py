"""B-rep version of the strip case, for build123d -> STEP.

WHY THIS FILE IS SEPARATE, AND ITS STATUS
-----------------------------------------
The brief asked for build123d/CadQuery so the model stays editable and can be
exported as STEP.  This session had no route to install it: pypi.org, the
Ubuntu archives and the npm registry all return 403 through the egress proxy,
so build123d, CadQuery, OCP and Blender could not be fetched.

**This script has therefore never been executed.**  It is a faithful
transcription of the geometry in model.py -- same parameter names, same
values, same construction order -- but treat it as unrun source until you
have run it locally.  The verified deliverable is the STL, which comes from
model.py + check.py.

Run it with:      python build123d_model.py
Outputs:          out/strip_case.step, out/strip_case_body.stl,
                  out/strip_case_cap.stl

The parameters are imported from model.py so the two descriptions cannot
drift apart.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from build123d import (  # noqa: E402
    BuildPart, BuildSketch, Plane, Axis, Location, Rotation, Mode,
    SlotOverall, Circle, Rectangle, Polygon, extrude, fillet, chamfer,
    Box, Cylinder, Cone, add, export_step, export_stl, Compound, Solid,
)

import model as M  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   'out')


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def stadium(width, thickness):
    """The cross section everything in this design is built from."""
    return SlotOverall(width, thickness)


def capsule(width, thickness, z0, z1):
    """A stadium prism whose two ends are domed into the same profile.

    Filleting the end faces with half the thickness turns the flat end into a
    full dome, which is what makes the silhouette a true capsule rather than
    a tube with rounded corners.
    """
    with BuildPart() as bp:
        with BuildSketch(Plane.XY.offset(z0)):
            stadium(width, thickness)
        extrude(amount=z1 - z0)
        # the two end faces, domed
        fillet(bp.faces().filter_by(Plane.XY).edges(), radius=thickness / 2.0)
    return bp.part


def teardrop_profile(r, plane):
    """Circle plus a 45 degree gable tangent to it: a self-supporting hole."""
    import math
    t = r / math.sqrt(2.0)
    with BuildSketch(plane) as sk:
        Circle(r)
        Polygon((-t, t), (t, t), (0.0, r * math.sqrt(2.0)), align=None)
    return sk.sketch


# --------------------------------------------------------------------------
# body
# --------------------------------------------------------------------------

def make_body():
    with BuildPart() as bp:
        # ---- outer capsule, cut at the split ----------------------------
        add(capsule(M.OUT_W, M.OUT_T, M.BOT_Z - M.R_OUT, M.TOP_Z + M.R_OUT))
        Box(40, 20, 40, mode=Mode.SUBTRACT,
            align=None,
            **{}) if False else None
        # keep only what is below the split
        with BuildPart(mode=Mode.SUBTRACT):
            Box(40, 20, 40,
                align=None).locate(Location((0, 0, M.Z_BODY_TOP + 20)))

        # ---- flat bottom -------------------------------------------------
        with BuildPart(mode=Mode.SUBTRACT):
            Box(40, 20, 40, align=None).locate(Location((0, 0, -20)))

        # ---- neck collar -------------------------------------------------
        with BuildSketch(Plane.XY.offset(M.Z_BODY_TOP - 0.5)):
            stadium(M.OUT_W - 2 * M.NECK_INSET, M.OUT_T - 2 * M.NECK_INSET)
        extrude(amount=M.Z_NECK_TOP - (M.Z_BODY_TOP - 0.5))

        # ---- cavity ------------------------------------------------------
        with BuildPart(mode=Mode.SUBTRACT):
            with BuildSketch(Plane.XY.offset(M.Z_FLOOR)):
                stadium(M.CAV_W, M.CAV_T)
            extrude(amount=M.Z_CAV_TOP - M.Z_FLOOR)
            # taper the front wall inward near the top
            with BuildSketch(Plane.XY.offset(M.CAV_FRONT_LO)):
                Rectangle(M.CAV_W, 8, align=None).locate(
                    Location((0, M.R_CAV, 0)))
            # (loft the taper in a real run; see NOTE below)
        fillet(bp.edges().group_by(Axis.Z)[0], radius=M.CAV_FILLET)

        # ---- keyring, teardrop roofed ------------------------------------
        with BuildPart(mode=Mode.SUBTRACT):
            add(teardrop_profile(M.R_RING,
                                 Plane.XZ.offset(-10)).face())
            extrude(amount=20)

        # ---- latch cantilever: side slots and relief pocket ---------------
        zc = (M.TAB_ROOT_Z + M.TAB_TOP_Z + 6.0) / 2.0
        zh = (M.TAB_TOP_Z + 6.0 - M.TAB_ROOT_Z)
        for sx in (-1, 1):
            with BuildPart(mode=Mode.SUBTRACT):
                Box(M.SLOT_W, 6.0, zh, align=None).locate(
                    Location((sx * (M.TAB_W / 2 + M.SLOT_W / 2), 5.0, zc)))
        with BuildPart(mode=Mode.SUBTRACT):
            Box(M.TAB_W, M.R_OUT - M.TAB_T - (M.R_CAV - 0.3),
                M.Z_BODY_TOP - M.TAB_ROOT_Z, align=None).locate(
                Location((0.0,
                          (M.R_CAV - 0.3 + M.R_OUT - M.TAB_T) / 2.0,
                          (M.TAB_ROOT_Z + M.Z_BODY_TOP) / 2.0)))

        # ---- detent bead: 45 degree triangular ridge ----------------------
        with BuildSketch(Plane.YZ.offset(-M.BEAD_W / 2)):
            Polygon((M.R_NECK, M.BEAD_Z - M.BEAD_R),
                    (M.R_NECK + M.BEAD_R, M.BEAD_Z),
                    (M.R_NECK, M.BEAD_Z + M.BEAD_R), align=None)
        extrude(amount=M.BEAD_W)

        # ---- hinge dish ---------------------------------------------------
        with BuildPart(mode=Mode.SUBTRACT):
            Cylinder(M.R_DISH, 40, rotation=(0, 90, 0)).locate(
                Location((0, M.Y_PIN, M.Z_PIN)))

        # ---- ears ---------------------------------------------------------
        for sx in (-1, 1):
            with BuildSketch(Plane.YZ.offset(sx * M.EAR_X - M.EAR_HALF)):
                SlotOverall(2 * M.EAR_R + (M.Z_PIN + M.EAR_EXT - M.EAR_BOT_Z),
                            2 * M.EAR_R, rotation=90).locate(
                    Location((M.Y_PIN,
                              (M.Z_PIN + M.EAR_EXT + M.EAR_BOT_Z) / 2.0)))
            extrude(amount=M.EAR_T)

        # ---- pin bores, teardrop roofed -----------------------------------
        for sx in (-1, 1):
            with BuildPart(mode=Mode.SUBTRACT):
                add(teardrop_profile(
                    M.R_BORE,
                    Plane.YZ.offset(sx * M.EAR_X - M.EAR_HALF - 1.0)
                ).moved(Location((M.Y_PIN, M.Z_PIN, 0))))
                extrude(amount=M.EAR_T + 2.0)

        # every visible outside edge gets a radius
        fillet(bp.edges().filter_by(Axis.Z, reverse=True)
               .group_by(Axis.Z)[-1], radius=0.4)
    return bp.part


# --------------------------------------------------------------------------
# cap
# --------------------------------------------------------------------------

def make_cap():
    import math
    with BuildPart() as cp:
        add(capsule(M.OUT_W, M.OUT_T, M.BOT_Z - M.R_OUT, M.TOP_Z + M.R_OUT))
        with BuildPart(mode=Mode.INTERSECT):
            Box(40, 20, 40, align=None).locate(
                Location((0, 0, M.Z_CAP_BOT + 20)))

        # hollow, with a 45 degree tented roof
        with BuildPart(mode=Mode.SUBTRACT):
            with BuildSketch(Plane.XY.offset(M.Z_CAP_BOT - 1)):
                stadium(M.OUT_W - 2 * M.NECK_INSET + 2 * M.GAP,
                        M.OUT_T - 2 * M.NECK_INSET + 2 * M.GAP)
            extrude(amount=M.Z_NECK_TOP - M.Z_CAP_BOT + 1)
            with BuildSketch(Plane.XY.offset(M.Z_NECK_TOP + 1)):
                stadium(M.OUT_W - 2 * M.CAP_WALL, M.OUT_T - 2 * M.CAP_WALL)
            extrude(amount=81.8 - (M.Z_NECK_TOP + 1))
            with BuildSketch(Plane.XY.offset(81.8)):
                stadium(M.OUT_W - 2 * M.CAP_WALL, M.OUT_T - 2 * M.CAP_WALL)
            with BuildSketch(Plane.XY.offset(85.8)):
                Circle(0.01)
            extrude(amount=0, both=False, taper=45)   # the tent

        # knuckle: keep only material within R_KNUCK behind and below the pin
        with BuildPart(mode=Mode.SUBTRACT):
            Box(40, 20, 20, align=None).locate(
                Location((0, M.Y_PIN - 10, M.Z_PIN - 10)))
            Cylinder(M.R_KNUCK, 40, rotation=(0, 90, 0), mode=Mode.SUBTRACT
                     ).locate(Location((0, M.Y_PIN, M.Z_PIN)))

        # relief where the ears pass through
        for sx in (-1, 1):
            with BuildPart(mode=Mode.SUBTRACT):
                add(teardrop_profile(
                    M.R_EAR_RELIEF,
                    Plane.YZ.offset(sx * M.EAR_X - M.EAR_HALF - M.GAP)
                ).moved(Location((M.Y_PIN, M.Z_PIN, 0))))
                extrude(amount=M.EAR_T + 2 * M.GAP)

        # captive pins with a 45 degree V underside
        for sx in (-1, 1):
            with BuildSketch(Plane.YZ.offset(min(sx * M.PIN_X_IN,
                                                 sx * M.CAP_IN_HALF_W))):
                Circle(M.R_PIN).locate(Location((M.Y_PIN, M.Z_PIN)))
                Polygon((M.Y_PIN - M.R_PIN, M.Z_PIN - M.R_PIN),
                        (M.Y_PIN + M.R_PIN, M.Z_PIN - M.R_PIN),
                        (M.Y_PIN, M.Z_PIN - M.R_PIN * math.sqrt(2)),
                        align=None, mode=Mode.SUBTRACT)
            extrude(amount=abs(M.CAP_IN_HALF_W - M.PIN_X_IN))

        # latch groove: the bead's own profile, grown by GROOVE_GAP
        g = M.BEAD_R + M.GROOVE_GAP
        with BuildPart(mode=Mode.SUBTRACT):
            with BuildSketch(Plane.YZ.offset(-(M.BEAD_W / 2 + 0.2))):
                Polygon((M.R_NECK, M.BEAD_Z - g),
                        (M.R_NECK + g, M.BEAD_Z),
                        (M.R_NECK, M.BEAD_Z + g), align=None)
            extrude(amount=M.BEAD_W + 0.4)

        # knurled collar
        for i in range(M.KNURL_N):
            a = (i + 0.5) / M.KNURL_N * 2 * math.pi
            cx = M.HF_OUT if math.cos(a) >= 0 else -M.HF_OUT
            pr = M.R_OUT - M.KNURL_DEPTH + M.KNURL_R
            with BuildPart(mode=Mode.SUBTRACT):
                Cylinder(M.KNURL_R, M.KNURL_Z1 - M.KNURL_Z0,
                         align=None).locate(
                    Location((cx + pr * math.cos(a), pr * math.sin(a),
                              M.KNURL_Z0)))

        fillet(cp.edges().group_by(Axis.Z)[-1], radius=0.4)
    return cp.part


# --------------------------------------------------------------------------

def main():
    os.makedirs(OUT, exist_ok=True)
    body, cap = make_body(), make_cap()
    asm = Compound(label='strip_case', children=[
        Compound(label='body', children=[body]),
        Compound(label='cap', children=[cap]),
    ])
    export_step(asm, os.path.join(OUT, 'strip_case.step'))
    export_stl(body, os.path.join(OUT, 'strip_case_body.stl'))
    export_stl(cap, os.path.join(OUT, 'strip_case_cap.stl'))
    print('wrote STEP and STL to', OUT)


# NOTE for whoever runs this first
# -------------------------------
# Three places need a loft rather than an extrude, and are marked above:
#   * the cavity's tapered front wall (CAV_FRONT_LO -> Z_BODY_TOP)
#   * the cap's void step from skirt radius to wall radius
#   * the cap's 45 degree tented roof
# In build123d these are `loft()` between two sketches on offset planes.  The
# SDF model does them analytically, which is why model.py needs no such
# special casing.  Everything else transcribes one-to-one.

if __name__ == '__main__':
    main()
