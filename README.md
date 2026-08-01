# Strip case — print-in-place flip-top

A flat-oval capsule that holds 10 test strips (60 × 8 × 0.5 mm) and reads as a
lip balm stick. Prints upright, closed, in one piece, no supports.

```
src/geom.py             SDF kernel, marching tetrahedra, STL writer
src/model.py            the model — all geometry and every parameter
src/build.py            polygonise and export
src/check.py            verification; every number in the report comes from here
src/diag.py             locates what check.py counts
src/build123d_model.py  B-rep transcription for STEP export (NEVER RUN — see below)
src/blender_setup.py    builds the .blend inside Blender (NEVER RUN — see below)
out/strip_case.stl      both parts, closed, print-in-place — this is the print
out/strip_case_body.stl,
out/strip_case_cap.stl  the same two solids split out, world coordinates
out/strip_case.obj      both solids as named objects, for Blender
out/verification.txt    the full check output
out/render.png          shaded 3/4 view, closed and open at 105 deg
out/views.png           flat silhouettes: front, side, side open
out/strip_case_stl.zip  all three STLs, compressed
```

Regenerate: `cd src && python3 check.py 0.3` (no dependencies, stdlib only).
The argument is the marching grid in mm; it builds the meshes, writes `out/`,
and runs every check.

## What it is

| | |
|---|---|
| Envelope | 21.1 W × 12.4 T × 96 H mm — stadium section, 1.70 : 1 |
| Ends | half-ellipsoids, R10.55 in front view — a true U, not a swept flat |
| Cavity | 10.8 × 7.0 stadium, 66.5 deep |
| Wall | 2.7 mm around the cavity, 2.2 mm cap, 2.0 mm skirt |
| Hinge | Ø2.7 captive pins on the cap, Ø3.5 bores in two body ears at x = ±6 |
| Latch | 8 × 1.8 × 11.65 mm cantilever, 45° bead into a matching 45° groove |
| Keyring | Ø7.0 through the bottom dome, axis along the thickness |
| Mass | 17.0 g in PLA at 100% infill |

Coordinates: Z is the long axis with Z=0 on the bed, X is the wide direction,
−Y is the back (hinge), +Y is the front (latch).

## The end shape — the thing that decides whether it reads as a vape

Both ends scale the *whole* stadium cross section to zero along a circular
profile whose radius is the body's half width. The front view of each end is
therefore an exact semicircle spanning the full width: a **U**.

The obvious alternative — sweeping the stadium, i.e. offsetting it in 3D —
gives a bottom that is flat across the middle with rounded corners. That is
what the first version did, and it read unmistakably as a tube. It is a
different surface, not a smaller radius, and no amount of filleting fixes it.

**The cost: a true U cannot print support-free standing up.** The bottom is
tangent to the bed at a single point, and the last few millimetres overhang
past 45°. Measured on the R10.55 bottom:

| cut off the tip | bed contact | wall angle at the bed |
|---|---|---|
| 0.00 mm (a true U) | a point | 90° |
| 1.00 mm | 37 mm² | 65° |
| 2.00 mm | 71 mm² | 54° |
| **3.09 mm** | **103 mm²** | **45° — the support-free limit** |

3.09 mm off a 10.55 mm dome leaves a round-over, not a U — which is the shape
that was rejected. So the model ships at `BOT_FLAT = 0.0`: fully round, print
with a brim and support under the lowest ~3 mm. One line in `model.py`
(`BOT_FLAT = 3.09`) switches to the support-free bottom if you change your
mind.

This is the one place this design knowingly breaks the original brief's "no
support material".

## How it is modelled, and why not with booleans

Everything is a signed distance function composed analytically. A fillet is a
smooth-min, evaluated on the surface itself — a real constant-radius blend, not
a mesh operation. Two things follow:

- **Collision tests are exact.** Asking whether a point is inside the body is
  one function evaluation. No ray casting, no tessellation error in the
  answer. The swing test and the intersection test are done on the functions,
  not on the meshes.
- **The rear hinge relief is swept, not guessed.** `cap_swept()` returns the
  cap's closest approach to a point over the whole 0–105° opening range, and
  the body subtracts it. Three hand-derived relief radii failed before this;
  the swept envelope is correct by construction.

Marching tetrahedra (not cubes) polygonises the result: 16 cases instead of
256, watertight by construction because vertices are keyed by lattice edge.

It then gets one pass of `smooth_project()`. Marching tetrahedra puts each
vertex wherever a lattice edge happens to cross the surface, which leaves
slivers and an uneven distribution — and *that*, not the grid pitch, is what
makes a curved surface look faceted under shading. The pass relaxes each
vertex toward its neighbours and then Newton-steps it back onto the exact
isosurface along the SDF gradient, so the triangles improve while every
vertex still lies on the true surface. Vertices on a crease are held still,
so the split line, the latch slot and the flats stay sharp.

Creases are detected from *area-weighted* face normals. Weighting matters: an
unweighted average lets sliver triangles, whose normals are numerical noise,
masquerade as creases — the first attempt flagged 85 % of the surface and the
smoothing did nothing.

## Deviations from the brief, and why

**Cavity 10.8 × 7.0, not 10 × 6.5.** A 10 × 6.5 stadium does not actually
contain an 8 × 5 bundle: the bundle's corners sit 3.363 mm from the section's
arc centre against a 3.25 mm radius, so they foul by 0.11 mm. 10.8 × 7.0 is the
smallest stadium of that proportion that swallows the bundle, with 0.21 mm to
spare.

**Wall 2.7 mm around the cavity, not 2.0–2.5.** Forced by the other two
constraints: a 1.70 : 1 outer oval around a 10.8 × 7.0 cavity, with the pin
bore buried in the back wall, cannot have a thinner wall. The brief says the
strips and the wall win over the envelope, so the envelope moved instead —
21.1 × 12.4 rather than 24 × 14.

**The cavity's front wall tapers** from 3.5 to 2.75 mm off centre over the top
7.85 mm. That is what buys the neck enough material to carry a 2.0 mm skirt,
which is what gives the latch groove 1.45 mm of cap wall behind it. Slope is
15° from vertical. The opening still measures 6.25 mm front to back.

**Ears at x = ±6.0, pin at y = −2.56.** Not centred in the back wall. At the
ears' outer edge the oval's skin is only ±5.51 mm from the centre plane, so a
pin further back would break the bore out through the outside face — the exact
failure mode the brief called out. This position gives exactly 1.20 mm of
material all round the bore.

**The back of the neck is open** from z ≈ 74.6 to 78.0, where the swept relief
clears the knuckle. The bore is fully enclosed below that, and the strips reach
only z = 71.5, so they are 3.1 mm clear of the opening. The cap covers it when
closed.

**The split line is flush at the front and sides but not across the back.**
Making it flush there wants a knuckle of radius 4.39 about the pin, and the
ears root through exactly that annulus — a knuckle that large sweeps into them
at about 95°. The knuckle is cut back to 3.84 instead, which leaves a 1.4 mm
gap across the back of the split. It reads as an ordinary hinge gap, and it is
the one place the brief's "no step, one continuous surface" is not met.

## Deliverables that could not be produced here

Outbound network access is blocked by this environment's egress policy —
`pypi.org`, `archive.ubuntu.com` and `registry.npmjs.org` all return 403 through
the proxy. So:

- **No STEP.** build123d, CadQuery and OCP could not be installed.
  `src/build123d_model.py` is a faithful transcription that exports STEP, but
  **it has never been executed** and should be treated as unrun source. It
  imports its parameters from `model.py` so the two cannot drift.
- **No .blend.** Blender and `bpy` are not installed and could not be fetched,
  and a `.blend` is a binary dump of Blender's internal structures — not
  something worth hand-writing blind. Instead, one command on a machine that
  has Blender produces it:

  ```
  blender --background --python src/blender_setup.py     # writes out/strip_case.blend
  ```

  That script does more than import the mesh. It sets the scene to
  millimetres, and it moves the cap's origin onto the hinge pin axis
  (y = −2.56, z = 78.60) with a rotation limit of 0–105° on X — so the cap
  articulates about the real hinge instead of the world origin, over exactly
  the range the swing test cleared. It has also never been run; it is
  defensive about the operators Blender renamed in 3.3 and 4.2.

  If you would rather do it by hand, `out/strip_case.obj` carries the two
  solids as named objects (`body`, `cap`) and imports in one step — but the
  cap will pivot about the world origin until you move its origin yourself.

The STL is the verified deliverable. It comes from `model.py`, and `check.py`
runs against the same functions the STL was generated from.

## Verified numbers

From `out/verification.txt`, 0.30 mm marching grid:

| | |
|---|---|
| Topology | both solids watertight, 0 non-manifold edges, 1 component each |
| Closed-position intersection | **0.0000 mm³** (0 of 340k lattice cells inside both) |
| Swing 0–105° | **no hinge contact at any station** |
| Swing 110°+ | contact — this is the end stop, past the 100° requirement |
| Latch at 5° | 0.12 mm contact on the bead only: the detent releasing |
| Bottom | full R10.55 round, 0 mm² bed contact — needs a brim + tip support |
| Wall round the pin hole | 1.20 mm (need ≥ 1.2) |
| Bundle clearance | +0.21 mm at the corners |
| Cavity enclosed to | z = 74.6; strips reach z = 71.5 → 3.1 mm margin |
| Latch force | 10.1 N (target 8–12) |
| Latch peak strain | 0.91 % (PLA yields near 2–3 %) |
| Mass | 16.95 g at 100 % infill |
| Steep overhangs | body 2.16 %, cap 4.17 % — most of the body figure is the round bottom, by design |

The swing test rotates the cap about the pin axis in 5° steps and evaluates
the body's SDF at ~6000 cap-surface points per station. Contact at 5° is
confined to the latch bead and is reported separately: the cantilever is
meant to deflect there. Everything else is hinge interference and there is
none of it below 110°.

### The overhangs that remain

Counted from SDF gradients, not facet normals, at 46° (1° of slack on the
45° limit). Surfaces with material within 0.6 mm directly below are counted
separately — they bridge a gap rather than print over air, which is the
premise of a print-in-place mechanism.

What is left is genuinely unsupported, and it is worth knowing where:

- **The cap's bottom rim over the hinge relief at the back** (worst case,
  90°). The body's rear shoulder is scalloped away so the cap can rotate, so
  the cap's rim there overhangs 1–2 mm of air. Unavoidable on a hinge that
  opens past 100°; it bridges.
- **The ears' undersides** where the swept relief rounds them off, about
  1–2 mm across each.
- **The teardrop ridges** over the keyring and pin bores. The two 45° faces
  meet in a ridge carrying a 0.3 mm radius, so the very top of that radius is
  briefly horizontal. Not a real overhang.
- **The pin's first layer** floats 0.4 mm above the bore's floor — that is
  the mandated Ø2.7-in-Ø3.5 radial clearance, and it is how every
  print-in-place pin hinge behaves.

## Printing

Upright, closed, as modelled. Z = 0 is the bed.

- **Brim, plus support under the bottom ~3 mm.** The round bottom demands it.
- 0.3 mm clearance on every print-in-place gap; 0.4 mm radial on the pin
- 0.2 mm layers, 3 perimeters. The thinnest structural section is the 1.19 mm
  tip of the latch cantilever.
- No supports. Every hole with a horizontal axis — both pin bores and the
  keyring — is teardropped with a 45° gable roof rather than left round.
- The pin's underside is a 45° V rather than a cylinder, so it self-supports.
