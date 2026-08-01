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
out/                    STL, OBJ, verification log
```

Regenerate: `cd src && python3 check.py 0.3` (no dependencies, stdlib only).
The argument is the marching grid in mm; it builds the meshes, writes `out/`,
and runs every check.

## What it is

| | |
|---|---|
| Envelope | 21.1 W × 12.4 T × 88 H mm — stadium section, 1.70 : 1, 7.1 : 1 long |
| Cavity | 10.8 × 7.0 stadium, 68.65 deep |
| Wall | 2.7 mm around the cavity, 2.2 mm cap, 2.0 mm skirt |
| Hinge | Ø2.7 captive pins on the cap, Ø3.5 bores in two body ears at x = ±6 |
| Latch | 8 × 1.8 × 11.65 mm cantilever, 45° bead into a matching 45° groove |
| Keyring | Ø4.5 through the bottom dome, axis along the thickness |
| Mass | ~17 g in PLA at 100% infill |

Coordinates: Z is the long axis with Z=0 on the bed, X is the wide direction,
−Y is the back (hinge), +Y is the front (latch).

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

**The back of the neck is open** from z ≈ 74.1 to 77.65, where the swept relief
clears the knuckle. The bore is fully enclosed below that, and the strips reach
only z = 69, so they are 5.1 mm clear of the opening. The cap covers it when
closed.

## Deliverables that could not be produced here

Outbound network access is blocked by this environment's egress policy —
`pypi.org`, `archive.ubuntu.com` and `registry.npmjs.org` all return 403 through
the proxy. So:

- **No STEP.** build123d, CadQuery and OCP could not be installed.
  `src/build123d_model.py` is a faithful transcription that exports STEP, but
  **it has never been executed** and should be treated as unrun source. It
  imports its parameters from `model.py` so the two cannot drift.
- **No .blend.** Blender and `bpy` are not installed and could not be fetched.
  `out/strip_case.obj` carries the two solids as named objects (`body`, `cap`)
  and imports into Blender in one step.

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
| Bed contact | 150.3 mm² (need ≥ 100) |
| Wall round the pin hole | 1.20 mm (need ≥ 1.2) |
| Bundle clearance | +0.21 mm at the corners |
| Cavity enclosed to | z = 74.3; strips reach z = 69.0 → 5.3 mm margin |
| Latch force | 10.1 N (target 8–12) |
| Latch peak strain | 0.91 % (PLA yields near 2–3 %) |
| Mass | 17.16 g at 100 % infill |
| Steep overhangs | body 0.76 %, cap 3.16 % of sampled surface — see below |

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

- 0.3 mm clearance on every print-in-place gap; 0.4 mm radial on the pin
- 0.2 mm layers, 3 perimeters. The thinnest structural section is the 1.19 mm
  tip of the latch cantilever.
- No supports. Every hole with a horizontal axis — both pin bores and the
  keyring — is teardropped with a 45° gable roof rather than left round.
- The pin's underside is a 45° V rather than a cylinder, so it self-supports.
