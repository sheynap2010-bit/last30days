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
| Envelope | **Ø18 × 94 mm** — a cylinder with hemispherical ends |
| Cavity | 10.8 × 7.0 stadium, 66.5 deep |
| Wall | 3.6 mm at the sides, 5.5 mm front and back, 2.2 mm cap |
| Hinge | **external** — Ø2.7 captive pins on a cap lug, Ø3.5 bores in two body ears |
| Latch | 8 × 1.8 × 11.65 mm cantilever, 45° bead into a matching 45° groove |
| Keyring | Ø7.0 through the bottom dome |
| Mass | 19.9 g in PLA at 100% infill |

Coordinates: Z is the long axis with Z=0 on the bed, −Y is the back (hinge),
+Y is the front (latch).

The cross-section machinery is the same one that drew the earlier flat oval —
a stadium scaled to zero over each end along a circular profile. Setting the
width equal to the thickness collapses the stadium to a circle, so one code
path gives a true cylinder with domed ends and no special casing.

## Why the hinge is outside the skin — and why that seals the case

Opening rotates the cap about the pin, and **every point of the cap behind
the pin sinks as it rotates.** With the pin buried in the back wall the cap's
whole rear skin is behind it, so the body had to be hollowed out of its way —
and that hollow was a 2 mm slot from the cavity to open air, wide enough for
a 0.5 mm strip to escape. Four fixes were built and measured against that,
and all four failed:

| attempted fix | result |
|---|---|
| Full-size knuckle, cut back only at the ears | closes the split line, not the pocket — 2.45 mm |
| Cap tongue reaching into the pocket | seals it, but collides with the rim from 40° to 70° |
| Internal plug entering the bore | cannot print — its first layer is over 66 mm of air |
| Less opening angle | 100°: 2.02, 85°: 1.18, then **floors at 1.09 mm** down to 50° |

Moving the pin 0.15 mm **behind** the skin puts every point of the cap in
front of it. Everything rises on opening, nothing sweeps into the body, and
the split closes to its 0.3 mm print clearance the whole way round. The neck
survives as a complete 360° collar, which is what actually seals the bore.

That leaves only the cap's own lug, which hangs down behind the pin. It is
bounded by a cylinder on the pin axis, and radius is invariant under
rotation, so its swept envelope *is* that cylinder at every angle — the body
just carries a matching cylindrical recess. The recess is 3.0 mm deep into a
5.5 mm back wall: a dish in the outer surface, not a hole into the cavity.

**Consequences.** The model needs no swept-envelope subtraction, no dish and
no knuckle cut; every hinge clearance is now an explicit cylinder on the
axis. The swing test is clean across the *entire* 0–120° sweep rather than
stopping at 100°, and the cap's steep-overhang area dropped from 4.2 % to
0.9 % because all that relief geometry is gone. The cost is visual: the ears
stand 3.1 mm proud of the barrel and the cap's lug 2.7 mm, like a real
flip-top lip balm.

## Printing

Upright, closed, as modelled. Z = 0 is the bed. **Layer height 0.20 mm** —
the clearances are sized as whole multiples of it.

### The clearances are deliberately not all the same

| interface | gap | direction |
|---|---|---|
| split line, body top → cap bottom | **0.60 mm** | horizontal |
| lug → its recess in the body | **0.50 mm** | horizontal underneath |
| pin → bore | 0.40 mm radial | |
| neck → skirt, ear → lug | 0.30 mm | vertical |
| latch slot | 0.45 mm | vertical |

A **vertical** gap is safe at 0.3 mm because the nozzle never enters it. A
**horizontal** one is not: the layer above it is laid down over air, and at
0.3 mm with 0.20 mm layers that is 1.5 layers — the cap's first layer lands
almost on the body's top face and welds to it. 0.60 mm is three layers, so at
least two are completely empty. This is the single change most likely to
decide whether the cap comes free.

The latch slot runs the other way: it has to stay **under** the 0.5 mm strip
thickness, because it leads from the cap's interior, down past the tab and
out through the skin. At 0.6 mm a strip walks straight out of it. 0.45 mm
blocks the strip, and being vertical it still prints cleanly.

### Slicer settings that matter

- **0.20 mm layers.** Not 0.25 or 0.3 — the split stops being a whole number
  of layers and the gap starts closing.
- **Gap fill off** (PrusaSlicer/Orca: "Fill gaps" / minimum feature size).
  Otherwise the slicer may lay a thin extrusion into the 0.45 mm latch slot
  and weld the tab solid.
- 3 perimeters. The thinnest structural section is the 1.0 mm membrane
  between the cavity and the latch pocket.
- **No supports anywhere except under the bottom dome.**

### Two bottoms are supplied

| file | bottom | printing |
|---|---|---|
| `strip_case.stl` | full hemisphere | touches the bed at a point — **needs a brim and support under the lowest ~3 mm** |
| `strip_case_flatbottom.stl` | truncated 2.64 mm | **124 mm² of bed contact, 45° at the bed, no support at all** |

The flat-bottom one is the safer print by a wide margin, and 2.64 mm off a
9 mm dome is hard to see on a Ø18 barrel. Everything above the bottom dome is
identical between the two.

### First open

The two parts are never fused in the model — intersection is 0.0000 mm³ — but
first-layer squish or a little stringing can tack them. Open it slowly the
first time; it should break free under light thumb pressure. If it resists,
run a thin blade round the split line rather than forcing it.

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

## OPEN ISSUE: the closed case is not sealed at the back

`check.py` section 4 measures this directly, and it currently reports that an
8 × 0.5 mm strip **can** escape. Do not treat this design as strip-tight yet.

Measured gap between the body's top and the cap's bottom, around the
perimeter of the closed case:

| bearing | gap |
|---|---|
| front and sides (0–75°, 285–345°) | **0.31 mm** — just the print clearance |
| 90° / 270° (sides) | 0.80 mm |
| 105° / 255° | 1.69 mm |
| 150° / 210° (back) | **2.02 mm** |

### Why it cannot be closed with the pin where it is

Opening rotates the cap about the pin. Every point of the cap *behind* the
pin sinks as it rotates, so the body has to be hollowed out of its way — and
that hollow is the opening. The pin sits inside the back wall, so the cap's
whole rear skin is behind it.

Four fixes were built and measured; all failed:

- **Full-size knuckle** (radius 5.02 instead of 4.54, cut back only in the
  two ear bands): closes the split line but not the pocket — 2.45 mm.
- **Cap tongue** reaching below the split into the pocket: seals it (the
  1.2 mm ball test passes), but the tongue swings forward and collides with
  the body's rim from 40° to 70°.
- **Internal plug** entering the bore: cannot print. Anything that enters the
  bore starts its first layer over 66 mm of open air.
- **Less opening angle**: 100° → 2.02 mm, 85° → 1.18 mm, and then it floors
  at **1.09 mm** all the way down to 50°. The knuckle sets that floor, not
  the tail. Still passes a 0.5 mm strip.

### The fix that does work

Move the pin *behind* the rear skin, at y ≈ −6.35. Then every point of the
cap is in front of the pin, everything rises on opening, nothing sweeps into
the body at all — the neck stays a full 360° collar, the rim stays intact,
and the split closes to 0.3 mm the whole way round.

The cost is visual: the hinge becomes an external boss. A Ø3.5 bore needs
1.2 mm of material round it, so the ears stand about 3 mm proud of the back,
with a matching strap on the cap. That is what a real flip-top lip balm looks
like, but it is not the flush silhouette this design was asked for. **It has
not been built** — it needs a decision first, because it trades the
appearance requirement against the containment requirement.

## Deviations from the brief, and why (continued)

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
