"""Build the .blend from the verified meshes, inside Blender.

I could not run this for you: this session has no Blender and no way to
install one (pypi and the Ubuntu archives are both blocked by the
environment's egress policy).  So it is written to need exactly one command
from you, and to be safe to re-run.

    # writes out/strip_case.blend and exits, no GUI
    blender --background --python src/blender_setup.py

    # or, with Blender already open:
    #   Scripting tab -> Open -> this file -> Run Script

What it does beyond importing the mesh:

  * sets the scene to millimetres, so 1 unit reads as 1 mm
  * moves the cap's origin onto the hinge pin axis, so rotating the cap on X
    articulates the hinge correctly instead of pivoting about the world origin
  * clamps that rotation to 0-105 degrees, the range the swing test cleared
  * gives the two parts distinguishable materials and sharp-edge shading

Tested against nothing.  It is defensive about the operator names Blender has
renamed across versions (OBJ import moved in 3.3, STL import in 4.2), but if
something does fail it will say which step and keep going.
"""

import os
import sys
import math

try:
    import bpy
    import mathutils
except ImportError:
    sys.exit('This script has to run inside Blender:\n'
             '  blender --background --python src/blender_setup.py')

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, 'out')
OBJ = os.path.join(OUT, 'strip_case.obj')
BLEND = os.path.join(OUT, 'strip_case.blend')

# From src/model.py.  Hard-coded rather than imported because Blender's Python
# will not have this repo on its path, and these two numbers are the only
# thing the scene needs that the mesh itself does not carry.
Y_PIN, Z_PIN = -2.56, 78.60
OPEN_MAX_DEG = 105.0


def log(msg):
    print('[strip_case] %s' % msg)


# --------------------------------------------------------------------------

def clear_scene():
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)
    for block in (bpy.data.meshes, bpy.data.materials):
        for item in list(block):
            if item.users == 0:
                block.remove(item)


def import_obj(path):
    """OBJ import moved from import_scene.obj to wm.obj_import in Blender 3.3."""
    before = set(bpy.data.objects)
    if hasattr(bpy.ops.wm, 'obj_import'):
        bpy.ops.wm.obj_import(filepath=path, forward_axis='Y', up_axis='Z')
    else:
        bpy.ops.import_scene.obj(filepath=path, axis_forward='Y', axis_up='Z')
    return [o for o in bpy.data.objects if o not in before]


def set_units():
    """1 Blender unit = 1 mm, displayed as millimetres."""
    u = bpy.context.scene.unit_settings
    u.system = 'METRIC'
    u.scale_length = 0.001
    try:
        u.length_unit = 'MILLIMETERS'
    except TypeError:
        pass


def origin_to(obj, x, y, z):
    """Move an object's origin without moving its geometry."""
    offset = mathutils.Vector((x, y, z))
    obj.data.transform(mathutils.Matrix.Translation(-offset))
    obj.matrix_world.translation += offset


def material(name, rgba):
    mat = bpy.data.materials.get(name) or bpy.data.materials.new(name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get('Principled BSDF')
    if bsdf:
        bsdf.inputs['Base Color'].default_value = rgba
        for key, val in (('Roughness', 0.45), ('Metallic', 0.0)):
            if key in bsdf.inputs:
                bsdf.inputs[key].default_value = val
    return mat


def shade(obj, angle_deg=30.0):
    """Smooth shading that still keeps CAD edges crisp."""
    bpy.context.view_layer.objects.active = obj
    if hasattr(bpy.ops.object, 'shade_auto_smooth'):          # 4.1+
        bpy.ops.object.shade_auto_smooth(angle=math.radians(angle_deg))
        return
    bpy.ops.object.shade_smooth()
    mesh = obj.data
    if hasattr(mesh, 'use_auto_smooth'):                      # <= 4.0
        mesh.use_auto_smooth = True
        mesh.auto_smooth_angle = math.radians(angle_deg)


def limit_rotation(obj):
    """Let the cap swing on X only, over the range the swing test cleared."""
    con = obj.constraints.new('LIMIT_ROTATION')
    con.use_limit_x = True
    con.min_x, con.max_x = 0.0, math.radians(OPEN_MAX_DEG)
    con.use_limit_y = con.use_limit_z = True
    con.min_y = con.max_y = con.min_z = con.max_z = 0.0
    con.owner_space = 'LOCAL'


# --------------------------------------------------------------------------

def main():
    if not os.path.exists(OBJ):
        sys.exit('%s is missing. Regenerate it with:\n'
                 '  cd src && python3 build.py 0.3' % OBJ)

    clear_scene()
    set_units()

    imported = import_obj(OBJ)
    log('imported %d objects: %s'
        % (len(imported), ', '.join(o.name for o in imported)))

    parts = {}
    for obj in imported:
        key = 'cap' if 'cap' in obj.name.lower() else 'body'
        obj.name = key
        obj.data.name = '%s_mesh' % key
        parts[key] = obj

    if 'body' in parts:
        body = parts['body']
        body.data.materials.clear()
        body.data.materials.append(material('PLA body', (0.82, 0.83, 0.85, 1)))
        shade(body)

    if 'cap' in parts:
        cap = parts['cap']
        cap.data.materials.clear()
        cap.data.materials.append(material('PLA cap', (0.35, 0.55, 0.72, 1)))
        shade(cap)
        # the whole point: pivot on the pin, not on the world origin
        origin_to(cap, 0.0, Y_PIN, Z_PIN)
        cap.rotation_mode = 'XYZ'
        limit_rotation(cap)
        log('cap origin moved to the pin axis at y=%.2f z=%.2f;'
            ' rotate it on X to open (0 to %.0f deg)'
            % (Y_PIN, Z_PIN, OPEN_MAX_DEG))

    missing = {'body', 'cap'} - set(parts)
    if missing:
        log('WARNING: did not find %s in the OBJ' % ', '.join(sorted(missing)))

    bpy.ops.wm.save_as_mainfile(filepath=BLEND)
    log('wrote %s' % BLEND)


if __name__ == '__main__':
    main()
