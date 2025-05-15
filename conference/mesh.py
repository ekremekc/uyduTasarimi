import gmsh
import os
import sys
from helezon.solver_utils import start_time, execution_time
import params
t0 = start_time()

dir_path = os.path.dirname(os.path.realpath(__file__))
print(dir_path)

filename = "MIKRO_UYDU_VC_REORIENTED"

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add(filename)
gmsh.option.setString("Geometry.OCCTargetUnit", "M")

path = os.path.dirname(os.path.abspath(__file__))

gmsh.model.occ.importShapes(os.path.join(path, "GeomDir/" + filename + ".stp"))

thermal_pad = params.thermal_pad
thickness = 0.005

if thermal_pad:
    # gmsh.model.occ.removeAllDuplicates()

    # add pad geometry
    px, py, pz = -0.08, -0.12, 0.005
    pad = gmsh.model.occ.addBox(px, py, pz, -2*px , 0.01+0.12, thickness)
    bounding_box_tag = pad
    # shift the rest
    gmsh.model.occ.translate([(3,2),(3,3)],0, 0, thickness)
    gmsh.model.occ.removeAllDuplicates()


if not thermal_pad:
    gmsh.model.occ.removeAllDuplicates()


gmsh.model.occ.synchronize()

print(gmsh.model.getEntities(dim=3))


lc = 5E-2

# led_tag = 2

# # Mesh refinement
# gmsh.model.mesh.field.add("Constant", 1)
# gmsh.model.mesh.field.setNumbers(1, "VolumesList", [led_tag])
# gmsh.model.mesh.field.setNumber(1, "VIn", lc / 10)
# gmsh.model.mesh.field.setNumber(1, "VOut", lc)

# gmsh.model.mesh.field.setAsBackgroundMesh(1)

# Get bounding box in Z
entities = gmsh.model.getEntities(dim=3)
volumes = [e[1] for e in entities if e[0] == 3]
# if not thermal_pad:
#     gmsh.model.occ.removeAllDuplicates()
#     bounding_box_tag = volumes[0]
bbox = gmsh.model.getBoundingBox(3, volumes[0])
zmin = bbox[2]
zmax = bbox[5]
dz = (zmax - zmin) / 2
# xmin = bbox[0]
# xmax = bbox[3]
# dx = (xmax - xmin) / 2

# Define a mesh size field that creates 10 layers in Z
field_id = gmsh.model.mesh.field.add("MathEval")
# Mesh size is small in Z, forces Z-layering
gmsh.model.mesh.field.setString(field_id, "F", f"{dz}")

# Set it as the background mesh field
gmsh.model.mesh.field.setAsBackgroundMesh(field_id)

gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
gmsh.option.setNumber("Mesh.Algorithm", 6)
gmsh.option.setNumber("Mesh.Algorithm3D", 10)
gmsh.option.setNumber("Mesh.Optimize", 1)
gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
gmsh.model.mesh.generate(3)

sur_tags = gmsh.model.getEntities(dim=2)

vol_tags = gmsh.model.getEntities(dim=3)
# print(sur_tags)
print(vol_tags)



if thermal_pad:
    gmsh.model.addPhysicalGroup(3, [1,4], tag=1) # chasis and solar panel
    gmsh.model.addPhysicalGroup(3, [2,3], tag=2) # equipment
    gmsh.model.addPhysicalGroup(3, [5], tag=3) # thermal pad

if not thermal_pad:
    gmsh.model.addPhysicalGroup(3, [1,4], tag=1) # chasis and solar panel
    gmsh.model.addPhysicalGroup(3, [2,3], tag=2) # equipment

for surface in sur_tags:
    gmsh.model.addPhysicalGroup(2, [surface[1]], tag=surface[1])

# for volume in vol_tags:
#     gmsh.model.addPhysicalGroup(3, [volume[1]], tag=volume[1])


gmsh.model.occ.synchronize()

if "-nopopup" not in sys.argv:
    gmsh.fltk.run()

gmsh.write("{}.msh".format(dir_path + "/MeshDir/" + filename))
gmsh.write("{}.stl".format(dir_path + "/MeshDir/" + filename))

gmsh.finalize()

from helezon.io_utils import write_xdmf_mesh

write_xdmf_mesh(dir_path + "/MeshDir/" + filename, dimension=3)

execution_time(t0, string="\nTotal execution time for mesh generation: ")
