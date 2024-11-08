import trimesh
# User Input Data
objname = 'armadillo.obj'
# Load the mesh
mesh = trimesh.load(objname)
# First fix the normals to point outwards
mesh.fix_normals()
# Export the geometry with normals information
mesh.export('armadillo_withnormals.obj',include_normals=True)
