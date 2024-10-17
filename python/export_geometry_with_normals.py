import trimesh

# USER INPUT DATA
objname = 'torus_withoutnormals.obj'
# Load the mesh
mesh = trimesh.load(objname)
# First fix the normals to point outwards
mesh.fix_normals()
# Export the geometry with normals information
mesh.export('torus_withnormals.obj',include_normals=True)
