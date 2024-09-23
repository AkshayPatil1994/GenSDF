import numpy as np
import trimesh

def create_torus(r_major, r_minor, center, scale_factor, bounding_box_min, bounding_box_max):
    """
    Create a torus and fit it within a given bounding box.
    
    Parameters:
    - r_major: Major radius of the torus (radius from center to center of tube)
    - r_minor: Minor radius of the torus (radius of the tube itself)
    - center: Center of the torus (before transformation)
    - scale_factor: Factor by which to scale the torus to fit into the bounding box
    - bounding_box_min: Minimum point of the bounding box
    - bounding_box_max: Maximum point of the bounding box
    
    Returns:
    - A scaled and transformed trimesh torus object
    """
    # Create the torus
    torus = trimesh.creation.torus(major_radius=r_major, minor_radius=r_minor)
    
    # Scale the torus to fit inside the bounding box
    torus.apply_scale(scale_factor)
    
    # Translate the torus to the center of the bounding box
    bounding_box_center = (bounding_box_min + bounding_box_max) / 2.0
    translation = bounding_box_center - center
    torus.apply_translation(translation)
    
    return torus

# Define the bounding box (from the input)
bounding_box_min = np.array([15.50003566, 3.50008895, 3.5])
bounding_box_max = np.array([16.5, 4.49997042, 4.5])

# Torus parameters
r_major = 0.3  # Major radius of the torus
r_minor = 0.1  # Minor radius of the torus
center = np.array([0.0, 0.0, 0.0])  # Center of the torus before transformation

# Calculate the scale factor to fit the torus inside the bounding box
# We use the distance between the min and max bounding box points to set a rough scale.
bounding_box_size = bounding_box_max - bounding_box_min
scale_factor = min(bounding_box_size) / (2 * (r_major + r_minor))

# Create the torus and fit it within the bounding box
torus = create_torus(r_major, r_minor, center, scale_factor, bounding_box_min, bounding_box_max)

# 0.01 is the target size of each triangle after subdivision
torus = torus.subdivide_to_size(max_edge=0.005)

# Export the torus mesh to an OBJ file
torus.export('torus.obj')

print("Torus mesh has been saved to 'torus.obj'.")
