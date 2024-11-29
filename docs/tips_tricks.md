## Tips & Tricks

With any software there are some tips and tricks that allow the user to efficiently and effectively get the results as desired. In this document, some details are provided to avoid the common pitfalls in the workflow.

1. Make sure that the geometry has outward vertex normals contained in its definition. If this information is not available, then simply pass the OBJ file through the python script where you load the OBJ and export it again using the trimesh native export function as 
   `mymesh.export("mygeometry.obj",include_normals=True)`
2. In some cases you will observe that the computational grid is relatively finer than the surface triangulation on the geometry. Sometimes this can lead to holes in the generated SDF. To resolve this issue, you can use the `wrapwrap` [https://github.com/ipadjen/wrapwrap] tool to shrink wrap the geometry with a finer surface triangulation. This is especially true for geometry with large triangular faces such as cities, square blocks, and other flat or box like geometry.