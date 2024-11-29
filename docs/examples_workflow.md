# Workflow using GenSDF

Once you have sucessfully compiled the software, you can follow the instructions below to run the example case listed in the `example` directory of the repository.

1. Fetch the "Armadillo" OBJ file from the openly available location. Credits to `alecjacobson` making the OBJ file openly available. 

   ```
   cd example/data      
   chmod +x fetchArmadillo.sh
   ./fetchArmadillo.sh
   ```
   This will download the required geometry to the data folder.
2. Since the geometry is not scaled to fit on grid provided in the example, we must scale it to the appropriate size.
   ```
   python rescale_geometry.py
   ```
    In this step the python script will rescale the geometry to fit within the grid as well as output the vertex outward normals in the OBJ file. This step is crucial as the software will raise an error if the geometry does not contain vertex normals information.
3. Navigate to the serial code and compile it and return to the examples folder
   ```
   cd ../../src/serial
   make
   cp gensdf_serial ../../example/
   cd ../../example/
   ```
4. Execute the code
   ```
   ./gensdf_serial
   ```
   You will get output that looks as shown below.
   ```
     ░▒▓██████▓▒░░▒▓████████▓▒░▒▓███████▓▒░ ░▒▓███████▓▒░▒▓███████▓▒░░▒▓████████▓▒░ 
     ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        
     ░▒▓█▓▒░      ░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        
     ░▒▓█▓▒▒▓███▓▒░▒▓██████▓▒░ ░▒▓█▓▒░░▒▓█▓▒░░▒▓██████▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓██████▓▒░   
     ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        
     ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        
     ░▒▓██████▓▒░░▒▓████████▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓███████▓▒░░▒▓███████▓▒░░▒▓█▓▒░        
    *** Input file sucessfully read ***
    *** Successfully read the CaNS grid ***
    *** Sucessfully finished setting up the grid spacing ***
    Successfully read OBJ file: data/armadillo_withnormals.obj
    Number of vertices:        49990
    Number of normals:        49990
    Number of faces:        99976
    Geometry is bounded by (minimum)   12.000000000000000       0.50000000000000000       0.50000000000000000     
    Geometry is bounded by (maximum)   16.000000000000000        5.2650256700000000        4.1345718800000002     
    *** Min-Max Index-Value pair ***
    Min-Max x:         182   11.375000000000000      |         266   16.625000000000000     
    Min-Max y:           1   3.1250000000000000E-002 |          94   5.8437500000000000     
    Min-Max z:           1   3.1250000000000000E-002 |          76   4.7187500000000000     
    -- Finished pre-processing geometry in   0.35796300000000003      seconds...
    -- Estimated Minimum Memory usage:  0.14 GiB(s)...
    *** Calculating the signed-distance-field | u-faces ***
    ||||||||||||||||||||||100.00%    99976/   99976 Elapsed:    11.92s Remaining:     0.00s
    *** Writing output data to file ***
    -- Done with file write in    2.6056000000000523E-002 seconds...
   ```
4. Once the SDF is generated for the four components, you can visualise the results in paraview using the `read_xmf_constant_grid.xmf` script. Please use the default XMF reader in paraview as the other two options will crash the session.
5. To visualise the interface, plot the contour with a value of 0.

## MPI-version

To use the MPI version most of the steps are similar except that you will need to compile the mpi version of the code hosted in the `src/mpi` directory.

- Important to make sure that the MPI library is correctly sourced in the terminal environment
- MPI-version of the code does not run with `mpirun -np 1 ./gensdf_mpi` as the decomposition relies on atleast 2 divisions in the streamwise direction.
- A note of warning that since the algorithm only computes the SDF for a bounding box based on the geometry, using increasing MPI ranks for SDF computations for a small bounding box may lead to problems. 