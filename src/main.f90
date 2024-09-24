program sdfgenerator
    !
    ! Signed-Distance-Field Generator 
    !
    use utils, only : read_inputfile, read_cans_grid, read_obj, getbbox, &
                      tagminmax, generate_directions, compute_scalar_distance, &
                      tag_narrowband_points, get_signed_distance
    use utils, only : inputfilename, nx, ny, nz, numberofrays, lx, ly, lz, dx, dy, dz, &
                      sdfresolution, scalarvalue, pbarwidth, dp, sp, r0, ng, &
                      non_uniform_grid, dx_inverse, dy_inverse, dz_inverse, &
                      xp, yp, zp, xf, yf, zf, gpu_threads
#if defined(_ISCUDA)                      
    use utils, only : get_signed_distance_cuda
    use cudafor
#endif
    implicit none

    ! Input geometry related data
    real(dp), allocatable, dimension(:,:) :: vertices
    integer, allocatable, dimension(:,:) :: faces
    integer :: num_faces, num_vertices
    ! Bounding related data
    real(dp), dimension(3) :: bbox_min, bbox_max  
    ! Min-Max bounds for x, y, and z
    integer :: sx, ex, sy, ey, sz, ez
    ! Arrays that store narrow band indices
    logical, allocatable, dimension(:,:,:) :: narrowmask
    integer, allocatable, dimension(:,:) :: narrowbandindices
    integer :: numberofnarrowpoints
    ! Data structure for masking
    real(dp), dimension(3, 19) :: directions
    real(dp), allocatable, dimension(:,:,:) :: pointinside, finaloutput
    ! Time loggers
    real(dp) :: time1, time2, elapsedTime
#if defined(_ISCUDA)
    ! CUDA related data
    ! -- thread and block information
    type(dim3) :: threads_per_block, num_blocks
    ! -- input arguments for get_signed_distance
    real(dp), allocatable, dimension(:), device :: d_xf, d_yp, d_zp !d_yf, d_zf, d_xp, d_yp, d_zp
    integer, device :: d_numrays, d_num_faces, d_numberofnarrowpoints
    integer, allocatable, dimension(:,:), device :: d_faces
    real(dp), allocatable, dimension(:,:), device :: d_vertices
    real(dp), dimension(3,19), device :: d_directions
    integer, allocatable, dimension(:,:), device :: d_narrowbandindices
    real(dp), allocatable, dimension(:,:,:), device :: d_pointinside
    type(cudaDeviceProp) :: prop
    ! Device details
    integer :: device, ierr, max_threads_per_block, max_blocks_x, max_blocks_y, max_blocks_z
#endif
    !
    ! PROGRAM BEGINS
    !
    ! Log CPU time at the start
    call cpu_time(time1)
    ! Read the input file (all ranks read inputfile)
    call read_inputfile()
    ! Read the grid (all ranks read the grid)
    call read_cans_grid('data/', dp, ng, r0, non_uniform_grid, xp, yp, zp, xf, yf, zf,dz)
    ! Allocate memory for the grid
    allocate(pointinside(nx,ny,nz))
    ! Read the OBJ file
    call read_obj(inputfilename, vertices, faces, num_vertices, num_faces)
    ! Query the bounding box
    call getbbox(vertices, num_vertices, bbox_min, bbox_max)    
    ! Print checks
    print *, "Geometry is bounded by (minimum)", bbox_min
    print *, "Geometry is bounded by (maximum)", bbox_max
    print *, "Geometry has ",num_faces, "number of faces.."
    print *, "Geometry has ",num_vertices, "number of vertices.."
    !Use bounding box extents plus a small value (order 5*grid size) as control points
    !Note: For z the buffer is 5*smallest grid point
    call tagminmax(xf,yp,zp,bbox_min,bbox_max,Nx,Ny,Nz,dx,dy,dz(2),sx,ex,sy,ey,sz,ez)
    ! Now allocate the narrow band indices search array
    allocate(narrowmask(nx,ny,nz))
    ! Generate the coorindate directions
    call generate_directions(directions)
    ! Log CPU time at the end of pre-processing
    call cpu_time(time2)
    elapsedTime = (time2-time1)
    print *, "- - - - - - - - Finished Pre-Processing in ", int(elapsedTime), "seconds..."
    ! Set a required value outside the bounding box (User specified)
    pointinside = scalarvalue
    ! Loop over all data 
    print *, "*** 1st - Pass: Distance calculation ***"
    call cpu_time(time1)        ! Log CPU time at the start of the operation    
    ! First compute the scalar distance
    call compute_scalar_distance(sx,ex,sy,ey,sz,ez,xf,yp,zp,num_vertices,vertices,sdfresolution,pointinside)
    print *, "*** 2nd - Pass: Narrow band tagging - Ld <",4.0_dp*min(dx,dy,dz(2)),"  ***"
    ! Check indices for where values are < sdfresolution*min(dx,dy,dz)    
    narrowmask = pointinside < 4.0_dp*min(dx,dy,dz(2))
    numberofnarrowpoints = count(narrowmask)
    print *, "Narrow band point: ", numberofnarrowpoints, "of", (ex-sx)*(ey-sy)*(ez-sz), "total points."
    ! Allocate the values and indices to the narrowbandindices array
    allocate(narrowbandindices(3,numberofnarrowpoints))
    call tag_narrowband_points(sx,ex,sy,ey,sz,ez,narrowmask,narrowbandindices)
#if defined(_ISCUDA)
    ! Now compute the ray intersection only for the narrowband points
    print *, "*** 3rd - Pass: Computing the sign for the distance | Computed on the GPU ***" 
    ! Allocate memory for auxiliary data
    allocate(d_xf(nx),d_yp(ny),d_zp(nz))
    allocate(d_faces(3,num_faces), d_vertices(3,num_vertices))
    allocate(d_narrowbandindices(3,numberofnarrowpoints),d_pointinside(nx,ny,nz))
    ! Query the GPU being used
    ierr = cudaGetDeviceProperties(prop, 0)
    if (ierr /= 0) then
        print *, "Error: Unable to get device properties."
        stop
    end if
    max_threads_per_block = prop%maxThreadsPerBlock
    max_blocks_x = prop%maxGridSize(1)
    max_blocks_y = prop%maxGridSize(2)
    max_blocks_z = prop%maxGridSize(3)
    ! Print information to screen
    print *, "Device Name: ", trim(prop%name)
    print *, "Max Threads Per Block: ", max_threads_per_block
    print *, "  X-dimension: ", max_blocks_x
    print *, "  Y-dimension: ", max_blocks_y
    print *, "  Z-dimension: ", max_blocks_z
    if (gpu_threads > max_threads_per_block) then
        print *, "** WARNING: User attempting to use more threads than allowed **"
    endif
    ! Define the block_size and grid_size
    threads_per_block = dim3(gpu_threads, 1, 1)
    num_blocks = dim3(ceiling(real(numberofnarrowpoints) / gpu_threads, kind=dp), 1, 1)
    print *, "GPU - Threads: ", threads_per_block
    print *, "GPU - Blocks: ", num_blocks
    ! First copy all data from host (CPU) to device (GPU)
    d_xf = xf
    d_yp = yp
    d_zp = zp
    d_numrays = numberofrays
    d_num_faces = num_faces
    d_numberofnarrowpoints = numberofnarrowpoints
    d_faces = faces
    d_vertices = vertices
    d_directions = directions
    d_narrowbandindices = narrowbandindices
    d_pointinside = pointinside
    ! Now call the cuda function with the specified threads   
    call get_signed_distance_cuda<<<num_blocks, threads_per_block>>>(d_xf, d_yp, d_zp, d_numrays, d_num_faces, d_numberofnarrowpoints, d_faces, d_vertices, d_directions, d_narrowbandindices, d_pointinside)
    ! Copy only the results back
    allocate(finaloutput(nx,ny,nz))
    finaloutput = d_pointinside
#else
    ! Now compute the ray intersection only for the narrowband points
    print *, "*** 3rd - Pass: Computing the sign for the distance | Computed on the CPU ***" 
    call get_signed_distance(xf,yp,zp,num_faces,numberofnarrowpoints,faces,vertices,directions,narrowbandindices,pointinside)    
    ! Assign the results to the filewrite array
    allocate(finaloutput(nx,ny,nz))
    finaloutput = pointinside
#endif
    ! Log CPU time at the end of the 
    call cpu_time(time2)
    print *, "*** Signed-Distance-Field (SDF) computed in ",int(time2-time1), "seconds ***"
    elapsedTime = elapsedTime + (time2-time1)
    print *, "*** Writing SDF to file ***"
    call cpu_time(time1)
    ! Write the mask to file
    open(unit=10, file='data/mask.bin', status='replace', form='unformatted')
    write(10) finaloutput
    close(10)
    elapsedTime = elapsedTime + (time2-time1)
    print *, "*** Finished in ",int(elapsedTime), "seconds ***"

end program sdfgenerator
