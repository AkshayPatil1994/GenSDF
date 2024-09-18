program sdfgenerator
    !
    ! Signed-Distance-Field Generator 
    !
    use mod_mpi, only :  myid, ierr, size
    use utils, only : read_inputfile, read_cans_grid, read_obj, getbbox, &
                      tagminmax, generate_directions, compute_scalar_distance, &
                      tag_narrowband_points, get_signed_distance 
    use utils, only : inputfilename, nx, ny, nz, numberofrays, lx, ly, lz, dx, dy, dz, &
                      sdfresolution, scalarvalue, pbarwidth, dp, sp, r0, ng, non_uniform_grid, &
                      xp, yp, zp, xf, yf, zf
    use mpi 

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
    real(dp), allocatable, dimension(:,:,:) :: pointinside
    ! Time loggers
    real(dp) :: time1, time2, elapsedTime
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
    ! Use bounding box extents plus a small value (order 5*grid size) as control points
    ! Note: For z the buffer is 5*smallest grid point
    call tagminmax(xf,yp,zp,bbox_min,bbox_max,Nx,Ny,Nz,dx,dy,dz(2),sx,ex,sy,ey,sz,ez)        
    ! Now allocate the narrow band indices search array
    allocate(narrowmask(nx,ny,nz))
    ! Generate the coorindate directions
    call generate_directions(directions)
    ! Log CPU time at the end of pre-processing
    call cpu_time(time2)
    elapsedTime = (time2-time1)
    print *, "- - - - - - - - Finished Pre-Processing in ", elapsedTime, "seconds..."
    ! Set a required value outside the bounding box (User specified)
    pointinside = scalarvalue
    ! Loop over all data 
    print *, "*** 1st - Pass: Distance calculation ***"
    call cpu_time(time1)        ! Log CPU time at the start of the operation
    ! First compute the scalar distance
    call compute_scalar_distance(sx,ex,sy,ey,sz,ez,xf,yp,zp,num_vertices,vertices,pointinside)  
    print *, "*** 2nd - Pass: Narrow band tagging ***"
    ! Check indices for where values are < sdfresolution*min(dx,dy,dz)
    narrowmask = pointinside < sdfresolution*min(dx,dy,dz(2))
    numberofnarrowpoints = count(narrowmask)
    print *, "Narrow band point: ", numberofnarrowpoints, "of", (ex-sx)*(ey-sy)*(ez-sz), "total points."
    ! Allocate the values and indices to the narrowbandindices array
    allocate(narrowbandindices(3,numberofnarrowpoints))
    call tag_narrowband_points(sx,ex,sy,ey,sz,ez,narrowmask,narrowbandindices)
    ! Now compute the ray intersection only for the narrowband points
    print *, "*** 3rd - Pass: Computing the sign for the distance ***" 
    call get_signed_distance(xf,yp,zp,num_faces,numberofnarrowpoints,faces,vertices,directions,narrowbandindices,pointinside)
    ! Log CPU time at the end of the 
    call cpu_time(time2)
    print *, "*** Signed-Distance-Field (SDF) computed in ",time2-time1, "seconds ***"
    elapsedTime = elapsedTime + (time2-time1)
    print *, "*** Writing SDF to file ***"
    call cpu_time(time1)
    ! Write the mask to file
    open(unit=10, file='data/mask.bin', status='replace', form='unformatted')
    write(10) pointinside
    close(10)
    elapsedTime = elapsedTime + (time2-time1)
    print *, "*** Finished in ",elapsedTime, "seconds ***"

end program sdfgenerator
