program generatesdf
    !
    ! Signed-Distance-Field Generator 
    !
    ! Subroutines from utils
    use utils, only : printlogo, read_inputfile, read_cans_grid, setup_grid_spacing, &
                      read_obj, getbbox, tagminmax, compute_scalar_distance_face, write2binary
    ! Data from utils
    use utils, only : dp, inputfilename, nx, ny, nz, lx, ly, lz, r0, ng, non_uniform_grid, &
                      xp, yp, zp, xf, yf, zf, dx, dx_inverse, dy, dy_inverse, dz, dz_inverse, &
                      buffer_points, scalarvalue                    
    implicit none
    
    ! Input geometry related data (Host i.e., CPU)
    integer :: nfaces, nvertices, nnormals                          ! Number of faces, vertices, and normals
    real(dp), allocatable, dimension(:,:) :: vertices, normals      ! Array that stores vertices & vertex normals information
    integer, allocatable, dimension(:,:) :: faces, face_normals     ! Array that stroes the face (vertex) and face_normal (vertex normals) ID
    real(dp), dimension(3) :: bbox_min, bbox_max                    ! Bounding box of the geometry
    integer :: sx, ex, sy, ey, sz, ez                               ! Tagged min-max indices on the grid
    ! Signed-Distance-Field array (Host i.e., CPU)
    real(dp), allocatable, dimension(:,:,:) :: sdf  ! Signed distance field array
    
    ! -- AUXILIARY DATA -- !
    real(dp) :: startTime, endTime, totalTime
    real(dp) :: time1, time2
!-------------------------------------------------------------------------------------------!
!                                   PROGRAM BEGINS                                          !
!-------------------------------------------------------------------------------------------!
    ! Log CPU time at the start of the program
    call cpu_time(startTime)
    ! Print logo for the program
    call printlogo()
    ! Input data    
    call read_inputfile()
    ! Read the grid
    call read_cans_grid('data/', dp, ng, r0, non_uniform_grid, xp, yp, zp, xf, yf, zf, dz,lx,ly,lz)
    ! Setup the grid spacing and other grid parameters
    allocate(dz_inverse(ng(3)))
    call setup_grid_spacing(xf,yp,zp,ng(3),dx,dy,dz,dx_inverse,dy_inverse,dz_inverse)
    ! Allocate memory for the SDF field
    allocate(sdf(nx,ny,nz))
    ! Read the OBJ file with the normals information
    call read_obj(inputfilename,vertices,normals,faces,face_normals,nvertices,nnormals,nfaces)
    ! Query the bounding box of the geometry
    call getbbox(vertices, nvertices, bbox_min, bbox_max)
    ! Tag the min-max bounding box indices on the grid
    call tagminmax(xf,yp,zp,bbox_min,bbox_max,nx,ny,nz,dx,dy,dz(2),buffer_points,sx,ex,sy,ey,sz,ez)
    ! Assign a large - user specified positive value everywhere in the domain
    sdf = scalarvalue
    ! Log CPU time at the end of pre-processing step
    call cpu_time(time1)
    print *, "-- Finished pre-processing geometry in ",time1-startTime,"seconds..."
    ! Compute the signed-distance-field
    print *, "*** Calculating the signed-distance-field ***"
    call compute_scalar_distance_face(sz,ez,xf,yp,zp,nfaces,faces,face_normals,vertices,normals,buffer_points,sdf)
    ! Write file to binary format
    call cpu_time(time1)
    print *, "*** Writing output data to file ***"
    call write2binary('data/sdf.bin',sdf)
    call cpu_time(time2)
    print *, "-- Done with file write in ", time2-time1,"seconds..."
    ! Log CPU time at the end of the program
    call cpu_time(endTime)
    totalTime = totalTime + (endTime-startTime)
    print *, "*** Calculation for SDF completed in ",totalTime, "seconds ***"
    
end program generatesdf
