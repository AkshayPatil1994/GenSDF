program generatesdf
    !
    ! Signed-Distance-Field Generator 
    !
    ! Subroutines from utils
    use utils, only : printlogo, read_inputfile, read_cans_grid, setup_grid_spacing, &
                      read_obj, getbbox, tagminmax, compute_scalar_distance_face, write2binary, &
                      fill_internal, gather_array, estimated_memoryusage
    ! Data from utils
    use utils, only : dp, inputfilename, nx, ny, nz, lx, ly, lz, r0, ng, non_uniform_grid, &
                      xp, yp, zp, xf, yf, zf, dx, dx_inverse, dy, dy_inverse, dz, dz_inverse, &
                      buffer_points, scalarvalue    
    ! MPI module
    use mpi
    ! Data from utils
    use utils, only : nprocs

    implicit none
    
    ! Input geometry related data 
    integer :: nfaces, nvertices, nnormals                          ! Number of faces, vertices, and normals
    real(dp), allocatable, dimension(:,:) :: vertices, normals      ! Array that stores vertices & vertex normals information
    integer, allocatable, dimension(:,:) :: faces, face_normals     ! Array that stroes the face (vertex) and face_normal (vertex normals) ID
    real(dp), dimension(3) :: bbox_min, bbox_max                    ! Bounding box of the geometry
    integer :: sx, ex, sy, ey, sz, ez                               ! Tagged min-max indices on the grid
    ! Signed-Distance-Field array 
    real(dp), allocatable, dimension(:,:,:) :: sdf                  ! Signed distance field array  
    ! Flood-Fill Temporary array
    real(dp), allocatable, dimension(:,:,:) :: flood_fill          ! Temporary Flood-Fill array on Processor - 0
    ! MPI data
    integer :: myid, ierror, mpi_dx, amode, fh
    integer, allocatable, dimension(:) :: decomp_x_start, decomp_x_end, decomp_size
    integer, parameter :: MPI_REAL_DP = MPI_DOUBLE_PRECISION
    ! -- AUXILIARY DATA -- !
    real(dp) :: startTime, endTime, totalTime
    real(dp) :: time1, time2
    integer :: ii, jj, kk
    logical :: debug = .False.     
!-------------------------------------------------------------------------------------------!
!                                   PROGRAM BEGINS                                          !
!-------------------------------------------------------------------------------------------!
    ! Start the MPI communication world
    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierror)    
    ! Allocate decomp_x array
    allocate(decomp_x_start(0:nprocs-1),decomp_x_end(0:nprocs-1),decomp_size(0:nprocs-1))       
    ! Print logo for the program
    if (myid == 0) then
        ! Log CPU time at the start of the program
        call cpu_time(startTime)
        call printlogo()        
        print *, "*** Starting with ",nprocs, "MPI ranks ***"
        if(nprocs == 1) then
            error stop ": Please use more than 1 processors to run the code!"        
        endif
    endif
    ! Add communication barrier for sanity check
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    ! Input data  (all MPI ranks read input file)  
    call read_inputfile(myid)
    ! Read the grid (all MPI ranks read the grid)
    call read_cans_grid(myid,'data/', dp, ng, r0, non_uniform_grid, xp, yp, zp, xf, yf, zf, dz,lx,ly,lz)
    ! Setup the grid spacing and other grid parameters
    allocate(dz_inverse(ng(3)))
    call setup_grid_spacing(myid,xf,yp,zp,ng(3),dx,dy,dz,dx_inverse,dy_inverse,dz_inverse)
    ! Read the OBJ file with the normals information
    call read_obj(myid,inputfilename,vertices,normals,faces,face_normals,nvertices,nnormals,nfaces)
    ! Query the bounding box of the geometry
    call getbbox(myid,vertices, nvertices, bbox_min, bbox_max)
    ! Tag the min-max bounding box indices on the grid
    call tagminmax(myid,xf,yp,zp,bbox_min,bbox_max,nx,ny,nz,dx,dy,dz(2),buffer_points,sx,ex,sy,ey,sz,ez)
    ! Domain decomposition based on the extent of the bounding box (in x)
    mpi_dx = ceiling(real((ex-sx+1), kind=dp) / real(nprocs, kind=dp))  ! Compute chunk size
    decomp_x_start(0) = sx
    decomp_x_end(0) = min(sx + mpi_dx - 1, ex)                        ! Ensure within bounds
    decomp_size(0) = decomp_x_end(0) - decomp_x_start(0) + 1
    
    do ii = 1, nprocs-2
        decomp_x_start(ii) = decomp_x_end(ii-1) + 1
        decomp_x_end(ii) = min(decomp_x_start(ii) + mpi_dx - 1, ex)   ! Ensure within bounds
        decomp_size(ii) = decomp_x_end(ii) - decomp_x_start(ii) + 1
    end do
    
    decomp_x_start(nprocs-1) = decomp_x_end(nprocs-2) + 1
    decomp_x_end(nprocs-1) = ex                                       ! Ensure last processor ends at ex
    decomp_size(nprocs-1) = decomp_x_end(nprocs-1) - decomp_x_start(nprocs-1) + 1
    ! Print the domain decomposition to screen (If DEBUG = .True.)
    if(debug) then
        if(myid == 0) print *, "*** Domain decomposition Details ***"
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        do ii=0,nprocs-1
            if(myid==ii) print *, "PID: ",myid," | xstart: ",decomp_x_start(myid)," | xend: ",decomp_x_end(myid), " | Size: ",decomp_size(myid) 
        enddo
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    endif    
    ! Assign a large - user specified positive value everywhere in the domain
    allocate(sdf(decomp_x_start(myid):decomp_x_end(myid),ny,nz))    
    sdf = scalarvalue
    ! Log CPU time at the end of pre-processing step
    if(myid == 0) then
        call cpu_time(time1)
        print *, "-- Finished pre-processing geometry in ",time1-startTime,"seconds..."
        ! Check estimated minimum memory required
        call estimated_memoryusage(nprocs,nfaces,nvertices,nx,ny,nz)
        print *, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "   
        print *, "*** Calculating the signed-distance-field | u-faces ***"
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)    
    ! Calculate the signed distance field - U FACES
    call compute_scalar_distance_face(myid,decomp_x_start(myid),decomp_x_end(myid),sy,ey,sz,ez,xf,yp,zp,nfaces,faces,face_normals,vertices,normals,buffer_points,sdf)            
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call gather_array(myid, nprocs, sdf, nx, ny, nz, decomp_x_start, decomp_x_end, decomp_size, sy, ey, sz, ez, scalarvalue, flood_fill, ierror)    
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    if (myid == 0) then          
        call fill_internal(flood_fill,size(xf),size(yp),size(zp),sx,sy,sz,ex,ey,ez,-100000.0_dp)
        call cpu_time(time1)
        print *, "*** Writing output data to file ***"    
        call write2binary('data/sdfu.bin',flood_fill)                
        call cpu_time(time2)
        print *, "-- Done with file write in ", time2-time1,"seconds... | u-faces "
        print *, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
        print *, "*** Calculating the signed-distance-field | v-faces ***"
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)  
      
    ! Calculate the signed distance field - V FACES
    sdf = scalarvalue               ! Reset the value for the sdf array on all processors to the user input scalar value
    call compute_scalar_distance_face(myid,decomp_x_start(myid),decomp_x_end(myid),sy,ey,sz,ez,xp,yf,zp,nfaces,faces,face_normals,vertices,normals,buffer_points,sdf)            
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call gather_array(myid, nprocs, sdf, nx, ny, nz, decomp_x_start, decomp_x_end, decomp_size, sy, ey, sz, ez, scalarvalue, flood_fill, ierror)    
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    if (myid == 0) then          
        call fill_internal(flood_fill,size(xp),size(yf),size(zp),sx,sy,sz,ex,ey,ez,-100000.0_dp)
        call cpu_time(time1)
        print *, "*** Writing output data to file ***"    
        call write2binary('data/sdfv.bin',flood_fill)                
        call cpu_time(time2)
        print *, "-- Done with file write in ", time2-time1,"seconds... | v-faces "
        print *, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
        print *, "*** Calculating the signed-distance-field | w-faces ***"
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)    

    ! Calculate the signed distance field - W FACES
    sdf = scalarvalue               ! Reset the value for the sdf array on all processors to the user input scalar value
    call compute_scalar_distance_face(myid,decomp_x_start(myid),decomp_x_end(myid),sy,ey,sz,ez,xp,yp,zf,nfaces,faces,face_normals,vertices,normals,buffer_points,sdf)            
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call gather_array(myid, nprocs, sdf, nx, ny, nz, decomp_x_start, decomp_x_end, decomp_size, sy, ey, sz, ez, scalarvalue, flood_fill, ierror)    
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    if (myid == 0) then          
        call fill_internal(flood_fill,size(xp),size(yp),size(zf),sx,sy,sz,ex,ey,ez,-100000.0_dp)
        call cpu_time(time1)
        print *, "*** Writing output data to file ***"    
        call write2binary('data/sdfw.bin',flood_fill)                
        call cpu_time(time2)
        print *, "-- Done with file write in ", time2-time1,"seconds... | w-faces "
        print *, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
        print *, "*** Calculating the signed-distance-field | Cell-Center ***"
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror) 

    ! Calculate the signed distance field - Cell Centers
    sdf = scalarvalue               ! Reset the value for the sdf array on all processors to the user input scalar value
    call compute_scalar_distance_face(myid,decomp_x_start(myid),decomp_x_end(myid),sy,ey,sz,ez,xp,yp,zp,nfaces,faces,face_normals,vertices,normals,buffer_points,sdf)            
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call gather_array(myid, nprocs, sdf, nx, ny, nz, decomp_x_start, decomp_x_end, decomp_size, sy, ey, sz, ez, scalarvalue, flood_fill, ierror)    
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    if (myid == 0) then          
        call fill_internal(flood_fill,size(xp),size(yp),size(zf),sx,sy,sz,ex,ey,ez,-100000.0_dp)
        call cpu_time(time1)
        print *, "*** Writing output data to file ***"    
        call write2binary('data/sdfp.bin',flood_fill)                
        call cpu_time(time2)
        print *, "-- Done with file write in ", time2-time1,"seconds... | Cell-Center "
        call cpu_time(endTime)
        totalTime = totalTime + (endTime-startTime)
        print *, "*** Calculation for SDF completed in ",totalTime, "seconds ***"
    endif
    ! Finalise the MPI communication
    call MPI_FINALIZE(ierror)
end program generatesdf
