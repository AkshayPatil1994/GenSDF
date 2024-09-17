program sdfgenerator

    use utils

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
    integer :: v1, v2, v3, intersection_count
    real(dp), dimension(3) :: querypoint, controlpoint
    logical :: intersect
    real(dp) :: t, u, v 
    ! Data for distance
    real(dp) :: setsign, sdfresolution
    real(dp) :: distance_v1, testdata, previoustestdata    
    ! Iterator
    integer :: i, j, k, faceid, cid, point_tracker, vertexid, pbarwidth
    ! Time loggers
    real(dp) :: time1, time2, elapsedTime
    !
    ! PROGRAM BEGINS
    !
    ! Log CPU time at the start
    call cpu_time(time1)
    ! Read the input file
    call read_inputfile()
    ! Read the grid
    call read_cans_grid()
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
    call tagminmax(x,y,z,bbox_min,bbox_max,Nx,Ny,Nz,dx,dy,dz(2),sx,ex,sy,ey,sz,ez)
    ! Now allocate the narrow band indices search array
    allocate(narrowmask(nx,ny,nz))
    ! Generate the coorindate directions
    call generate_directions(directions)
    ! Log CPU time at the end of pre-processing
    call cpu_time(time2)
    elapsedTime = (time2-time1)
    print *, "- - - - - - - - Finished Pre-Processing in ", elapsedTime, "seconds..."
    ! Set a required value outside the bounding box
    pointinside = 100.0
    sdfresolution = 5.0
    pbarwidth = 30
    ! Loop over all data 
    print *, "*** 1st - Pass: Distance calculation ***"
    call cpu_time(time1)        ! Log CPU time at the start of the operation
    ! First compute the scalar distance
    do k=sz,ez
        call show_progress((k-sz),(ez-sz),pbarwidth)
        do j=sy,ey
            do i=sx,ex
                querypoint = (/x(i),y(j),z(k)/)
                previoustestdata = 1e10
                testdata = 1e10
                do vertexid=1,num_vertices
                    ! Each face has three vertices
                    ! Compute the distance
                    distance_v1 = sqrt( (x(i) - vertices(1,vertexid))*(x(i) - vertices(1,vertexid)) &
                                    +   (y(j) - vertices(2,vertexid))*(y(j) - vertices(2,vertexid)) &
                                    +   (z(k) - vertices(3,vertexid))*(z(k) - vertices(3,vertexid)) & 
                                        )                
                    testdata = min(distance_v1,previoustestdata)
                    ! Assign current testdata to previoustestdata
                    previoustestdata = testdata
                end do
                pointinside(i,j,k) = testdata                
            end do
        end do            
    end do   

    print *, "*** 2nd - Pass: Narrow band tagging"
    ! Check indices for where values are < sdfresolution*min(dx,dy,dz)
    narrowmask = pointinside < sdfresolution*min(dx,dy,dz(2))
    numberofnarrowpoints = count(narrowmask)
    print *, "Total number of narrow band point: ", numberofnarrowpoints, "of", (ex-sx)*(ey-sy)*(ez-sz), "total points."
    ! Allocate the values and indices to the narrowbandindices array
    allocate(narrowbandindices(3,numberofnarrowpoints))
    point_tracker = 1
    do k=sz,ez
        call show_progress((k-sz),(ez-sz),pbarwidth)
        do j=sy,ey
            do i=sx,ex
                if(narrowmask(i,j,k) .eqv. .True.) then
                    narrowbandindices(:,point_tracker) = (/i,j,k/)
                    point_tracker = point_tracker + 1
                endif
            end do
        end do
    end do
    ! Now compute the ray intersection only for the narrowband points
    print *, "*** 3rd - Pass: Computing the sign for the distance ***"
    ! Reset 
    do point_tracker=1,numberofnarrowpoints
        call show_progress(point_tracker,numberofnarrowpoints,pbarwidth)
        querypoint = narrowbandindices(:,point_tracker)
        intersection_count = 0
        previoustestdata = 1e10
        testdata = 1e10
        do cid=1,numberofrays                                       
            do faceid=1,num_faces
                ! Each face has three vertices
                v1 = faces(1,faceid)
                v2 = faces(2,faceid)
                v3 = faces(3,faceid)
                call ray_triangle_intersection(querypoint, directions(:,cid), vertices(:,v1), vertices(:,v2), vertices(:,v3), t, u, v, intersect)        
                if (intersect) then
                    intersection_count = intersection_count + 1
                end if        
            end do
        end do
        if(mod(intersection_count,2)>0) then
            setsign = -1.0
        else
            setsign = 1.0
        endif
        pointinside(i,j,k) = setsign*pointinside(i,j,k)                          
    end do
    ! Log CPU time at the end of the 
    call cpu_time(time2)
    print *, "*** Signed-Distance-Field (SDF) computed in ",time2-time1, "seconds ***"
    elapsedTime = elapsedTime + (time2-time1)
    print *, "*** Writing SDF to file ***"
    call cpu_time(time1)
    ! Write the mask to file
    open(unit=10, file='mask.bin', status='replace', form='unformatted')
    write(10) pointinside
    close(10)
    ! Write debugging index files
    open(unit=10, file='mask_index.bin', status='replace', form='unformatted')
    write(10) narrowbandindices
    close(10)
    call cpu_time(time2)
    elapsedTime = elapsedTime + (time2-time1)
    print *, "*** Finished in ",elapsedTime, "seconds ***"

end program sdfgenerator
