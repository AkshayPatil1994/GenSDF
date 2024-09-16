program sdfgenerator

    use utils

    implicit none

    ! Input geometry related data
    real(dp), allocatable, dimension(:,:) :: vertices
    integer, allocatable, dimension(:,:) :: faces
    integer :: num_faces, num_vertices
    
    ! Bounding related data
    real(dp), dimension(3) :: bbox_min, bbox_max
    ! Cartesian grid
    real(dp), allocatable :: x(:), y(:), z(:)
    real(dp) :: dx, dy, dz
    ! Min-Max bounds for x, y, and z
    integer :: sx, ex, sy, ey, sz, ez
    ! Data structure for masking
    real(dp), dimension(3, 19) :: directions
    real(dp), allocatable, dimension(:,:,:) :: pointinside
    integer :: v1, v2, v3, intersection_count
    real(dp), dimension(3) :: querypoint, controlpoint
    logical :: intersect
    real(dp) :: t, u, v 
    ! Data for distance
    real(dp) :: setsign
    real(dp) :: distance_v1, distance_v2, distance_v3, testdata, previoustestdata    
    ! Iterator
    integer :: i, j, k, faceid, cid
    ! Time loggers
    real(dp) :: time1, time2, elapsedTime
    !
    ! PROGRAM BEGINS
    !
    ! Log CPU time at the start
    call cpu_time(time1)
    ! Define the input file name
    inputfilename = 'sphere.obj'
    ! Define the grid points
    Nx = 512
    Ny = 256
    Nz = 128
    Lx = 30.0
    Ly = 10.0
    Lz = 8.0
    ! Allocate memory for the grid
    allocate(x(Nx),y(Ny),z(Nz))
    allocate(pointinside(Nx,Ny,Nz))
    print *, "*** Setting up grid ***"
    dx = Lx/Nx
    dy = Ly/Ny
    dz = Lz/Nz
    x(1) = 0.0
    do i=2,Nx
        x(i) = x(i-1) + dx
    end do
    y(1) = 0.0
    do i=2,Ny
        y(i) = y(i-1) + dy 
    end do
    z(1) = 0.0
    do i=2,Nz
        z(i) = z(i-1) + dz 
    end do
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
    call tagminmax(x,y,z,bbox_min,bbox_max,Nx,Ny,Nz,dx,dy,dz,sx,ex,sy,ey,sz,ez)
    ! Generate the coorindate directions
    call generate_directions(directions)
    ! Log CPU time at the end of pre-processing
    call cpu_time(time2)
    elapsedTime = (time2-time1)
    print *, "- - - - - - - - Finished Pre-Processing in ", elapsedTime, "seconds..."
    ! Set a required value outside the bounding box
    pointinside = 100.0
    ! Loop over all data 
    call cpu_time(time1)        ! Log CPU time at the start of the operation
    do k=sz,ez
        call show_progress((k-sz),(ez-sz),80)
        do j=sy,ey
            do i=sx,ex
                querypoint = (/x(i),y(j),z(k)/)
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
                        ! Compute the distance
                        distance_v1 = sqrt( (x(i) - vertices(1,v1))*(x(i) - vertices(1,v1)) &
                                        +   (y(j) - vertices(2,v1))*(y(j) - vertices(2,v1)) &
                                        +   (z(k) - vertices(3,v1))*(z(k) - vertices(3,v1)) & 
                                          )
                        distance_v2 = sqrt( (x(i) - vertices(1,v2))*(x(i) - vertices(1,v2)) &
                                        +   (y(j) - vertices(2,v2))*(y(j) - vertices(2,v2)) &
                                        +   (z(k) - vertices(3,v2))*(z(k) - vertices(3,v2)) & 
                                          )
                        distance_v3 = sqrt( (x(i) - vertices(1,v3))*(x(i) - vertices(1,v3)) &
                                        +   (y(j) - vertices(2,v3))*(y(j) - vertices(2,v3)) &
                                        +   (z(k) - vertices(3,v3))*(z(k) - vertices(3,v3)) & 
                                          )
                        testdata = min(distance_v1,distance_v2,distance_v3,previoustestdata)
                        ! Assign current testdata to previoustestdata
                        previoustestdata = testdata
                    end do
                end do
                if(mod(intersection_count,2)>0) then
                    setsign = -1.0
                else
                    setsign = 1.0
                endif
                pointinside(i,j,k) = setsign*testdata                
            end do
        end do            
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
    call cpu_time(time2)
    elapsedTime = elapsedTime + (time2-time1)
    print *, "*** Finished in ",elapsedTime, "seconds ***"

end program sdfgenerator
