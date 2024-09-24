module utils

    implicit none
    ! Real data type selector
    integer, parameter, public :: dp = selected_real_kind(15,307), &
                                  sp = selected_real_kind(6 , 37)
    !                              
    ! Input data 
    !                              
    ! -- Geometry data --
    character(len=512),  public :: inputfilename
    ! -- Cartesian grid data --
    integer, public :: nx, ny, nz, numberofrays
    real(dp), public :: lx, ly, lz, dx, dy, dx_inverse, dy_inverse
    real(dp), allocatable, dimension(:), public :: dz, dz_inverse
    real(dp), allocatable, dimension(:), public :: xp, yp, zp, xf, yf, zf
    real(dp), dimension(3), public :: r0
    integer, dimension(3), public :: ng
    logical, public :: non_uniform_grid
    ! -- SDF data --
    real(dp), public :: scalarvalue
    integer, public :: sdfresolution
    ! -- Auxiliary data --
    integer, public :: pbarwidth
    ! -- GPU data --
    integer, public :: gpu_threads

contains

  subroutine read_inputfile()
    !
    ! This subroutine reads the input file
    !
    implicit none
    ! Local variables
    integer :: iunit, ierr              ! File I/O variables
    character(len=512) :: dummyline     ! dummyline to read help info in the inputfile
    
    iunit = 10
    open(newunit=iunit,file='parameters.in',status='old',action='read',iostat=ierr)
      if( ierr == 0 ) then
        read(iunit,*,iostat=ierr) dummyline
        read(iunit,*,iostat=ierr) inputfilename
        read(iunit,*,iostat=ierr) dummyline
        read(iunit,*,iostat=ierr) numberofrays
        read(iunit,*,iostat=ierr) dummyline
        read(iunit,*,iostat=ierr) scalarvalue, sdfresolution, pbarwidth
        read(iunit,*,iostat=ierr) dummyline
        read(iunit,*,iostat=ierr) nx, ny, nz
        read(iunit,*,iostat=ierr) dummyline
        read(iunit,*,iostat=ierr) r0(1), r0(2), r0(3)
        read(iunit,*,iostat=ierr) dummyline
        read(iunit,*,iostat=ierr) non_uniform_grid
        read(iunit,*,iostat=ierr) dummyline
        read(iunit,*,iostat=ierr) gpu_threads
        ! Force numberofrays is odd
        if(mod(numberofrays,2)==0) then
          numberofrays = numberofrays + 1
          if(numberofrays > 19) then
            numberofrays = 19
            print *, "Number of rays forced within limit to", numberofrays
          else
            print *, "Number of rays edited:", numberofrays
          endif
        endif

      else
        error stop "parameters.in file encountered a problem!"
      end if
  end subroutine read_inputfile

  subroutine read_cans_grid(loc, iprecision, npoints, origin, non_uni_grid, xin_p, yin_p, zin_p, xin_f, yin_f, zin_f, dzin)
    !
    ! This subroutine reads the CaNS grid to memory
    !
    ! Input
    character(len=*), intent(in) :: loc                   ! Directory location of the grid file
    integer, intent(in) :: iprecision                     ! Precision level (4 (single) or 8 (double))
    real(dp), dimension(3), intent(in) :: origin          ! Location of the origin
    logical, intent(in) :: non_uni_grid                   ! Flag for non-uniform grid
    integer, dimension(3), intent(out) :: npoints         ! Grid size in x, y, z (modified by geometry file)
    ! Output
    real(dp), intent(out), allocatable, dimension(:) :: xin_p, yin_p, zin_p   ! Centered grid arrays
    real(dp), intent(out), allocatable, dimension(:) :: xin_f, yin_f, zin_f   ! Staggered grid arrays
    real(dp), intent(out), allocatable, dimension(:) :: dzin                  ! Grid spacing in z
    ! Local variables
    logical :: fexists                                    ! File exists boolean
    integer :: iter                                       ! Local Iterator
    real(dp) :: dl(3), l(3)                               ! Grid spacing and length
    real(sp), allocatable :: grid_z4(:,:)                 ! For non-uniform grid, single precision
    real(dp), allocatable :: grid_z8(:,:)                 ! For non-uniform grid, double precision
    character(len=512) :: geofile, grdfile                ! Filenames for geometry.out and grid.out

    ! Check the file directory
    inquire(file=trim(loc), exist=fexists)
    if (fexists .eqv. .False.) then
        print*, "The input directory does not exist: ", trim(loc)
        stop
    endif

    ! File paths
    geofile = trim(loc) // "geometry.out"

    ! Read geometry file
    open(unit=10, file=geofile, status='old')
    read(10, *) npoints      ! Read the grid size (x, y, z)
    read(10, *) l       ! Read the domain lengths
    close(10)

    ! Calculate grid spacing
    dl = l / real(npoints, 8)
    dx = dl(1)
    dx_inverse = 1.0_dp/dx
    dy = dl(2)
    dy_inverse = 1.0_dp/dy

    ! Generate centered grids
    allocate(xin_p(npoints(1)), yin_p(npoints(2)), zin_p(npoints(3)))
    allocate(xin_f(npoints(1)), yin_f(npoints(2)), zin_f(npoints(3)))
    allocate(dzin(npoints(3)),dz_inverse(npoints(3)))

    do iter = 1, npoints(1)
        xin_p(iter) = origin(1) + (iter - 0.5) * dl(1)  ! Centered x grid
        xin_f(iter) = xin_p(iter) + dl(1) / 2.0        ! Staggered x grid
    end do

    do iter = 1, npoints(2)
        yin_p(iter) = origin(2) + (iter - 0.5) * dl(2)  ! Centered y grid
        yin_f(iter) = yin_p(iter) + dl(2) / 2.0        ! Staggered y grid
    end do

    do iter = 1, npoints(3)
        zin_p(iter) = origin(3) + (iter - 0.5) * dl(3)  ! Centered z grid
        zin_f(iter) = zin_p(iter) + dl(3) / 2.0        ! Staggered z grid
    end do
    ! Compute dz
    dzin(1) = zin_f(1)
    dz_inverse(1) = 1.0_dp/dzin(1)
    do iter = 2,npoints(3)
      dzin(iter) = zin_f(iter) - zin_f(iter-1)
      dz_inverse(iter) = 1.0_dp/(zin_f(iter) - zin_f(iter-1))
    end do

    ! Non-uniform grid handling
    if (non_uni_grid) then
        grdfile = trim(loc) // "grid.bin"

        if (iprecision == 4) then
          open(unit=20, file=grdfile, form="unformatted", access="stream")
          allocate(grid_z4(npoints(3), 4))
          read(20) grid_z4  ! Read the grid_z binary file
          close(20)
          zin_p = origin(3) + grid_z4(:,2)
          zin_f = origin(3) + grid_z4(:,3)
        else if (iprecision == 8) then
            open(unit=20, file=grdfile, form="unformatted", access="stream")
            allocate(grid_z8(npoints(3), 4))
            read(20) grid_z8  ! Read the grid_z binary file
            close(20)
            zin_p = origin(3) + grid_z8(:,2)
            zin_f = origin(3) + grid_z8(:,3)
        endif
    endif


  end subroutine read_cans_grid

  subroutine read_obj(filename, vertices_in, faces_in, num_vertices_in, num_faces_in)
    !
    ! This subroutine reads an OBJ file and returns the vertices, faces, number of vertices, and number of faces
    !
    implicit none
    ! Inputs
    character(len=*), intent(in) :: filename  ! Name of the OBJ file
    ! Outputs
    real(dp), allocatable, intent(out) :: vertices_in(:,:)     ! 2D array to store vertex coordinates (3 x num_vertices)
    integer, allocatable, intent(out) :: faces_in(:,:)         ! 2D array to store faces (3 x num_faces)
    integer, intent(out) :: num_vertices_in                    ! Number of vertices
    integer, intent(out) :: num_faces_in                       ! Number of faces
    ! Local variables
    integer :: unit, ios, v_count, f_count
    character(len=512) :: line                              ! Line buffer to read the file
    real(dp) :: vx, vy, vz                                  ! Variables for vertex coordinates
    integer :: v1_in, v2_in, v3_in                          ! Variables for face indices
    real(dp) :: time1, time2                                ! Total CPU time spent on reading the file

    ! Log CPU time at the start
    call cpu_time(time1)
    ! Open the OBJ file
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening file: ', filename
      stop
    end if

    ! Initialize counters
    v_count = 0
    f_count = 0

    ! First pass: Count vertices and faces
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit

      ! Check the line for vertex data
      if (line(1:2) == 'v ') then
        v_count = v_count + 1
      end if

      ! Check the line for face data
      if (line(1:2) == 'f ') then
        f_count = f_count + 1
      end if
    end do

    ! Allocate arrays based on the counts
    allocate(vertices_in(3, v_count))
    allocate(faces_in(3, f_count))

    ! Store the number of vertices and faces
    num_vertices_in = v_count
    num_faces_in = f_count

    ! Rewind the file to read it again
    rewind(unit)

    ! Reset counters
    v_count = 0
    f_count = 0

    ! Second pass: Read vertex and face data
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit

      ! Read vertex data
      if (line(1:2) == 'v ') then
        read(line(3:), *) vx, vy, vz
        v_count = v_count + 1
        vertices_in(1, v_count) = vx
        vertices_in(2, v_count) = vy
        vertices_in(3, v_count) = vz
      end if

      ! Read face data (only triangular faces assumed)
      if (line(1:2) == 'f ') then
        read(line(3:), *) v1_in, v2_in, v3_in
        f_count = f_count + 1
        faces_in(1, f_count) = v1_in
        faces_in(2, f_count) = v2_in
        faces_in(3, f_count) = v3_in
      end if
    end do

    ! Close the file
    close(unit)

    ! Log CPU time at the end of the operation
    call cpu_time(time2)
    print *, "*** File: '",trim(filename), "' sucessfully read in ",time2-time1, "seconds ***"
  end subroutine read_obj

  subroutine getbbox(vertices_in, num_vertices_in, bbox_min, bbox_max)
    !
    ! This subroutine calculates the bounding box for the OBJ geometry
    !
    implicit none
    ! Inputs
    real(dp), intent(in) :: vertices_in(3, num_vertices_in)   ! 2D array to store vertex coordinates (3 x num_vertices)
    integer, intent(in) :: num_vertices_in                    ! Number of vertices  
    ! Outputs
    real(dp), intent(out) :: bbox_min(3)                ! Min values for x, y, z
    real(dp), intent(out) :: bbox_max(3)                ! Max values for x, y, z
    ! Local variables
    integer :: iter                                     ! Iterator
    real(dp) :: vx, vy, vz                              ! Vertex data
  
    ! Initialize bbox_min and bbox_max to the first vertex
    bbox_min = vertices_in(:, 1)
    bbox_max = vertices_in(:, 1)
  
    ! Iterate over all vertices to find the bounding box
    do iter = 2, num_vertices_in
      vx = vertices_in(1, iter)
      vy = vertices_in(2, iter)
      vz = vertices_in(3, iter)
  
      ! Update the minimum values for x, y, z
      if (vx < bbox_min(1)) bbox_min(1) = vx
      if (vy < bbox_min(2)) bbox_min(2) = vy
      if (vz < bbox_min(3)) bbox_min(3) = vz
  
      ! Update the maximum values for x, y, z
      if (vx > bbox_max(1)) bbox_max(1) = vx
      if (vy > bbox_max(2)) bbox_max(2) = vy
      if (vz > bbox_max(3)) bbox_max(3) = vz
    end do
  
  end subroutine getbbox

  subroutine getfacenormals(vertices_in, faces_in, num_faces_in, num_vertices_in,normals)
    !
    ! This subroutine computes the face normals for each face
    !
    implicit none  
    ! Inputs
    integer, intent(in) :: num_faces_in, num_vertices_in      ! Number of faces and vertices
    real(dp), intent(in) :: vertices_in(3,num_vertices_in)    ! Vertex coordinates (3 x num_vertices)
    integer, intent(in) :: faces_in(3,num_faces_in)           ! Face indices (3 x num_faces)  
    ! Outputs
    real(dp), allocatable, intent(out) :: normals(:,:)  ! Face normals (3 x num_faces)  
    ! Local variables
    integer :: iter, v1_in, v2_in, v3_in                ! Iterator and vertex data
    real(dp) :: e1(3), e2(3), n(3), norm                ! Face normals and edge vectors
    real(dp) :: time1, time2                            ! CPU time markers
  
    ! Log the CPU time at the start
    call cpu_time(time1)
    
    ! Allocate memory for normals
    allocate(normals(3, num_faces_in))
  
    ! Loop over each face to calculate the normal
    do iter = 1, num_faces_in
      ! Get the vertex indices for the face (OBJ indices are 1-based)
      v1_in = faces_in(1, iter)
      v2_in = faces_in(2, iter)
      v3_in = faces_in(3, iter)
  
      ! Calculate edge vectors e1 and e2
      e1(1) = vertices_in(1, v2_in) - vertices_in(1, v1_in)
      e1(2) = vertices_in(2, v2_in) - vertices_in(2, v1_in)
      e1(3) = vertices_in(3, v2_in) - vertices_in(3, v1_in)
  
      e2(1) = vertices_in(1, v3_in) - vertices_in(1, v1_in)
      e2(2) = vertices_in(2, v3_in) - vertices_in(2, v1_in)
      e2(3) = vertices_in(3, v3_in) - vertices_in(3, v1_in)
  
      ! Calculate cross product e1 x e2 to get the normal vector
      n(1) = e1(2) * e2(3) - e1(3) * e2(2)
      n(2) = e1(3) * e2(1) - e1(1) * e2(3)
      n(3) = e1(1) * e2(2) - e1(2) * e2(1)
  
      ! Calculate the magnitude (norm) of the normal vector
      norm = sqrt(n(1)**2 + n(2)**2 + n(3)**2)
  
      ! Normalize the normal vector
      if (norm > 0.0d0) then
        n(1) = n(1) / norm
        n(2) = n(2) / norm
        n(3) = n(3) / norm
      end if
  
      ! Store the normalized normal vector
      normals(1, iter) = n(1)
      normals(2, iter) = n(2)
      normals(3, iter) = n(3)
    end do
  
    ! Log the CPU time at the end 
    call cpu_time(time2)
    print *, "Normals computed in ",time2-time1,"seconds...."
  end subroutine getfacenormals

  subroutine tagminmax(x,y,z,bbox_min,bbox_max,Nx,Ny,Nz,dx,dy,dz,sx,ex,sy,ey,sz,ez)
    !
    ! This subroutine tags the location of the min and max indices based on bounding box
    !
    implicit none
    ! Input Argument
    real(dp), dimension(:), intent(in) :: x, y, z                 ! x, y, z coordinates
    real(dp), intent(in) :: dx, dy, dz                            ! grid spacing
    real(dp), dimension(:), intent(in) :: bbox_min, bbox_max      ! Bounding box min & max
    integer, intent(in) :: Nx, Ny, Nz                             ! Grid points
    ! Output
    integer, intent(out) :: sx,ex,sy,ey,sz,ez                     ! Start-end index
    ! Local variables
    integer :: i, j, k                                            ! Iterators

    ! Set min and max to domain extents
    sx = 1
    ex = nx
    sy = 1
    ey = ny
    sz = 1
    ez = nz
    
    do i=1,Nx
      ! Tag min
      if(x(i)<=bbox_min(1)-5*dx) then
          sx = i
      else if (x(i)<=bbox_max(1)+5*dx) then
          ex = i
      else
        exit 
      end if
    end do

    do j=1,Ny
      ! Tag max
      if(y(j)<=bbox_min(2)-5*dy) then
          sy = j
      else if(y(j)<=bbox_max(2)+5*dy) then
          ey = j
      else
        exit 
      end if
    end do

    do k=1,Nz
      ! Tag max
      if(z(k)<=bbox_min(3)-5*dz) then
          sz = k
      else if (z(k)<=bbox_max(3)+5*dz) then
          ez = k
      else
          exit
      end if
    end do

    ! Print to screen
    print *, "*** Min-Max Index-Value pair ***"
    print *, "Min-Max x:", sx, x(sx), "|", ex, x(ex)
    print *, "Min-Max y:", sy, y(sy), "|", ey, y(ey)
    print *, "Min-Max z:", sz, z(sz), "|", ez, z(ez)
  end subroutine tagminmax

  subroutine cross_product(result, u, v)
    !
    ! This subroutine computes the cross-product between two vectors
    !
    implicit none
    ! Input
    real(dp), dimension(3), intent(in) :: u, v            ! Input vectors
    ! Output
    real(dp), dimension(3), intent(out) :: result         ! Output vector

    result(1) = u(2) * v(3) - u(3) * v(2)
    result(2) = u(3) * v(1) - u(1) * v(3)
    result(3) = u(1) * v(2) - u(2) * v(1)
  end subroutine cross_product

  subroutine ray_triangle_intersection(orig, dir, v0, v1, v2, t, u, v, intersect)
    !
    ! This subroutine checks for ray-triangle intersection
    !
    implicit none
    ! Input parameters
    real(dp), dimension(3), intent(in) :: orig    ! Ray origin
    real(dp), dimension(3), intent(in) :: dir     ! Ray direction
    real(dp), dimension(3), intent(in) :: v0, v1, v2  ! Triangle vertices    
    ! Output parameters
    real(dp), intent(out) :: t, u, v              ! Intersection details
    logical, intent(out) :: intersect             ! Whether an intersection occurred    
    ! Local variables
    real(dp) :: epsilon, inv_det, det, a(3), b(3), pvec(3), tvec(3), qvec(3)
    
    ! Initialize
    epsilon = 1.0e-12     ! >> This could be put in as a subroutine input argument later
    intersect = .false.   

    ! Find vectors for the two edges sharing v0
    a = v1 - v0
    b = v2 - v0
    
    ! Begin calculating the determinant
    call cross_product(pvec, dir, b)
    det = dot_product(a, pvec)
    
    ! If the determinant is close to 0, the ray lies in the plane of the triangle
    if (abs(det) < epsilon) return
    
    inv_det = 1.0 / det
    
    ! Calculate distance from v0 to ray origin
    tvec = orig - v0
    
    ! Calculate u parameter and test bounds
    u = dot_product(tvec, pvec) * inv_det
    if (u < 0.0 .or. u > 1.0) return
    
    ! Prepare to test v parameter
    call cross_product(qvec, tvec, a)
    
    ! Calculate v parameter and test bounds
    v = dot_product(dir, qvec) * inv_det
    if (v < 0.0 .or. u + v > 1.0) return
    
    ! Calculate t to determine the intersection point
    t = dot_product(b, qvec) * inv_det
    
    ! If t is positive, an intersection has occurred
    intersect = t > epsilon
  end subroutine ray_triangle_intersection

  subroutine generate_directions(directions)
    !
    ! This subroutine generates the lattice-like directions for shooting rays
    !
    implicit none
    ! Output
    real(dp), dimension(3, 19), intent(out) :: directions
    ! Local variables
    integer :: i

    ! Define the 6 axial directions (along x, y, z)
    directions(:, 1) = (/ 1.0d0,  0.0d0,  0.0d0 /)   ! +x
    directions(:, 2) = (/ -1.0d0, 0.0d0,  0.0d0 /)   ! -x
    directions(:, 3) = (/ 0.0d0,  1.0d0,  0.0d0 /)   ! +y
    directions(:, 4) = (/ 0.0d0, -1.0d0,  0.0d0 /)   ! -y
    directions(:, 5) = (/ 0.0d0,  0.0d0,  1.0d0 /)   ! +z
    directions(:, 6) = (/ 0.0d0,  0.0d0, -1.0d0 /)   ! -z

    ! Define the 12 diagonal directions (along face diagonals of the cube)
    directions(:, 7)  = (/ 1.0d0,  1.0d0,  0.0d0 /)   ! +x, +y
    directions(:, 8)  = (/ -1.0d0, 1.0d0,  0.0d0 /)   ! -x, +y
    directions(:, 9)  = (/ 1.0d0, -1.0d0,  0.0d0 /)   ! +x, -y
    directions(:, 10) = (/ -1.0d0,-1.0d0,  0.0d0 /)   ! -x, -y

    directions(:, 11) = (/ 1.0d0,  0.0d0,  1.0d0 /)   ! +x, +z
    directions(:, 12) = (/ -1.0d0, 0.0d0,  1.0d0 /)   ! -x, +z
    directions(:, 13) = (/ 1.0d0,  0.0d0, -1.0d0 /)   ! +x, -z
    directions(:, 14) = (/ -1.0d0, 0.0d0, -1.0d0 /)   ! -x, -z

    directions(:, 15) = (/ 0.0d0,  1.0d0,  1.0d0 /)   ! +y, +z
    directions(:, 16) = (/ 0.0d0, -1.0d0,  1.0d0 /)   ! -y, +z
    directions(:, 17) = (/ 0.0d0,  1.0d0, -1.0d0 /)   ! +y, -z
    directions(:, 18) = (/ 0.0d0, -1.0d0, -1.0d0 /)   ! -y, -z
    
    do i = 7, 18
      directions(:, i) = directions(:, i) / sqrt(2.0d0)  ! Normalize diagonal vectors
    end do

    directions(:, 19) = (/ 0.0d0, 0.0d0, 0.0d0 /)  ! Optional zero vector for completeness

  end subroutine generate_directions

  subroutine compute_scalar_distance(sx,ex,sy,ey,sz,ez,x,y,z,num_vertices,vertices,sdf_index_add,pointinside)
    !
    ! This subroutine computes the distance between the vertex of the geometry and 'n' query points
    ! closest to this vertex where 'n' is a user input
    !
    implicit none

    ! Input
    integer, intent(in) :: sx, ex, sy, ey, sz, ez               ! Bounding box based on geometry
    integer, intent(in) :: num_vertices, sdf_index_add          ! Number of vertices and surrounding number of grid points to query
    real(dp), intent(in), dimension(:,:) :: vertices            ! Array of vertices
    real(dp), intent(in), dimension(:) :: x, y, z               ! Vectors of cartesian grid
    ! Output
    real(dp), intent(out), dimension(:,:,:) :: pointinside      ! Array of distance without sign
    ! Local variable declarations
    integer :: k, vertexid                                      
    integer :: sidx, sidy, sidz
    integer :: ii, jj, kk
    integer :: xloc, yloc, zloc
    real(dp) :: local_distance(sdf_index_add*2,sdf_index_add*2,sdf_index_add*2)

    ! Loop over all vertices
    do vertexid = 1, num_vertices
        ! Display progress
        call show_progress(vertexid, num_vertices, 50)
        ! Calculate xloc, yloc, zloc (floor based on vertex location)
        xloc = floor(vertices(1,vertexid) * dx_inverse) - sdf_index_add
        yloc = floor(vertices(2,vertexid) * dy_inverse) - sdf_index_add
        ! Find zloc such that z(zloc) <= vertices(3,vertexid)
        k = sz
        do while (z(k) <= vertices(3,vertexid))
            zloc = k
            k = k + 1
        end do
        zloc = zloc - sdf_index_add    
        ! Calculate the local distances
        do kk = 1, sdf_index_add * 2
            do jj = 1, sdf_index_add * 2
                do ii = 1, sdf_index_add * 2
                    sidx = xloc + ii
                    sidy = yloc + jj
                    sidz = zloc + kk

                    ! Check if indices are within bounds
                    if (sidx >= sx .and. sidx <= ex .and. sidy >= sy .and. sidy <= ey .and. sidz >= sz .and. sidz <= ez) then
                        local_distance(ii,jj,kk) =  (x(sidx) - vertices(1,vertexid))*(x(sidx) - vertices(1,vertexid)) + &
                                                    (y(sidy) - vertices(2,vertexid))*(y(sidy) - vertices(2,vertexid)) + &
                                                    (z(sidz) - vertices(3,vertexid))*(z(sidz) - vertices(3,vertexid))  
                        ! Assign the distance if its smaller than previous value
                        pointinside(sidx,sidy,sidz) = min(pointinside(sidx,sidy,sidz),local_distance(ii,jj,kk))                    
                    end if
                end do
            end do
        end do
    end do

    ! Square root of the distance
    pointinside(sx:ex,sy:ey,sz:ez) = sqrt(pointinside(sx:ex,sy:ey,sz:ez))

  end subroutine compute_scalar_distance

  subroutine tag_narrowband_points(sx,ex,sy,ey,sz,ez,narrowmask,narrowbandindices)
    !
    ! This subroutine tags the narrowband points
    !
    implicit none
    ! Input
    integer, intent(in) :: sx,ex,sy,ey,sz,ez
    logical, intent(in), dimension(:,:,:) :: narrowmask
    ! Output
    integer, intent(inout), dimension(:,:) :: narrowbandindices
    ! Local variables
    integer :: i, j, k
    integer :: point_tracker

    point_tracker = 1
    do k=sz,ez
      call show_progress((k-sz),(ez-sz),50)
      do j=sy,ey
          do i=sx,ex
              if(narrowmask(i,j,k) .eqv. .True.) then
                  narrowbandindices(:,point_tracker) = (/i,j,k/)
                  point_tracker = point_tracker + 1
              endif
          end do
      end do
    end do

  end subroutine tag_narrowband_points

  subroutine get_signed_distance(x,y,z,num_faces,numberofnarrowpoints,faces,vertices,directions,narrowbandindices,pointinside)
    ! 
    ! This subroutine tags a sign for the distance that is pre-computed
    !
    implicit none
    ! Input
    real(dp), intent(in), dimension(:) :: x, y, z
    integer, intent(in) :: numberofnarrowpoints, num_faces
    integer, intent(in), dimension(:,:) :: narrowbandindices
    integer, intent(in), dimension(:,:) :: faces
    real(dp), intent(in), dimension(:,:) :: vertices
    real(dp), intent(in), dimension(:,:) :: directions 
    ! Output
    real(dp), intent(inout), dimension(:,:,:) :: pointinside
    ! Local variables
    integer :: point_tracker, cid, faceid, i, j, k
    integer :: v1, v2, v3
    integer :: intersection_count
    real(dp), dimension(3) :: querypoint
    logical :: intersect
    real(dp) :: setsign
    real(dp) :: t, u, v
  
    do point_tracker=1,numberofnarrowpoints
        call show_progress(point_tracker,numberofnarrowpoints,pbarwidth)
        i = narrowbandindices(1,point_tracker)
        j = narrowbandindices(2,point_tracker)
        k = narrowbandindices(3,point_tracker)
        querypoint = (/x(i),y(j),z(k)/)
        intersection_count = 0
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
  
  end subroutine get_signed_distance  

  subroutine show_progress(current, total, width)
    ! 
    ! This subroutine aims to display a progress bar and overwrite the same line
    !
    implicit none
    ! Input
    integer, intent(in) :: current, total, width
    ! Local variables
    integer(dp) :: progress
    real(dp) :: percent
    character(len=width) :: bar

    ! Calculate percentage and progress bar length
    percent = real(current) / real(total) * 100.0
    progress = int(real(current) / real(total) * width)

    ! Construct the progress bar string
    bar = repeat('|', progress) // repeat(' ', width - progress)

    ! Print the progress bar with a carriage return to stay on the same line
    write(*,'(a, "|", a, "|", f6.2, "%", 1x, i8, "/", i8)', advance='no') char(13), bar, percent, current, total

    ! If the current iteration is the last, move to a new line
    if (current == total) then
        print *  ! Move to the next line after the final iteration
    end if
  end subroutine show_progress
  
#if defined(_ISCUDA)  
  ! --- CUDA RELATED SUBROUTINES & FUNCTIONS --- !
  attributes(global) subroutine get_signed_distance_cuda(x, y, z, numrays, num_faces, numberofnarrowpoints, faces, vertices, directions, narrowbandindices, pointinside)
    implicit none
    ! Input
    real(dp), intent(in), dimension(:), device :: x, y, z
    integer, intent(in), device :: numberofnarrowpoints, num_faces, numrays
    integer, intent(in), dimension(:,:), device :: narrowbandindices
    integer, intent(in), dimension(:,:), device :: faces
    real(dp), intent(in), dimension(:,:), device :: vertices
    real(dp), intent(in), dimension(:,:), device :: directions 
    ! Output
    real(dp), intent(inout), dimension(:,:,:), device :: pointinside
    ! Local variables
    integer :: point_tracker, cid, faceid, i, j, k
    integer :: v1, v2, v3
    integer :: intersection_count
    real(dp), dimension(3) :: querypoint
    logical :: intersect
    real(dp) :: setsign
    real(dp) :: t, u, v

    ! Get the thread index for parallel execution
    point_tracker = blockIdx%x * blockDim%x + threadIdx%x + 1

    if (point_tracker <= numberofnarrowpoints) then
        i = narrowbandindices(1, point_tracker)
        j = narrowbandindices(2, point_tracker)
        k = narrowbandindices(3, point_tracker)
        querypoint = (/x(i), y(j), z(k)/)
        intersection_count = 0

        do cid = 1, numrays
            do faceid = 1, num_faces
                v1 = faces(1, faceid)
                v2 = faces(2, faceid)
                v3 = faces(3, faceid)
                call ray_triangle_intersection_cuda(querypoint, directions(:,cid), vertices(:,v1), vertices(:,v2), vertices(:,v3), t, u, v, intersect)
                if (intersect) then
                    intersection_count = intersection_count + 1
                end if
            end do
        end do

        if (mod(intersection_count, 2) > 0) then
            setsign = -1.0
        else
            setsign = 1.0
        endif
        pointinside(i, j, k) = setsign * pointinside(i, j, k)
    endif
  end subroutine get_signed_distance_cuda

  attributes(device)  subroutine ray_triangle_intersection_cuda(orig, dir, v0, v1, v2, t, u, v, intersect)
    !
    ! This subroutine checks for ray-triangle intersection
    !
    implicit none
    ! Input parameters
    real(dp), dimension(3), intent(in), device :: orig        ! Ray origin
    real(dp), dimension(3), intent(in), device :: dir         ! Ray direction
    real(dp), dimension(3), intent(in), device :: v0, v1, v2  ! Triangle vertices    
    ! Output parameters
    real(dp), intent(out), device :: t, u, v              ! Intersection details
    logical, intent(out), device :: intersect             ! Whether an intersection occurred    
    ! Local variables
    real(dp) :: epsilon, inv_det, det, a(3), b(3), pvec(3), tvec(3), qvec(3)
    
    ! Initialize
    epsilon = 1.0e-12     ! >> This could be put in as a subroutine input argument later
    intersect = .false.   

    ! Find vectors for the two edges sharing v0
    a = v1 - v0
    b = v2 - v0
    
    ! Begin calculating the determinant
    call cross_product_cuda(pvec, dir, b)
    det = dot_product_cuda(a, pvec, 3)
    
    ! If the determinant is close to 0, the ray lies in the plane of the triangle
    if (abs(det) < epsilon) return
    
    inv_det = 1.0 / det
    
    ! Calculate distance from v0 to ray origin
    tvec = orig - v0
    
    ! Calculate u parameter and test bounds
    u = dot_product_cuda(tvec, pvec, 3) * inv_det
    if (u < 0.0 .or. u > 1.0) return
    
    ! Prepare to test v parameter
    call cross_product_cuda(qvec, tvec, a)
    
    ! Calculate v parameter and test bounds
    v = dot_product_cuda(dir, qvec, 3) * inv_det
    if (v < 0.0 .or. u + v > 1.0) return
    
    ! Calculate t to determine the intersection point
    t = dot_product_cuda(b, qvec, 3) * inv_det
    
    ! If t is positive, an intersection has occurred
    intersect = t > epsilon

  end subroutine ray_triangle_intersection_cuda

  attributes(device) subroutine cross_product_cuda(result, u, v)
    !
    ! This subroutine computes the cross-product between two vectors
    !
    implicit none
    ! Input
    real(dp), dimension(3), intent(in), device :: u, v            ! Input vectors
    ! Output
    real(dp), dimension(3), intent(out), device :: result         ! Output vector

    result(1) = u(2) * v(3) - u(3) * v(2)
    result(2) = u(3) * v(1) - u(1) * v(3)
    result(3) = u(1) * v(2) - u(2) * v(1)
  end subroutine cross_product_cuda

  attributes(device) function dot_product_cuda(vec1, vec2, n) result(outvec)
    implicit none
    integer, intent(in), device :: n              ! Length of the vectors
    real(dp), intent(in), device :: vec1(n)       ! First vector
    real(dp), intent(in), device :: vec2(n)       ! Second vector
    real(dp), device :: outvec                    ! Result of the dot product

    integer :: i

    outvec = 0.0_dp

    ! Compute the dot product
    do i = 1, n
      outvec = outvec + vec1(i) * vec2(i)
    end do
  end function dot_product_cuda
#endif

end module utils
