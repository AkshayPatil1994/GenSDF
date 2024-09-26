module utils

    implicit none
    ! Real data type selector
    integer, parameter, public :: dp = selected_real_kind(15,307), &
                                  sp = selected_real_kind(6 , 37)
    !                              
    ! Input data [Public data declared first then the rest]
    !                              
    ! -- Geometry data --
    character(len=512), public :: inputfilename
    ! -- Cartesian grid data --
    integer, public :: numberofrays
    integer, public :: nx, ny, nz
    real(dp), dimension(3), public :: r0
    integer, dimension(3), public :: ng
    logical, public :: non_uniform_grid
    real(dp) :: lx, ly, lz, dx, dy, dx_inverse, dy_inverse
    real(dp), allocatable, dimension(:) :: dz, dz_inverse
    real(dp), allocatable, dimension(:) :: xp, yp, zp, xf, yf, zf
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

    ! Generate centered grids
    allocate(xin_p(npoints(1)), yin_p(npoints(2)), zin_p(npoints(3)))
    allocate(xin_f(npoints(1)), yin_f(npoints(2)), zin_f(npoints(3)))
    allocate(dzin(npoints(3)))

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
        zin_f(iter) = zin_p(iter) + dl(3) / 2.0         ! Staggered z grid
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

  subroutine setup_grid_spacing(xin,yin,zin,nz_in,dx_out,dy_out,dz_out,dx_inv_out,dy_inv_out,dz_inv_out)
    !
    ! This subroutine setups up the grid spacing and its inverse
    !
    implicit none
    ! Input
    real(dp), intent(in), dimension(:) :: xin, yin, zin
    integer, intent(in) :: nz_in
    ! Output
    real(dp), intent(out) :: dx_out, dy_out, dx_inv_out, dy_inv_out
    real(dp), intent(out), dimension(:) :: dz_out, dz_inv_out
    ! Local variable
    integer :: ii
    
    ! Compute grid spacing below
    dx_out = xin(10) - xin(9)
    dy_out = yin(10) - yin(9)
    ! Inverse of grid spacing
    dx_inv_out = 1.0_dp/dx_out
    dy_inv_out = 1.0_dp/dy_out
    ! In the vertical directions
    dz_out(1) = zp(1)
    do ii=2,nz_in
        ! dz 
        dz_out(ii) = zin(ii) - zin(ii-1)
        ! inverse of dz
        dz_inv_out(ii) = 1.0_dp/dz_out(ii) 
    end do

  end subroutine setup_grid_spacing

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

  subroutine tagminmax(xin,yin,zin,bbox_min,bbox_max,nx_in,ny_in,nz_in,dx_in,dy_in,dz_in,sx,ex,sy,ey,sz,ez)
    !
    ! This subroutine tags the location of the min and max indices based on bounding box
    !
    implicit none
    ! Input Argument
    real(dp), dimension(:), intent(in) :: xin, yin, zin                 ! x, y, z coordinates
    real(dp), intent(in) :: dx_in, dy_in, dz_in                            ! grid spacing
    real(dp), dimension(:), intent(in) :: bbox_min, bbox_max      ! Bounding box min & max
    integer, intent(in) :: nx_in, ny_in, nz_in                             ! Grid points
    ! Output
    integer, intent(out) :: sx,ex,sy,ey,sz,ez                     ! Start-end index
    ! Local variables
    integer :: iteri, iterj, iterk                                            ! Iterators
    real(dp) :: buffer_distance

    ! Set min and max to domain extents
    sx = 1
    ex = nx_in
    sy = 1
    ey = ny_in
    sz = 1
    ez = nz_in
    ! Compute the buffer distance based on used unput
    buffer_distance = real(sdfresolution,kind=dp)

    do iteri=1,nx_in
      ! Tag min
      if(xin(iteri)<=bbox_min(1)-buffer_distance*dx_in) then
          sx = iteri
      else if (xin(iteri)<=bbox_max(1)+buffer_distance*dx_in) then
          ex = iteri
      else
        exit 
      end if
    end do

    do iterj=1,ny_in
      ! Tag max
      if(yin(iterj)<=bbox_min(2)-buffer_distance*dy_in) then
          sy = iterj
      else if(yin(iterj)<=bbox_max(2)+buffer_distance*dy_in) then
          ey = iterj
      else
        exit 
      end if
    end do

    do iterk=1,nz_in
      ! Tag max
      if(zin(iterk)<=bbox_min(3)-buffer_distance*dz_in) then
          sz = iterk
      else if (zin(iterk)<=bbox_max(3)+buffer_distance*dz_in) then
          ez = iterk
      else
          exit
      end if
    end do

    ! Print to screen
    print *, "*** Min-Max Index-Value pair ***"
    print *, "Min-Max x:", sx, xin(sx), "|", ex, xin(ex)
    print *, "Min-Max y:", sy, yin(sy), "|", ey, yin(ey)
    print *, "Min-Max z:", sz, zin(sz), "|", ez, zin(ez)
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

  subroutine compute_scalar_distance_face(zstart,zend,xin,yin,zin,nfaces,input_faces,input_vertices,buffer_point_size,distance)
    !
    ! This subroutine computes the shortest distance between the input geometry face and the query points on the grid
    !
    implicit none
    ! Input
    integer, intent(in) :: zstart, zend
    real(dp), intent(in), dimension(:) :: xin, yin, zin
    integer, intent(in) :: nfaces
    real(dp), intent(in), dimension(:,:) :: input_vertices
    integer, intent(in), dimension(:,:) :: input_faces
    integer, intent(in) :: buffer_point_size
    ! Output
    real(dp), intent(out), dimension(:,:,:) :: distance
    ! Local variables
    integer :: face_id, ii, jj, kk 
    real(dp), dimension(3) :: vertex_1, vertex_2, vertex_3
    real(dp), dimension(3) :: min_query, max_query
    integer, dimension(3) :: min_index, max_index
    real(dp), dimension(3) :: query_point
    real(dp) :: temp_distance
    real(dp) :: deltax_inverse, deltay_inverse

    ! Compute the grid spacing and inverse
    deltax_inverse = 1.0_dp/(xin(10) - xin(9))
    deltay_inverse = 1.0_dp/(yin(10) - yin(9))

    ! Loop over all the faces of the geometry
    do face_id=1,nfaces
      call show_progress(face_id,nfaces,50)
      ! Get the vertices corresponding to this face
      vertex_1 = input_vertices(:,input_faces(1,face_id))
      vertex_2 = input_vertices(:,input_faces(2,face_id))
      vertex_3 = input_vertices(:,input_faces(3,face_id))
      ! Compute the bounding box based on the min-max of these three vertices to shrink the query points
      do ii=1,3
        min_query(ii) = min(vertex_1(ii),vertex_2(ii),vertex_3(ii))
        max_query(ii) = max(vertex_1(ii),vertex_2(ii),vertex_3(ii))
      end do
      ! Now locate the min and max indices to loop over (face-local bounding box)
      min_index(1) = floor(min_query(1) * deltax_inverse) - buffer_point_size
      min_index(2) = floor(min_query(2) * deltay_inverse) - buffer_point_size
      max_index(1) = floor(max_query(1) * deltax_inverse) + buffer_point_size
      max_index(2) = floor(max_query(2) * deltay_inverse) + buffer_point_size
      ! Min locator
      kk = zstart
      min_index(3) = zstart       ! Assign first query point as minimum if condition is not met 
      do while (zin(kk) <= min_query(3))
        min_index(3) = kk - buffer_point_size
        kk = kk + 1
      end do
      ! Max locator
      kk = zstart
      max_index(3) = zend       ! Assign last query point as minimum if condition is not met 
      do while (zin(kk) <= max_query(3))
        max_index(3) = kk + buffer_point_size
        kk = kk + 1
      end do
      ! Now loop over the face-local bounding box and compute the distance to the triangle
      if (min_index(1) >= 1 .and. max_index(1) <= nx .and. min_index(2) >= 1 .and. max_index(2) <= ny .and. min_index(3) >= 1 .and. max_index(3) <= nz) then
        do kk=min_index(3),max_index(3)
          do jj=min_index(2),max_index(2)
            do ii=min_index(1),max_index(1)
              query_point = [xin(ii),yin(jj),zin(kk)]
              call distance_point_to_triangle(query_point, vertex_1, vertex_2, vertex_3, temp_distance)
              distance(ii,jj,kk) = min(temp_distance,distance(ii,jj,kk))
            end do
          end do
        end do
      endif

    end do

  end subroutine compute_scalar_distance_face

  subroutine distance_point_to_triangle(point, vert0, vert1, vert2, dist_to_face)
    !
    ! This subroutine computes the distance between a point and a triangle
    !
    implicit none
    ! Input
    real(dp), intent(in) :: point(3)                        ! Query point coordinates
    real(dp), intent(in) :: vert0(3), vert1(3), vert2(3)    ! Triangle vertices
    ! Output
    real(dp), intent(out) :: dist_to_face         ! Distance between point and triangle
    ! Local variables
    real(dp) :: edge0(3), edge1(3), v0p(3), n(3), dist_vec(3), proj(3)
    real(dp) :: d, denom, param_s, param_t

    ! Compute the edges of the triangle
    edge0 = vert1 - vert0                               ! Edge v0->v1
    edge1 = vert2 - vert0                               ! Edge v0->v2
    v0p = point - vert0                                 ! Vector from v0 to point p

    ! Compute the normal of the triangle
    call cross_product(n, edge0, edge1)

    ! Projection of the point onto the triangle plane
    d = dot_product(n, vert0)
    denom = sqrt(dot_product(n, n))

    if (denom > 1.0e-12_dp) then
        dist_to_face = abs(dot_product(n, point) - d) / denom
    else
        dist_to_face = 0.0_dp
    end if

    ! Check if the point is inside the triangle (using barycentric coordinates)
    call project_point_to_triangle(point, vert0, vert1, vert2, param_s, param_t)

    ! Compute the closest point on the triangle using s and t
    proj = vert0 + param_s * edge0 + param_t * edge1

    ! Calculate the distance between the query point and the closest point on the triangle
    dist_vec = point - proj
    dist_to_face = sqrt(dot_product(dist_vec, dist_vec))

  end subroutine distance_point_to_triangle

  subroutine project_point_to_triangle(point, vert0, vert1, vert2, param_s, param_t)
    !
    ! This subroutine projects a point onto a triangle and calculates the barycentric coordinates (s, t).
    !
    implicit none
    ! Input
    real(dp), intent(in) :: point(3)                        ! Query point coordinates
    real(dp), intent(in) :: vert0(3), vert1(3), vert2(3)    ! Triangle vertices
    ! Output
    real(dp), intent(out) :: param_s, param_t               ! Barycentric coordinates of the point projection

    ! Local variables
    real(dp) :: edge0(3), edge1(3), v0p(3)
    real(dp) :: a, b, c, d, e, det

    ! Compute the edges of the triangle
    edge0 = vert1 - vert0   ! Edge v0->v1
    edge1 = vert2 - vert0   ! Edge v0->v2
    v0p = point - vert0     ! Vector from v0 to point p

    ! Compute dot products for the barycentric coordinates
    a = dot_product(edge0, edge0)
    b = dot_product(edge0, edge1)
    c = dot_product(edge1, edge1)
    d = dot_product(edge0, v0p)
    e = dot_product(edge1, v0p)

    ! Determinant of the matrix
    det = a * c - b * b

    if (det /= 0.0_dp) then
        ! Compute barycentric coordinates (s, t) of the projection
        param_s = (b * e - c * d) / det
        param_t = (b * d - a * e) / det

        ! Clamp s and t to ensure the point lies inside the triangle or on the edges
        if (param_s < 0.0_dp) then
            param_s = 0.0_dp
        else if (param_s > 1.0_dp) then
            param_s = 1.0_dp
        end if

        if (param_t < 0.0_dp) then
            param_t = 0.0_dp
        else if (param_t > 1.0_dp) then
            param_t = 1.0_dp
        end if
    else
        ! Handle degenerate triangles (when det is zero, the points are collinear)
        param_s = 0.0_dp
        param_t = 0.0_dp
    end if

  end subroutine project_point_to_triangle

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
  attributes(global) subroutine get_signed_distance_cuda(xd, yd, zd, dnumrays, dnum_faces, dnumberofnarrowpoints, dfaces, dvertices, ddirections, dnarrowbandindices, dpointinside)
    implicit none
    ! Input
    real(dp), intent(in), dimension(:), device :: xd, yd, zd
    integer, intent(in), device :: dnumberofnarrowpoints, dnum_faces, dnumrays
    integer, intent(in), dimension(:,:), device :: dnarrowbandindices
    integer, intent(in), dimension(:,:), device :: dfaces
    real(dp), intent(in), dimension(:,:), device :: dvertices
    real(dp), intent(in), dimension(:,:), device :: ddirections 
    ! Output
    real(dp), intent(inout), dimension(:,:,:), device :: dpointinside
    ! Local variables
    integer :: dpoint_tracker, dcid, dfaceid, di, dj, dk
    integer :: v1, v2, v3
    integer :: intersection_count
    real(dp), dimension(3) :: querypoint
    logical :: intersect
    real(dp) :: setsign
    real(dp) :: t, u, v

    ! Get the thread index for parallel execution
    dpoint_tracker = blockIdx%x * blockDim%x + threadIdx%x + 1

    ! Ensure query is within bounds for the narrowband indices
    if (dpoint_tracker <= dnumberofnarrowpoints) then
        ! Extract di, dj, dk from narrowband indices (add boundary checks)
        di = dnarrowbandindices(1, dpoint_tracker)
        dj = dnarrowbandindices(2, dpoint_tracker)
        dk = dnarrowbandindices(3, dpoint_tracker)

        ! Ensure di, dj, dk are within valid bounds of xd, yd, zd, and dpointinside
        if (di > 0 .and. di <= size(xd) .and. dj > 0 .and. dj <= size(yd) .and. dk > 0 .and. dk <= size(zd)) then
            querypoint = [xd(di), yd(dj), zd(dk)]
            intersection_count = 0

            ! Loop over rays and faces
            do dcid = 1, dnumrays
                do dfaceid = 1, dnum_faces
                    ! Extract vertices for the face
                    v1 = dfaces(1, dfaceid)
                    v2 = dfaces(2, dfaceid)
                    v3 = dfaces(3, dfaceid)

                    ! Perform ray-triangle intersection test
                    call ray_triangle_intersection_cuda(querypoint, ddirections(:, dcid), &
                                                        dvertices(:, v1), dvertices(:, v2), dvertices(:, v3), &
                                                        t, u, v, intersect)
                    if (intersect) then
                        intersection_count = intersection_count + 1
                    end if
                end do
            end do

            ! Set the sign based on the number of intersections
            if (mod(intersection_count, 2) > 0) then
                setsign = -1.0_dp
            else
                setsign = 1.0_dp
            end if

            ! Update dpointinside array
            if (di > 0 .and. di <= size(dpointinside, 1) .and. &
                dj > 0 .and. dj <= size(dpointinside, 2) .and. &
                dk > 0 .and. dk <= size(dpointinside, 3)) then
                dpointinside(di, dj, dk) = setsign * dpointinside(di, dj, dk)
            end if
        end if
    end if
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
