module utils

    integer, parameter, public :: dp = selected_real_kind(15,307)
    character(len=512),  public :: inputfilename
    integer, public :: nx, ny, nz, numberofrays
    real(dp), public :: lx, ly, lz

contains

  subroutine read_inputfile()
    !
    ! This subroutine reads the input file
    !
    implicit none
    ! Local variables
    integer :: iunit, ierr
    character(len=512) :: dummyline
    
    iunit = 10
    open(newunit=iunit,file='parameters.in',status='old',action='read',iostat=ierr)
      if( ierr == 0 ) then
        read(iunit,*,iostat=ierr) dummyline
        read(iunit,*,iostat=ierr) inputfilename
        read(iunit,*,iostat=ierr) dummyline
        read(iunit,*,iostat=ierr) nx, ny, nz
        read(iunit,*,iostat=ierr) dummyline
        read(iunit,*,iostat=ierr) lx, ly, lz
        read(iunit,*,iostat=ierr) dummyline
        read(iunit,*,iostat=ierr) numberofrays
      else
        error stop "parameters.in file encountered a problem!"
      end if
  end subroutine read_inputfile

  subroutine read_obj(filename, vertices, faces, num_vertices, num_faces)
    !
    ! This subroutine reads an OBJ file and returns the vertices, faces, number of vertices, and number of faces
    !
    implicit none
    ! Inputs
    character(len=*), intent(in) :: filename  ! Name of the OBJ file
    ! Outputs
    real(dp), allocatable, intent(out) :: vertices(:,:)     ! 2D array to store vertex coordinates (3 x num_vertices)
    integer, allocatable, intent(out) :: faces(:,:)         ! 2D array to store faces (3 x num_faces)
    integer, intent(out) :: num_vertices                    ! Number of vertices
    integer, intent(out) :: num_faces                       ! Number of faces
    ! Local variables
    integer :: unit, ios, v_count, f_count
    character(len=512) :: line                              ! Line buffer to read the file
    real(dp) :: vx, vy, vz                                  ! Variables for vertex coordinates
    integer :: v1, v2, v3                                   ! Variables for face indices
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
    allocate(vertices(3, v_count))
    allocate(faces(3, f_count))

    ! Store the number of vertices and faces
    num_vertices = v_count
    num_faces = f_count

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
        vertices(1, v_count) = vx
        vertices(2, v_count) = vy
        vertices(3, v_count) = vz
      end if

      ! Read face data (only triangular faces assumed)
      if (line(1:2) == 'f ') then
        read(line(3:), *) v1, v2, v3
        f_count = f_count + 1
        faces(1, f_count) = v1
        faces(2, f_count) = v2
        faces(3, f_count) = v3
      end if
    end do

    ! Close the file
    close(unit)

    ! Log CPU time at the end of the operation
    call cpu_time(time2)
    print *, "*** File: '",trim(filename), "' sucessfully read in ",time2-time1, "seconds ***"
  end subroutine read_obj

  subroutine getbbox(vertices, num_vertices, bbox_min, bbox_max)
    !
    ! This subroutine calculates the bounding box for the OBJ geometry
    !
    implicit none
    ! Inputs
    real(dp), intent(in) :: vertices(3, num_vertices)  ! 2D array to store vertex coordinates (3 x num_vertices)
    integer, intent(in) :: num_vertices               ! Number of vertices  
    ! Outputs
    real(dp), intent(out) :: bbox_min(3)               ! Min values for x, y, z
    real(dp), intent(out) :: bbox_max(3)               ! Max values for x, y, z
    ! Local variables
    integer :: i
    real(dp) :: vx, vy, vz
  
    ! Initialize bbox_min and bbox_max to the first vertex
    bbox_min = vertices(:, 1)
    bbox_max = vertices(:, 1)
  
    ! Iterate over all vertices to find the bounding box
    do i = 2, num_vertices
      vx = vertices(1, i)
      vy = vertices(2, i)
      vz = vertices(3, i)
  
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

  subroutine getfacenormals(vertices, faces, num_faces, num_vertices,normals)
    !
    ! This subroutine computes the face normals for each face
    !
    implicit none  
    ! Inputs
    integer, intent(in) :: num_faces, num_vertices      ! Number of faces and vertices
    real(dp), intent(in) :: vertices(3,num_vertices)    ! Vertex coordinates (3 x num_vertices)
    integer, intent(in) :: faces(3,num_faces)           ! Face indices (3 x num_faces)  
    ! Outputs
    real(dp), allocatable, intent(out) :: normals(:,:)  ! Face normals (3 x num_faces)  
    ! Local variables
    integer :: i, v1, v2, v3
    real(dp) :: e1(3), e2(3), n(3), norm
    real(dp) :: time1, time2
  
    ! Log the CPU time at the start
    call cpu_time(time1)
    
    ! Allocate memory for normals
    allocate(normals(3, num_faces))
  
    ! Loop over each face to calculate the normal
    do i = 1, num_faces
      ! Get the vertex indices for the face (OBJ indices are 1-based)
      v1 = faces(1, i)
      v2 = faces(2, i)
      v3 = faces(3, i)
  
      ! Calculate edge vectors e1 and e2
      e1(1) = vertices(1, v2) - vertices(1, v1)
      e1(2) = vertices(2, v2) - vertices(2, v1)
      e1(3) = vertices(3, v2) - vertices(3, v1)
  
      e2(1) = vertices(1, v3) - vertices(1, v1)
      e2(2) = vertices(2, v3) - vertices(2, v1)
      e2(3) = vertices(3, v3) - vertices(3, v1)
  
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
      normals(1, i) = n(1)
      normals(2, i) = n(2)
      normals(3, i) = n(3)
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
    real(dp), dimension(:), intent(in) :: x, y, z
    real(dp), intent(in) :: dx, dy, dz
    real(dp), dimension(:), intent(in) :: bbox_min, bbox_max
    integer, intent(in) :: Nx, Ny, Nz
    ! Output
    integer, intent(out) :: sx,ex,sy,ey,sz,ez
    ! Local variables
    integer :: i, j, k

    do i=1,Nx
      ! Tag min
      if(x(i)<=bbox_min(1)-5*dx) then
        sx = i
      end if
      ! Tag max
      if(x(i)<=bbox_max(1)+5*dx) then
          ex = i
      end if
    end do

    do j=1,Ny
      ! Tag max
      if(y(j)<=bbox_min(2)-5*dy) then
          sy = j
      end if
      ! Tag min
      if(y(j)<=bbox_max(2)+5*dy) then
          ey = j
      end if
    end do

    do k=1,Nz
      ! Tag max
      if(z(k)<=bbox_min(3)-5*dz) then
          sz = k
      end if
      ! Tag min
      if(z(k)<=bbox_max(3)+5*dz) then
          ez = k
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
    real(dp), dimension(3), intent(in) :: u, v
    ! Output
    real(dp), dimension(3), intent(out) :: result

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
    epsilon = 1.0e-12
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
    write(*,'(a, "|", a, "|", f6.2, "%", 1x, i4, "/", i4)', advance='no') char(13), bar, percent, current, total

    ! If the current iteration is the last, move to a new line
    if (current == total) then
        print *  ! Move to the next line after the final iteration
    end if
  end subroutine show_progress
  

end module utils
