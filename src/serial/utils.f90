module utils
    implicit none
    ! Real data type selector
    integer, parameter, public :: dp = selected_real_kind(15,307), &
                                  sp = selected_real_kind(6 , 37)
    !                              
    ! Input data [Public data declared first then the rest]
    !    
    ! -- Geometry data --
    character(len=512) :: inputfilename
    ! -- Cartesian grid data --
    integer, protected :: nx, ny, nz
    real(dp), dimension(3), protected :: r0
    integer, dimension(3) :: ng
    logical, protected :: non_uniform_grid
    real(dp), allocatable, dimension(:) :: xp, yp, zp, xf, yf, zf
    real(dp) :: lx, ly, lz, dx, dy, dx_inverse, dy_inverse
    real(dp), allocatable, dimension(:) :: dz, dz_inverse
    ! -- SDF data --
    real(dp), protected :: scalarvalue
    integer, protected :: buffer_points
    ! -- Auxiliar data --
    integer, protected :: pbarwidth

contains

    subroutine read_inputfile()
        !
        ! This subroutine reads the input file
        !
        implicit none
        ! Input
        integer :: iunit, ierr              ! File I/O variables
        character(len=512) :: dummyline     ! dummyline to read help info in the inputfile
        
        iunit = 10
        open(newunit=iunit,file='parameters.in',status='old',action='read',iostat=ierr)
        if( ierr == 0 ) then
            read(iunit,*,iostat=ierr) dummyline
            read(iunit,*,iostat=ierr) inputfilename
            read(iunit,*,iostat=ierr) dummyline
            read(iunit,*,iostat=ierr) scalarvalue, buffer_points, pbarwidth
            read(iunit,*,iostat=ierr) dummyline
            read(iunit,*,iostat=ierr) nx, ny, nz
            read(iunit,*,iostat=ierr) dummyline
            read(iunit,*,iostat=ierr) r0(1), r0(2), r0(3)
            read(iunit,*,iostat=ierr) dummyline
            read(iunit,*,iostat=ierr) non_uniform_grid
            ! Setup what ng is 
            ng(1) = nx
            ng(2) = ny
            ng(3) = nz    
            ! Print to screen
            print *, "*** Input file sucessfully read ***"
        else
            error stop "ERROR: parameters.in file encountered a problem! :("
        end if
        
    end subroutine read_inputfile

    subroutine read_cans_grid(procid,loc, iprecision, npoints, origin, non_uni_grid, xin_p, yin_p, zin_p, xin_f, yin_f, zin_f, dzin, inputlx, inputly, inputlz)
        !
        ! This subroutine reads the CaNS grid to memory
        !
        ! Input
        integer, intent(in) :: procid                         ! Processor ID
        character(len=*), intent(in) :: loc                   ! Directory location of the grid file
        integer, intent(in) :: iprecision                     ! Precision level (4 (single) or 8 (double))
        real(dp), dimension(3), intent(in) :: origin          ! Location of the origin
        logical, intent(in) :: non_uni_grid                   ! Flag for non-uniform grid
        integer, dimension(3), intent(out) :: npoints         ! Grid size in x, y, z (modified by geometry file)
        ! Output
        real(dp), intent(out), allocatable, dimension(:) :: xin_p, yin_p, zin_p   ! Centered grid arrays
        real(dp), intent(out), allocatable, dimension(:) :: xin_f, yin_f, zin_f   ! Staggered grid arrays
        real(dp), intent(out), allocatable, dimension(:) :: dzin                  ! Grid spacing in z
        real(dp), intent(out) :: inputlx, inputly, inputlz                        ! Domain size in x y and z
        ! Local variables
        logical :: fexists                                    ! File exists boolean
        integer :: iter, ii, jj                               ! Local Iterator
        real(dp) :: dl(3), l(3)                               ! Grid spacing and length
        real(sp), allocatable :: grid_z4(:,:)                 ! For non-uniform grid, single precision
        real(dp), allocatable :: grid_z8(:,:)                 ! For non-uniform grid, double precision
        character(len=512) :: geofile, grdfile                ! Filenames for geometry.out and grid.out
        integer :: filerr                                     ! File error handle
    
        ! Check the file directory
        inquire(file=trim(loc), exist=fexists)
        if (fexists .eqv. .false.) then
            print*, "The input directory does not exist: ", trim(loc)
            stop
        endif
    
        ! File paths
        geofile = trim(loc) // "geometry.out"
    
        ! Read geometry file
        open(unit=10, file=geofile, status='old')
        read(10, *) npoints      ! Read the grid size (x, y, z)
        read(10, *) l            ! Read the domain lengths
        close(10)
    
        ! Set the values for inputlx, inputly, inputlz
        inputlx = l(1)
        inputly = l(2)
        inputlz = l(3)
    
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
            grdfile = trim(loc) // "grid.out"
            
            if (iprecision == 4) then
                open(unit=20, file=grdfile, status = "old", action="read", iostat=filerr)
                allocate(grid_z4(npoints(3),5))
                if (filerr /= 0) then
                    error stop "Error opening grid file...."
                endif

                do ii=1,npoints(3)
                    read(20,*,iostat=filerr) (grid_z4(ii,jj), jj=1,5)
                    if (filerr /= 0) then
                        error stop "Error reading line in file...."
                    endif
                end do
                close(20)
                ! Allocate the grid to the output arrays
                zin_p = origin(3) + grid_z4(:,2)
                zin_f = origin(3) + grid_z4(:,3)
            else
                open(unit=20, file=grdfile, status = "old", action="read", iostat=filerr)
                allocate(grid_z8(npoints(3),5))
                if (filerr /= 0) then
                    error stop "Error opening grid file...."
                endif

                do ii=1,npoints(3)
                    read(20,*,iostat=filerr) (grid_z8(ii,jj), jj=1,5)
                    if (filerr /= 0) then
                        error stop "Error reading line in file...."
                    endif
                end do
                close(20)
                ! Allocate the grid to the output arrays
                zin_p = origin(3) + grid_z8(:,2)
                zin_f = origin(3) + grid_z8(:,3)
            endif                            
        endif
    
        if (procid == 0) print *, "*** Successfully read the CaNS grid ***"
    
    end subroutine read_cans_grid

    subroutine read_obj(filename, vertices, normals, faces, face_normals, num_vertices, num_normals, num_faces)
        !
        ! This subroutine reads an OBJ file and extracts vertices, normals, and faces.
        !
        implicit none
        ! Inputs
        character(len=*), intent(in) :: filename                ! Name of the OBJ file
        ! Outputs
        real(dp), allocatable, intent(out) :: vertices(:,:)   ! 2D array to store vertex coordinates (3 x num_vertices)
        real(dp), allocatable, intent(out) :: normals(:,:)    ! 2D array to store normal vectors (3 x num_normals)
        integer, allocatable, intent(out) :: faces(:,:), face_normals(:,:)       ! 2D array to store face vertex indices (3 x num_faces)    
        integer, intent(out) :: num_vertices, num_normals, num_faces  ! Number of vertices, normals, and faces    
        ! Local variables
        integer :: unit, ios, v_count, vn_count, f_count
        character(len=512) :: line                               ! Line buffer to read the file
        real(dp) :: verx, very, verz                             ! Variables for vertex coordinates
        real(dp) :: normx, normy, normz                          ! Variables for normal vector components
        integer :: verind1, verind2, verind3                     ! Variables for vertex indices in a face
        integer :: vn1, vn2, vn3                                 ! Variables for normal indices in a face
        character(len=512) :: dline
    
        ! Open the OBJ file
        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'Error opening file: ', filename
            stop
        end if
    
        ! Initialize counters
        v_count = 0
        vn_count = 0
        f_count = 0
    
        ! First pass: Count vertices, normals, and faces
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
    
            select case (trim(adjustl(line(1:2))))
                case ('v ')
                    v_count = v_count + 1
                case ('vn')
                    vn_count = vn_count + 1
                case ('f ')
                    f_count = f_count + 1
            end select
        end do
    
        if(vn_count == 0) then
            print *, "FATAL ERROR: Geometry does not have normals information!"
            print *, "Please use trimesh library based - python script to mesh.fix_normals() and output the geometry with normals to fix this error."
            error stop
        endif

        ! Allocate arrays
        allocate(vertices(3, v_count))
        allocate(normals(3, vn_count))
        allocate(faces(3, f_count))
        allocate(face_normals(3, f_count))
    
        ! Store the counts
        num_vertices = v_count
        num_normals = vn_count
        num_faces = f_count
    
        ! Rewind the file to read it again
        rewind(unit)
    
        ! Reset counters
        v_count = 0
        vn_count = 0
        f_count = 0
    
        ! Second pass: Read vertices, normals, and faces
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
    
            ! Read vertex data
            if (trim(adjustl(line(1:2))) == 'v ') then
                read(line(3:), *) verx, very, verz
                v_count = v_count + 1
                vertices(1, v_count) = verx
                vertices(2, v_count) = very
                vertices(3, v_count) = verz
    
            ! Read normal data
            else if (trim(adjustl(line(1:3))) == 'vn ') then
                read(line(4:), *) normx, normy, normz
                vn_count = vn_count + 1
                normals(1, vn_count) = normx
                normals(2, vn_count) = normy
                normals(3, vn_count) = normz
    
            ! Read face data (assuming the format f v1//vn1 v2//vn2 v3//vn3)
            else if (trim(adjustl(line(1:2))) == 'f ') then
                dline = trim(line(3:))
                call parse_face_line(dline, verind1, vn1, verind2, vn2, verind3, vn3)
                f_count = f_count + 1
                faces(1, f_count) = verind1
                faces(2, f_count) = verind2
                faces(3, f_count) = verind3
                ! Face normal logic (No assumption that vertex id and vertex normal id is the same order)
                face_normals(1, f_count) = vn1
                face_normals(2, f_count) = vn2
                face_normals(3, f_count) = vn3
            end if
        end do
    
        ! Close the file
        close(unit)
            
        print *, 'Successfully read OBJ file: ', trim(filename)
        print *, 'Number of vertices: ', num_vertices
        print *, 'Number of normals: ', num_normals
        print *, 'Number of faces: ', num_faces
    
    end subroutine read_obj

    subroutine parse_face_line(line, v1, vn1, v2, vn2, v3, vn3)
        !
        ! Parses a face line in the format f v1//vn1 v2//vn2 v3//vn3
        !
        implicit none
        ! Input
        character(len=512), intent(inout) :: line               ! Line to be parsed
        ! Output
        integer, intent(out) :: v1, vn1, v2, vn2, v3, vn3       ! Vertex and Vertex normal values
        ! Local Variables
        character(len=32) :: part                               ! Temp char variable to split line
        integer :: pos                                          ! Position counter
    
        ! Parse first vertex-normal pair (v1//vn1)
        pos = index(line, ' ')
        part = trim(line(1:pos-1))
        call parse_vertex_normal_pair(part, v1, vn1)
        line = trim(line(pos+1:))
    
        ! Parse second vertex-normal pair (v2//vn2)
        pos = index(line, ' ')
        part = trim(line(1:pos-1))
        call parse_vertex_normal_pair(part, v2, vn2)
        line = trim(line(pos+1:))
    
        ! Parse third vertex-normal pair (v3//vn3)
        part = trim(line)
        call parse_vertex_normal_pair(part, v3, vn3)
    
    end subroutine parse_face_line
    
    subroutine parse_vertex_normal_pair(pair, vertex, normal)
        !
        ! Parses a vertex/normal pair in the format v//vn
        !
        implicit none
        ! Input
        character(len=*), intent(inout) :: pair                     ! Character to parse pair
        ! Output
        integer, intent(out) :: vertex, normal                      ! Vertex and normal ID
        ! Local Variables
        integer :: double_slash_pos                                 ! Identifier for double slash to trim
    
        ! Find the position of the double slash (//)
        double_slash_pos = index(pair, '//')
        if (double_slash_pos > 0) then
            read(pair(1:double_slash_pos-1), *) vertex
            read(pair(double_slash_pos+2:), *) normal
        else
            ! If there is no normal index, handle as an error or assign default
            print *, "Warning: No normal index found in pair ", trim(pair)
            vertex = 0
            normal = 0
        end if
    
    end subroutine parse_vertex_normal_pair

    subroutine setup_grid_spacing(xin,yin,zin,nz_in,dx_out,dy_out,dz_out,dx_inv_out,dy_inv_out,dz_inv_out)
        !
        ! This subroutine setups up the grid spacing and its inverse
        !
        implicit none
        ! Input
        real(dp), intent(in), dimension(:) :: xin, yin, zin                     ! Input grid arrays
        integer, intent(in) :: nz_in                                            ! Number of points in z
        ! Output
        real(dp), intent(out) :: dx_out, dy_out, dx_inv_out, dy_inv_out         ! Grid size in x and y (normal + inverse)
        real(dp), intent(out), dimension(:) :: dz_out, dz_inv_out               ! Grid size in z (arrays, normal + inverse)
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
    
        print * , "*** Sucessfully finished setting up the grid spacing ***"
    end subroutine setup_grid_spacing

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

        ! Print the information to screen
        print *, "Geometry is bounded by (minimum)", bbox_min
        print *, "Geometry is bounded by (maximum)", bbox_max
    end subroutine getbbox

    subroutine tagminmax(xin,yin,zin,bbox_min,bbox_max,nx_in,ny_in,nz_in,dx_in,dy_in,dz_in,sdfresolution,sx,ex,sy,ey,sz,ez)
        !
        ! This subroutine tags the location of the min and max indices based on bounding box
        !
        implicit none
        ! Input Argument
        real(dp), dimension(:), intent(in) :: xin, yin, zin                     ! x, y, z coordinates
        real(dp), intent(in) :: dx_in, dy_in, dz_in                             ! grid spacing
        real(dp), dimension(:), intent(in) :: bbox_min, bbox_max                ! Bounding box min & max
        integer, intent(in) :: nx_in, ny_in, nz_in                              ! Grid points
        integer, intent(in) :: sdfresolution                                    ! Number of buffer grid points around the bounding box
        ! Output
        integer, intent(out) :: sx,ex,sy,ey,sz,ez                               ! Start-end index
        ! Local variables
        integer :: iteri, iterj, iterk                                          ! Iterators
        real(dp) :: buffer_distance                                             ! Buffer distance check
    
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

    subroutine compute_scalar_distance_face(xstart,xend,ystart,yend,zstart,zend,xin,yin,zin,nfaces,input_faces,input_face_normals,input_vertices,input_normals,buffer_point_size,distance)
        !
        ! This subroutine computes the shortest distance between the input geometry face and the query points on the grid
        !
        implicit none
        ! Input
        integer, intent(in) :: xstart,xend,ystart,yend, zstart, zend            ! Start and Ending indices for the bounding box loop
        real(dp), intent(in), dimension(:) :: xin, yin, zin                     ! Cartesian grid over which SDF is to be computed
        integer, intent(in) :: nfaces                                           ! Number of faces in the geometry
        real(dp), intent(in), dimension(:,:) :: input_vertices, input_normals   ! Vertex and vertex normal arrays
        integer, intent(in), dimension(:,:) :: input_faces, input_face_normals  ! Faces and face normal indices
        integer, intent(in) :: buffer_point_size                                ! Number of buffer points - Narrow band size
        ! Output
        real(dp), intent(out), dimension(:,:,:) :: distance                     ! Signed-Distance-Array
        ! Local variables
        integer :: face_id, ii, jj, kk                                          ! Grid and face iterators
        real(dp), dimension(3) :: vertex_1, vertex_2, vertex_3                  ! Vertices for nth face
        real(dp), dimension(3) :: norm_vert_1, norm_vert_2, norm_vert_3         ! Vertex normals for nth face
        real(dp), dimension(3) :: min_query, max_query                          ! Min-max query points
        integer, dimension(3) :: min_index, max_index                           ! Min-max query point indices
        real(dp), dimension(3) :: query_point                                   ! Local query point
        real(dp) :: temp_distance                                               ! nth face distance - temp comparison
        real(dp) :: avg_normal_mag_inv                                          ! Inverse of average normal L2 norm
        real(dp) :: deltax_inverse, deltay_inverse                              ! Grid size in x and y - Inverse 
        real(dp) :: avg_normal(3)                                               ! Average face normal (only for orientation)
        real(dp) :: temp_value(2)                                               ! Temporary signed value array for comparison (signed)                                            
        real(dp) :: stime                                                       ! Starting time of the calculation
            
        ! Query the start time to add to the progress bar
        call cpu_time(stime)        

        ! Compute the grid spacing and inverse
        deltax_inverse = 1.0_dp/(xin(10) - xin(9))
        deltay_inverse = 1.0_dp/(yin(10) - yin(9))
    
        ! Loop over all the faces of the geometry
        do face_id=1,nfaces
            call show_progress(face_id,nfaces,pbarwidth,stime)
            ! Get the vertices corresponding to this face
            vertex_1 = input_vertices(:,input_faces(1,face_id))
            vertex_2 = input_vertices(:,input_faces(2,face_id))
            vertex_3 = input_vertices(:,input_faces(3,face_id))
            ! Get the vertex normals corresponding to the this faces
            norm_vert_1 = input_normals(:,input_face_normals(1,face_id))
            norm_vert_2 = input_normals(:,input_face_normals(2,face_id))
            norm_vert_3 = input_normals(:,input_face_normals(3,face_id))
            ! Get the average normal cooresponding to this face
            avg_normal(1) = norm_vert_1(1) + norm_vert_2(1) + norm_vert_3(1)
            avg_normal(2) = norm_vert_1(2) + norm_vert_2(2) + norm_vert_3(2)
            avg_normal(3) = norm_vert_1(3) + norm_vert_2(3) + norm_vert_3(3)
            ! Get magnitude of the avg_normal
            avg_normal_mag_inv = 1.0_dp/sqrt(dot_product(avg_normal,avg_normal))
            ! Normalise the vector
            avg_normal = avg_normal*avg_normal_mag_inv
            ! Compute the bounding box based on the min-max of these three vertices to shrink the query points
            do ii=1,3
                min_query(ii) = min(vertex_1(ii),vertex_2(ii),vertex_3(ii)) 
                max_query(ii) = max(vertex_1(ii),vertex_2(ii),vertex_3(ii))       
            end do
            ! Now locate the min indices to loop over (face-local bounding box)
            min_index(1) = floor(min_query(1) * deltax_inverse) - buffer_point_size
            max_index(1) = floor(max_query(1) * deltax_inverse) + buffer_point_size
            min_index(2) = floor(min_query(2) * deltay_inverse) - buffer_point_size
            max_index(2) = floor(max_query(2) * deltay_inverse) + buffer_point_size                        
            ! Z direction can be non-homogeneous, hence special case loop
            kk = zstart
            min_index(3) = zend  ! Initialize to the maximum possible index
            do while (kk <= zend)
                if (zin(kk) > min_query(3)) exit
                min_index(3) = kk
                kk = kk + 1
            end do
            min_index(3) = max(min_index(3) - buffer_point_size, zstart)

            ! Max query for Z direction
            kk = zstart
            max_index(3) = zstart  ! Initialize to the minimum possible index
            do while (kk <= zend)
                if (zin(kk) > max_query(3)) exit
                max_index(3) = kk
                kk = kk + 1
            end do
            max_index(3) = min(max_index(3) + buffer_point_size, zend)
            
            ! Now loop and calculate distance 
            do kk=min_index(3),max_index(3)
                do jj=min_index(2),max_index(2)
                    do ii=min_index(1),max_index(1)                        
                        ! Check if indices are within bounds
                        if (ii >= xstart .and. ii <= xend .and. jj >= ystart .and. jj <= yend .and. kk >= zstart .and. kk <= zend) then
                            query_point = [xin(ii),yin(jj),zin(kk)]
                            call distance_point_to_triangle(query_point, vertex_1, vertex_2, vertex_3, avg_normal, temp_distance)
                            temp_value = [temp_distance,distance(ii,jj,kk)]
                            distance(ii,jj,kk) = temp_value(minloc(abs(temp_value),dim=1))
                        end if
                    end do
                end do
            end do                  
    
        end do
    
    end subroutine compute_scalar_distance_face

    subroutine fill_internal(grid, inx, iny, inz, sx, sy, sz, ex, ey, ez, large_negative_value)
        ! 
        ! This subroutine uses the flood-fill algorithm to fill in values inside the geometry
        !        
        implicit none
        ! Input
        integer, intent(in) :: inx, iny, inz                ! Number of grid points
        integer, intent(in) :: sx, sy, sz, ex, ey, ez       ! Starting and ending indices of the bounding box
        real(dp), intent(in) :: large_negative_value        ! Large negative value to be set inside
        ! Output
        real(dp), intent(inout) :: grid(inx, iny, inz)      ! Array with distance pre-computed            
        ! Local variables
        logical, allocatable :: visited(:,:,:)              ! Boolean array for points visited
        integer :: ii, jj, kk, ci, cj, ck                   ! Iterators
        integer :: top, i_dir                               ! Stack locators for visited points
        integer, allocatable :: stack(:,:)                  ! Array storing the visited points stack
        integer :: ni, nj, nk                               ! Local iterators            
        integer, parameter :: n_dirs = 6                    ! Directions for neighbors
        integer :: dirs(3, n_dirs) = reshape([ &            ! Array with directions i.e., 6
            1,  0,  0, -1,  0,  0, &    ! x+1, x-1
            0,  1,  0,  0, -1,  0, &    ! y+1, y-1
            0,  0,  1,  0,  0, -1], &   ! z+1, z-1
            shape=[3, n_dirs])
    
        ! Allocate arrays
        allocate(visited(inx, iny, inz))
        allocate(stack(3, (ex - sx + 1) * (ey - sy + 1) * (ez - sz + 1)))  ! Stack for flood fill, size large enough for the entire subregion
    
        ! Initialize the visited array to false
        visited = .false.
    
        ! Mark all points on the boundaries of the subregion as starting points for exterior flood fill
        top = 0
        do kk = sz, ez
            do jj = sy, ey
                do ii = sx, ex
                    if (ii == sx .or. ii == ex .or. jj == sy .or. jj == ey .or. kk == sz .or. kk == ez) then
                        if (grid(ii, jj, kk) >= 0.0_dp .and. .not. visited(ii, jj, kk)) then
                            ! Add boundary point to stack for flood fill
                            top = top + 1
                            stack(:, top) = [ii, jj, kk]
                            visited(ii, jj, kk) = .true.
                        end if
                    end if
                end do
            end do
        end do
    
        ! Flood fill to mark all exterior points that are reachable from the boundary
        do while (top > 0)
            ! Pop the current point from the stack
            ci = stack(1, top)
            cj = stack(2, top)
            ck = stack(3, top)
            top = top - 1
    
            ! Traverse all neighboring points
            do i_dir = 1, n_dirs
                ni = ci + dirs(1, i_dir)
                nj = cj + dirs(2, i_dir)
                nk = ck + dirs(3, i_dir)
    
                ! Check if the neighbor is within bounds of the subregion
                if (ni >= sx .and. ni <= ex .and. &
                    nj >= sy .and. nj <= ey .and. &
                    nk >= sz .and. nk <= ez) then
    
                    if (.not. visited(ni, nj, nk) .and. grid(ni, nj, nk) >= 0.0_dp) then
                        ! Mark as visited and add to stack
                        top = top + 1
                        stack(:, top) = [ni, nj, nk]
                        visited(ni, nj, nk) = .true.
                    end if
                end if
            end do
        end do
    
        ! Set a large value for points are never visited then require filling
        do kk = sz, ez
            do jj = sy, ey
                do ii = sx, ex                    
                    if (.not. visited(ii, jj, kk)) then
                        grid(ii, jj, kk) = large_negative_value
                    end if
                end do
            end do
        end do
    
        ! Deallocate arrays
        deallocate(visited)
        deallocate(stack)
    
    end subroutine fill_internal

    subroutine distance_point_to_triangle(point, vert0, vert1, vert2, avg_normal, dist_to_face)
        !
        ! This subroutine computes the distance between a point and a triangle
        !
        implicit none
        ! Input
        real(dp), intent(in) :: point(3)                        ! Query point coordinates
        real(dp), intent(in) :: vert0(3), vert1(3), vert2(3)    ! Triangle vertices
        real(dp), intent(in) :: avg_normal(3)                   ! Avg_normal (From OBJ data)
        ! Output
        real(dp), intent(out) :: dist_to_face                   ! Distance between point and triangle
        ! Local variables
        real(dp) :: edge0(3), edge1(3), v0p(3), n(3), dist_vec(3), proj(3)  ! Edge vectors 
        real(dp) :: d, denom, param_s, param_t                  ! Plane parameters on the triangle
        real(dp) :: sign_indicator, temp_dot_result             ! Temporary results 

        ! Compute the edges of the triangle
        edge0 = vert1 - vert0                               ! Edge v0->v1
        edge1 = vert2 - vert0                               ! Edge v0->v2
        v0p = point - vert0                                 ! Vector from v0 to point p    
        ! Compute the normal of the triangle
        call cross_product(n, edge0, edge1)        
        ! Projection of the point onto the triangle plane
        d = dot_product(n, vert0)        
        denom = sqrt(dot_product(n, n))        
        ! Note that n here may not be outward pointing since the face orientation is not known
        n = n/denom     ! Normalise the normal
        ! Check if the point is inside the triangle (using barycentric coordinates)
        call project_point_to_triangle(point, vert0, vert1, vert2, param_s, param_t)    
        ! Compute the closest point on the triangle using s and t
        proj = vert0 + param_s * edge0 + param_t * edge1    
        ! Calculate the distance between the query point and the closest point on the triangle
        dist_vec = point - proj
        temp_dot_result = dot_product(dist_vec,avg_normal)
        ! Perpendicularity check for query point
        sign_indicator = sign(1.0_dp,temp_dot_result)
        dist_to_face = sign_indicator*sqrt(dot_product(dist_vec, dist_vec))
        
    end subroutine distance_point_to_triangle

    subroutine project_point_to_triangle(point, vert0, vert1, vert2, param_s, param_t)
        !
        ! This subroutine projects a point onto a triangle and calculates the barycentric coordinates (s, t).
        ! Algorithm based on: https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
        ! and some ideas/understanding from : https://3d.bk.tudelft.nl/courses/backup/geo1004/2023/
        implicit none
        ! Input
        real(dp), intent(in) :: point(3)                        ! Query point coordinates
        real(dp), intent(in) :: vert0(3), vert1(3), vert2(3)    ! Triangle vertices
        ! Output
        real(dp), intent(out) :: param_s, param_t               ! Barycentric coordinates of the point projection    
        ! Local variables
        real(dp) :: edge0(3), edge1(3), v0p(3)                  ! Edge vectors
        real(dp) :: a, b, c, d, e, det                          ! Parameters for dot-product

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
        ! Check the location of the point in the s-t plane    
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

    subroutine write2binary(outputfilename,finaloutput)
        !
        ! This subroutine writes the output to binary file
        !
        implicit none
        ! Input
        character(len=*), intent(in) :: outputfilename              ! Outputfilename
        real(dp), intent(in), dimension(:,:,:) :: finaloutput       ! Array to be written

        ! Access=stream does not add 4 bytes at the start and end of the binary file
        open(unit=10, file=trim(outputfilename), access='stream', status='replace', form='unformatted')
        ! WARNING: Line below will add 4 bytes at the start and end of the file
        !open(unit=10, file=trim(outputfilename), status='replace', form='unformatted')
        write(10) finaloutput
        close(10)
    end subroutine write2binary

    subroutine show_progress(current, total, width, start_time)
        ! 
        ! This subroutine aims to display a progress bar and overwrite the same line.
        ! It also displays elapsed time and estimated time remaining.
        !
        implicit none
        ! Input
        integer, intent(in) :: current, total               ! Current and total values
        integer, intent(in) :: width                        ! Width of the progress bar
        real(dp), intent(in) :: start_time                  ! Start time of the process (computes estimated time)
        ! Local variables
        integer :: progress                                 ! Current progress counter
        real(dp) :: percent, elapsed_time, estimated_time   ! Percentage completition and time counters
        character(len=width) :: bar                         ! Character bar size
        real(dp) :: current_time                            ! Log for current time
    
        ! Get the current time
        call cpu_time(current_time)
    
        ! Calculate elapsed time
        elapsed_time = current_time - start_time
    
        ! Estimate remaining time
        if (current > 0) then
            estimated_time = elapsed_time / real(current) * (real(total) - real(current))
        else
            estimated_time = 0.0_dp
        end if
    
        ! Calculate percentage and progress bar length
        percent = real(current) / real(total) * 100.0
        progress = int(real(current) / real(total) * width)
    
        ! Construct the progress bar string
        bar = repeat('|', progress) // repeat(' ', width - progress)
    
        ! Print the progress bar with estimated time and elapsed time
        write(*,'(a, "|", a, "|", f6.2, "%", 1x, i8, "/", i8, " Elapsed: ", f8.2, "s Remaining: ", f8.2, "s")', advance='no') &
            char(13), bar, percent, current, total, elapsed_time, estimated_time
    
        ! If the current iteration is the last, move to a new line
        if (current == total) then
            print *  ! Move to the next line after the final iteration
        end if
    end subroutine show_progress

    subroutine printlogo()
        !
        ! This subroutine prints the logo of `genSDF`
        !
        implicit none
        print *, " ░▒▓██████▓▒░░▒▓████████▓▒░▒▓███████▓▒░ ░▒▓███████▓▒░▒▓███████▓▒░░▒▓████████▓▒░ "
        print *, "░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        "
        print *, "░▒▓█▓▒░      ░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        "
        print *, "░▒▓█▓▒▒▓███▓▒░▒▓██████▓▒░ ░▒▓█▓▒░░▒▓█▓▒░░▒▓██████▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓██████▓▒░   "
        print *, "░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        "
        print *, "░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        "
        print *, " ░▒▓██████▓▒░░▒▓████████▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓███████▓▒░░▒▓███████▓▒░░▒▓█▓▒░        "
                                                                                                                                                            
    end subroutine printlogo
end module utils
