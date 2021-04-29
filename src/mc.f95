module functions
    use numeric_kinds
    use atoms_module,            only: atoms, finalise, &
                                     set_cutoff, &
                                     calc_connect, print, &
                                     map_into_cell, has_property
    use dictionary_module,       only: dictionary, initialise, finalise, &
                                     set_value, STRING_LENGTH
    use cinoutput_module,        only : read, write


real(dp), parameter :: pi = 3.141592653589793

contains


pure function distance(a, b, lattice) result(r)
    real(dp), intent(in) :: a(3), b(3), lattice(3)
    real(dp)             :: dr(3), r ! separation

    dr = abs(a - b)
    dr = min(dr, abs(lattice - dr))
    r = norm2(dr)
end function

subroutine rho_calc(pos, lattice, disc, ma, bw, bmp, rho)
    implicit none
    real(dp), intent(in), pointer    :: pos(:,:) 
    real(dp), intent(in)    :: lattice(3) ! 
    integer, intent(in)     :: disc ! real space discretisation (bins)
    integer                 :: i, j, loc
    real(dp), intent(out)                 :: ma, bw ! maximum r values to consider
    real(dp), intent(out), allocatable    :: bmp(:), rho(:) ! bin mid-points
    real(dp)                              :: dr


    ma = minval(lattice)/2
    bw = ma/disc

    allocate( bmp(disc), rho(disc ) )

    do i = 1, disc
        bmp(i) = bw * (i + 0.5)
    end do

    rho = 0.0

    ! main section of routine
    do i = 1, ubound(pos, dim=2)
        do j = i+1, ubound(pos, dim=2)
            dr = distance( pos(:,i), pos(:,j), lattice)
            loc = floor(dr / bw) + 1

            if (loc < disc) then
               rho(loc) = rho(loc) + 1
            end if
        end do
    end do 
end subroutine

pure function gr_from_rho(rho, rs, dr, rho_av, N, disc) result(gr)
  real(dp), intent(in), allocatable :: rho(:), rs(:)
  real(dp), intent(in)              :: dr, rho_av
  integer, intent(in)               :: N, disc
  real(dp), allocatable :: gr(:)

  allocate( gr(disc) )

  gr = rho / N / 2 / pi / rho_av / dr / rs**2

end function

pure function trapz(x, y) result(r)
    !! Calculates the integral of an array y with respect to x using the trapezoid
    !! approximation. Note that the mesh spacing of x does not have to be uniform.
    real(dp), intent(in)  :: x(:)         !! Variable x
    real(dp), intent(in)  :: y(size(x))   !! Function y(x)
    real(dp)              :: r            !! Integral ∫y(x)·dx

    ! Integrate using the trapezoidal rule
    associate(n => size(x))
      r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
    end associate
end function


pure function Sq_debeye(gr, rs, rho_av, disc, q) result(Sq)
    real(dp), intent(in), allocatable         :: gr(:), rs(:), q(:)
    real(dp), intent(in)         							:: rho_av
    integer, intent(in)                       :: disc
    real(dp), allocatable        :: sq(:)
    real(dp), allocatable        :: sqs(:,:)
    integer                                   :: i, j 


    allocate( sq(disc), sqs(disc, ubound(rs, dim=1)) )

    do i = 1, disc
			Sq(i) = trapz(	rs, rs * (gr-1) * sin(q(i)*rs)/q(i)	)
		end do

		Sq = 1 + 4*pi*rho_av*Sq
end function

subroutine linspace(from, to, array)
    real(dp), intent(in) :: from, to
    real(dp), intent(out) :: array(:)
    real(dp) :: range
    integer :: n, i
    n = size(array)
    range = to - from

    if (n == 0) return

    if (n == 1) then
        array(1) = from
        return
    end if


    do i=1, n
        array(i) = from + range * (i - 1) / (n - 1)
    end do
end subroutine

subroutine print_long_vec(vec)
    implicit none
    real(dp), intent(in), allocatable :: vec(:)
    integer :: row, col

    do row=1, ubound(vec, dim=1)/10
           write(*,10)  (vec(((row-1)*10)+col),col=0,9)
    end do
        10  format(15f10.3,/)
    write(*,'(A,/)')

end subroutine 

!subroutine read_xyz(f, pos)
!    use numeric_kinds
!    use atoms_module,            only: atoms, finalise, &
!                                     set_cutoff, &
!                                     calc_connect, print, &
!                                     map_into_cell, has_property
!    use cinoutput_module,        only : read, write
!
!    implicit none
!    real(dp), intent(out), allocatable :: pos(:,:)
!    real(dp), intent(out) :: lattice
!    character :: f(100), dat_str(1000)
!    integer :: ios, N
!
!
!    open(1, file=f, status='old', iostat=ios)
!
!    read(1, *) N
!    read(1, *) lattice
!
!
!
!    close(1, f)
end module functions

program main
    use numeric_kinds
    use functions
    use cinoutput_module,        only : read, write
    use atoms_module,            only: atoms, finalise, &
                                     set_cutoff, &
                                     calc_connect, print, &
                                     map_into_cell, has_property
    use dictionary_module,       only: dictionary, initialise, finalise, &
                                     set_value, STRING_LENGTH
    implicit none

    character(len=STRING_LENGTH) :: Library, &
                                   xyz_file

    integer :: real_disc, k_disc, row, col
    real(dp) :: r_out, a(3), b(3), lattice(3), ma, bw, rho_av
    real(dp), allocatable    :: bmp(:), rho(:), gr(:), q(:), Sq(:)!, pos(:,:) ! bin mid-points
    type(Atoms) :: my_atoms

    !! Set input parameters !!
    xyz_file = 'Si1k.xyz'
    real_disc = 200
    k_disc = 200

    print *, 'Reading in XYZ file'
    call read(my_atoms,trim(xyz_file))
    print *, my_atoms%pos(:,2)
    !print *, my_atoms%pos(2,:)

    lattice = [my_atoms%lattice(1,1), my_atoms%lattice(2,2), my_atoms%lattice(3,3)]
		rho_av = my_atoms%N/lattice(1)/lattice(2)/lattice(3)

		allocate( gr(real_disc), bmp(real_disc), rho(real_disc), q(k_disc), sq(k_disc) )
		call linspace(0.1_dp, 1._dp, q)
		
    call rho_calc(my_atoms%pos, lattice, real_disc, ma, bw, bmp, rho)
		
		gr = gr_from_rho(rho, bmp, bw, rho_av, my_atoms%N, real_disc)
		sq = Sq_debeye(gr, bmp, rho_av, k_disc, q)

    call print_long_vec(rho)
    call print_long_vec(bmp)
    call print_long_vec(gr)
    call print_long_vec(sq)


end program

