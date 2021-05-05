
program main
    use numeric_kinds
    use mc_routines
		use mt_stream
		use utils
    use cinoutput_module,        only : read, write
    use atoms_module,            only: atoms, finalise, &
                                     set_cutoff, cell_volume, &
                                     calc_connect, print, &
                                     map_into_cell, has_property, set_lattice
    use dictionary_module,       only: dictionary, initialise, finalise, &
                                     set_value, STRING_LENGTH
    implicit none

    character(len=STRING_LENGTH) :: Library, &
                                    comp_atoms_file, init_atoms_file
    character(len=50)           :: FMT
    !FMT = 

    integer :: real_disc, k_disc, row, col, i,j,k, dump_freq, Nsteps, n, seed, &
								a_ct, a_up_ct, r_ct, r_hs_ct, q_fit(2)
    real(dp) :: r_out, a(3), b(3), lattice(3), init_lattice(3), comp_lattice(3), &
                ma, bw, rho_av, init_rho_av, ra_rat, &
                k_fit_region(2), step_size, hard_sphere, sigma_accept, &
                lattice_scale_factor, V, rand, chi, chi_t, new_pos(3)
    real(dp), allocatable    :: bmp(:), gr(:), ngr(:), gr_comp(:), q(:), Sq(:), nsq(:), &
																Sq_comp(:), rand_pos(:,:), rnorms(:), k_weights(:)
    integer, allocatable     :: rho(:), rho_comp(:), drho(:), rand_at(:), tmp(:), n_rho(:)
    type(Atoms) :: comp_atoms, init_atoms, curr_atoms
		logical :: flag
		logical, allocatable :: q_mask(:)

    !! Set input parameters !!
    comp_atoms_file = '/p_2/Documents/MoS/RMC/JDM_RMC/examples/Si1k.xyz'
    init_atoms_file = '/p_2/Documents/MoS/RMC/JDM_RMC/examples/liq_1500k_216_GAP18.xyz'
    real_disc = 500
    k_disc = 500
    k_fit_region = [0.5, 15.]
    step_size = 0.5/sqrt(3.)
    hard_sphere = 0.5
    sigma_accept = 0.06
    dump_freq = 1000
		Nsteps = 1e6
!		k_weights = 1
    !----------------------------------------------------------!



    !! System-specific initialisation !!

    print *, 'Reading in XYZ file'
    call read(comp_atoms, trim(comp_atoms_file))
    call read(init_atoms, trim(init_atoms_file))
    call read(curr_atoms, trim(init_atoms_file))
    print *, comp_atoms%pos(:,2)
    !print *, comp_atoms%pos(2,:)
		call write(curr_atoms, 'examples/ats.xyz', append=.False.) ! write initial to xyz (and delete any old)

    comp_lattice = [comp_atoms%lattice(1,1), comp_atoms%lattice(2,2), comp_atoms%lattice(3,3)]
    init_lattice = [init_atoms%lattice(1,1), init_atoms%lattice(2,2), init_atoms%lattice(3,3)]
		rho_av = comp_atoms%N / cell_volume(comp_atoms)

    lattice_scale_factor = (rho_av * cell_volume(init_atoms)/init_atoms%N)**(-1./3.)
    print *, 'lattice scale factor: ', lattice_scale_factor
!    print *, 'initial lattice:      ', init_atoms%lattice
!		print *, 'old cell volume: ', cell_volume(init_atoms)
		init_rho_av = init_atoms%N / cell_volume(init_atoms)
		print *, 'Init rho_av: ', init_rho_av
		print *, 'scaled lattice: ', (init_atoms%lattice * lattice_scale_factor)
    call set_lattice(init_atoms, (init_atoms%lattice * lattice_scale_factor), .True.)
		init_rho_av = init_atoms%N / cell_volume(init_atoms)
!    print *, 'scaled lattice:      ', init_atoms%lattice
!		print *, 'new cell volume: ', cell_volume(init_atoms)

		print *, 'Comp rho_av: ', rho_av, 'Init rho_av: ', init_rho_av

    lattice = [init_atoms%lattice(1,1), init_atoms%lattice(2,2), init_atoms%lattice(3,3)]
    comp_lattice = [comp_atoms%lattice(1,1), comp_atoms%lattice(2,2), comp_atoms%lattice(3,3)]

		allocate( gr_comp(real_disc), gr(real_disc), &
              rho(real_disc), & 
							q(k_disc), sq(k_disc), sq_comp(k_disc), k_weights(k_disc), q_mask(k_disc) )

		call linspace(0.05_dp, 20._dp, q)
		
		q_mask = .False.
		where (q > k_fit_region(1) .and. q < k_fit_region(2)) q_mask = .True.
			
		tmp =  minloc( q, mask=q_mask )
		q_fit(1) = tmp(1)
		tmp =  maxloc( q, mask=q_mask )
 		q_fit(2) = tmp(1)

		print *, 'q fit indices', q_fit
		k_weights = 1
		
    call rho_calc(comp_atoms%pos, comp_lattice, real_disc, ma, bw, bmp, rho_comp)
!     write (*,5) 'x', 'y', 'z'
!    do i =1, size(comp_atoms%pos, 2) 
!      write (*,11) comp_atoms%pos(:,i)
!    enddo
!		write (*,*) 'ma:', ma, 'bw:', bw, 'N: ', comp_atoms%N, 'rho_av: ', rho_av
		
		gr_comp = gr_from_rho(rho_comp, bmp, bw, rho_av, comp_atoms%N)
		sq_comp = Sq_debeye(gr_comp, bmp, rho_av, q)

    open (2, file='examples/gr_comp.dat')
      write (2,5) 'r/A', 'g(r)'
    do i =1, size(bmp) 
      write (2,7) bmp(i), gr_comp(i)
    enddo
    close(2)
    open (2, file='examples/rho_comp.dat')
      write (2,5) 'r/A', 'g(r)'
    do i =1, size(bmp) 
      write (2,9) bmp(i), rho_comp(i)
    enddo
    close(2)
    open (2, file='examples/sq_comp.dat')
      write (2,5) 'q/A^-1', 'S(q)'
    do i =1, size(q) 
      write (2,7) q(i), Sq_comp(i)
    enddo
    close(2)


    call rho_calc(init_atoms%pos, lattice, real_disc, ma, bw, bmp, rho)
		gr = gr_from_rho(rho, bmp, bw, rho_av, init_atoms%N)
		sq = Sq_debeye(gr, bmp, rho_av, q)
		chi = chi_squared(sq, sq_comp, k_weights, q_fit)
    print *, 'initial chi: ', chi

!    call print_long_vec(real(rho, 8))
!    print *, rho
!    call print_long_vec(bmp)
!    call print_long_vec(gr)
!    call print_long_vec(sq)
    open (2, file='examples/sq_init.dat')
      write (2,5) 'q/A^-1', 'S(q)'
    do i =1, size(q) 
      write (2,7) q(i), Sq(i)
    enddo
    close(2)

    open (2, file='examples/gr_init.dat')
      write (2,5) 'r/A', 'g(r)'
    do i =1, size(bmp) 
      write (2,7) bmp(i), gr(i)
    enddo
    close(2)
 5  format( (2(A15)) )
 7  format( (2(E15.5)) )
 9  format( ( (E15.5) (I15) ) )
 11  format( (3(E15.5)) )
 13  format( (E15.7) )

! ---------------- MAIN MC moves begin here ------------------- !	


!n = 1
!call random_seed(size=n)
!allocate (seed(n))
!allocate (rput(n))
!rput=1
!call random_seed(put=rput)
!call random_seed(get=seed)
rand_pos = random_positions(Nsteps, step_size, seed)
rand_at = random_ats(Nsteps, seed, curr_atoms%N)
rnorms = random_normals(Nsteps, seed, sigma_accept)
a_ct = 0
a_up_ct = 0
r_ct = 0
r_hs_ct = 0
ra_rat = 1
allocate( drho(real_disc) )
open(2, file='examples/grs.dat')
	write (2,*) 'r/A then g(r)s'
	write (2,*) bmp
close(2)
open(2, file='examples/sqs.dat')
	write (2,*) 'q/A then S(q)s'
close(2)
open(2, file='examples/chis.dat')
	write (2,*) 'Xs'
close(2)
open(2, file='examples/stats.dat')
	write (2, "(5A15)") 'a_ct', 'up_a_ct', 'r_ct', 'r_hs_ct', 'ratio' 
	write (2,"(4I15,F4.2)") a_ct, a_up_ct, r_ct, r_hs_ct, ra_rat 
close(2)

do i=1, Nsteps
 ! print *, 'step ', i	
	! write to dat files to keep track of simulation
	if ( mod(i, dump_freq) .eq. 0 ) then
		call write(curr_atoms, 'examples/ats.xyz', append=.True.)

		open(2, file='examples/grs.dat', access='append')
			write (2,*) gr
    close(2)

		open(2, file='examples/sqs.dat', access='append')
			write (2,*) sq
    close(2)

		open(2, file='examples/chis.dat', access='append')
			write (2,13) chi
    close(2)
    
    ra_rat = a_ct/(a_ct+r_ct)
		open(2, file='examples/stats.dat', access='append')
			write (2,"(5I15)") a_ct, a_up_ct, r_ct, r_hs_ct 
    close(2)
	endif
   
  new_pos = rand_pos(:,i) + curr_atoms%pos(:, rand_at(i))
!  print *, 'old pos: ', curr_atoms%pos(:, rand_at(i))
!  print *, 'new pos: ', new_pos
	call rho_diff(drho, flag, new_pos, rand_at(i), curr_atoms%pos, lattice, bmp, bw, hard_sphere)
  n_rho = rho+drho
  !write(*, '(I50)') ( n_rho(j), j=1,size(rho))
  if (.not. flag ) then
    r_hs_ct = r_hs_ct + 1
    cycle
  end if

	ngr = gr_from_rho(n_rho, bmp, bw, rho_av, curr_atoms%N)
	nsq = Sq_debeye(ngr, bmp, rho_av, q)
  
  if ( i .eq. 1 ) then
  open (2, file='examples/sq_moved.dat')
    write (2,5) 'q/A^-1', 'S(q)'
  do j =1, size(q) 
    write (2,7) q(j), nsq(j)
  enddo
  close(2)
  end if

	chi_t = chi_squared(nsq, sq_comp, k_weights, q_fit)
  print*, 'chi_t: ', chi_t

	if ( chi_t < chi ) then
    
    where ( new_pos > lattice ) ! map back into lattice if strayed outside
      new_pos = new_pos - lattice
    elsewhere ( new_pos < 0 )
      new_pos = new_pos + lattice
    end where

    curr_atoms%pos(:, rand_at(i)) = new_pos
    rho = rho + drho
    gr = ngr
    sq = nsq
    chi = chi_t
    a_ct = a_ct+1

  else if ( chi_t - chi < rnorms(i) ) then
    
    where ( new_pos > lattice ) ! map back into lattice if strayed outside
      new_pos = new_pos - lattice
    elsewhere ( new_pos < 0 )
      new_pos = new_pos + lattice
    end where

    curr_atoms%pos(:, rand_at(i)) = new_pos
    rho = rho + drho
    gr = gr
    sq = sq
    chi = chi_t
    a_ct = a_ct+1
    a_up_ct = a_up_ct+1

  else
    r_ct = r_ct + 1
		
  end if		


end do


contains
	function random_positions(Nsteps, step_size, seed) result(new_pos)
		
		real(dp), allocatable			:: new_pos(:,:)
		real(dp), intent(in)									:: step_size
		integer, intent(in)				:: Nsteps, seed
		type(mt_state) :: mts
		call set_mt19937
		call new(mts)
		call init(mts, seed)
		call print(mts)
		allocate ( new_pos(3, Nsteps) )
!		rand = genrand_double1(mts)

!		print *, rand

		do i=1,Nsteps
			new_pos(:,i) = [ (genrand_double1(mts), i=1,3) ] * step_size
		end do

		print *, new_pos(:,1:5)
	end function random_positions

	function random_ats(Nsteps, seed, N) result(ats)
		
		integer, allocatable			:: ats(:)
    real, allocatable         :: t_ats(:)
		integer, intent(in)				:: Nsteps, seed, N
    integer                   :: i
		type(mt_state) :: mts
		call set_mt19937
		call new(mts)
		call init(mts, seed)
		call print(mts)
		allocate ( ats(Nsteps) )
    
		t_ats = [ (genrand_double2(mts), i=1,Nsteps) ] ! r in [0,1)
    ats = floor(t_ats*(N)) + 1

	end function random_ats

	function random_normals(Nsteps, seed, sigma) result(normals)
		
		real(dp), allocatable			:: normals(:), rands(:,:)
		integer, intent(in)				:: Nsteps, seed
    real(dp), intent(in)      :: sigma
		type(mt_state) :: mts
		call set_mt19937
		call new(mts)
		call init(mts, seed)
		call print(mts)
		allocate ( normals(Nsteps), rands(2,Nsteps) )

		do i=1,Nsteps
			rands(:,i) = [ (genrand_double1(mts), i=1,2) ]
		end do
    normals = sin(2*pi*rands(1,:)) * sqrt(-2*log(rands(2,:))) * sigma

	end function random_normals

end program

