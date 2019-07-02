!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the Wannier90 code and !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the Wannier90      !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the Wannier90 code is www.wannier.org       !
!                                                            !
! The Wannier90 code is hosted on GitHub:                    !
!                                                            !
! https://github.com/wannier-developers/wannier90            !
!------------------------------------------------------------!

module w90_hamiltonian
  !! Module to obtain the Hamiltonian in a wannier basis
  !! This is a simplified routine, more sophisticated properties
  !! are found in postw90 (e.g. w90_get_oper)
  use w90_constants, only: dp
  use w90_comms, only: on_root

  implicit none

  private
  !
  complex(kind=dp), public, save, allocatable :: ham_r(:, :, :)
  !! Hamiltonian matrix in WF representation
  !
  integer, public, save, allocatable :: irvec(:, :)
  !!  The irpt-th Wigner-Seitz grid point has components
  !! irvec(1:3,irpt) in the basis of the lattice vectors
  !
  integer, public, save, allocatable :: shift_vec(:, :)
  !
  integer, public, save, allocatable :: ndegen(:)
  !! Weight of the irpt-th point is 1/ndegen(irpt)
  !
  integer, public, save              :: nrpts
  !! number of Wigner-Seitz grid points
  !
  integer, public, save              :: rpt_origin
  !! index of R=0
  !
  real(kind=dp), public, save, allocatable :: wannier_centres_translated(:, :)
  !! translated Wannier centres
  !
  complex(kind=dp), public, save, allocatable :: spn_r(:, :, :, :) !jmlihm
  !! Hamiltonian matrix in WF representation

  public :: hamiltonian_get_hr
  public :: hamiltonian_write_hr
  public :: hamiltonian_setup
  public :: hamiltonian_dealloc
  public :: hamiltonian_write_rmn
  public :: hamiltonian_write_tb
  public :: hamiltonian_get_spnr !jmlihm
  public :: hamiltonian_write_spnr !jmlihm

  ! Module variables
  logical, save :: ham_have_setup = .false.
  logical, save :: have_translated = .false.
  logical, save :: use_translation = .false.
  logical, save :: have_ham_r = .false.
  logical, save :: have_ham_k = .false.
  logical, save :: hr_written = .false.
  logical, save :: tb_written = .false.
  logical, save :: spnr_written = .false. !jmlihm

  complex(kind=dp), save, allocatable :: ham_k(:, :, :)

contains

  !============================================!
  subroutine hamiltonian_setup()
    !! Allocate arrays and setup data
    !============================================!

    use w90_constants, only: cmplx_0
    use w90_io, only: io_error
    use w90_parameters, only: num_wann, num_kpts, bands_plot, transport, &
      bands_plot_mode, transport_mode

    implicit none

    integer :: ierr

    if (ham_have_setup) return

    !
    ! Determine whether to use translation
    !
    if (bands_plot .and. (index(bands_plot_mode, 'cut') .ne. 0)) use_translation = .true.
    if (transport .and. (index(transport_mode, 'bulk') .ne. 0)) use_translation = .true.
    if (transport .and. (index(transport_mode, 'lcr') .ne. 0)) use_translation = .true.
    !
    ! Set up Wigner-Seitz vectors
    !
    call hamiltonian_wigner_seitz(count_pts=.true.)
    !
    allocate (irvec(3, nrpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating irvec in hamiltonian_setup')
    irvec = 0
    !
    allocate (ndegen(nrpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ndegen in hamiltonian_setup')
    ndegen = 0
    !
    allocate (ham_r(num_wann, num_wann, nrpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ham_r in hamiltonian_setup')
    ham_r = cmplx_0
    !
    allocate (ham_k(num_wann, num_wann, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ham_k in hamiltonian_setup')
    ham_k = cmplx_0
    !
    ! Set up the wigner_seitz vectors
    !
    call hamiltonian_wigner_seitz(count_pts=.false.)
    !
    allocate (wannier_centres_translated(3, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating wannier_centres_translated in hamiltonian_setup')
    wannier_centres_translated = 0.0_dp

    ham_have_setup = .true.

    return
  end subroutine hamiltonian_setup

  !============================================!
  subroutine hamiltonian_dealloc()
    !! Deallocate module data
    !============================================!

    use w90_io, only: io_error

    implicit none

    integer :: ierr

    if (allocated(ham_r)) then
      deallocate (ham_r, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating ham_r in hamiltonian_dealloc')
    end if
    if (allocated(ham_k)) then
      deallocate (ham_k, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating ham_k in hamiltonian_dealloc')
    end if
    if (allocated(irvec)) then
      deallocate (irvec, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating irvec in hamiltonian_dealloc')
    end if
    if (allocated(ndegen)) then
      deallocate (ndegen, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating ndegen in hamiltonian_dealloc')
    end if
    if (allocated(wannier_centres_translated)) then
      deallocate (wannier_centres_translated, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating wannier_centres_translated in param_dealloc')
    end if

    return
  end subroutine hamiltonian_dealloc

  !============================================!
  subroutine hamiltonian_get_hr()
    !============================================!
    !                                            !
    !!  Calculate the Hamiltonian in the WF basis
    !                                            !
    !============================================!

    use w90_constants, only: cmplx_0, cmplx_i, twopi
    use w90_io, only: io_error, io_stopwatch
    use w90_parameters, only: num_bands, num_kpts, num_wann, u_matrix, &
      eigval, kpt_latt, u_matrix_opt, lwindow, ndimwin, &
      have_disentangled, timing_level
    use w90_parameters, only: lsitesymmetry !YN:

    implicit none

    integer, allocatable :: shift_vec(:, :)
    complex(kind=dp)     :: fac
    real(kind=dp)        :: rdotk
    real(kind=dp)        :: eigval_opt(num_bands, num_kpts)
    real(kind=dp)        :: eigval2(num_wann, num_kpts)
    real(kind=dp)        :: irvec_tmp(3)
    integer              :: loop_kpt, i, j, m, irpt, ideg, ierr, counter
    complex(kind=dp)     :: utmp(num_bands, num_wann) !RS:

    if (timing_level > 1) call io_stopwatch('hamiltonian: get_hr', 1)

    if (have_ham_r) then
      if (have_translated .eqv. use_translation) then
        goto 200
      else
        goto 100
      endif
    end if

    if (have_ham_k) go to 100

!~    if (.not. allocated(ham_k)) then
!~       allocate(ham_k(num_wann,num_wann,num_kpts),stat=ierr)
!~       if (ierr/=0) call io_error('Error in allocating ham_k in hamiltonian_get_hr')
!~    end if

    ham_k = cmplx_0
    eigval_opt = 0.0_dp
    eigval2 = 0.0_dp

    if (have_disentangled) then

      ! slim down eigval to contain states within the outer window

      do loop_kpt = 1, num_kpts
        counter = 0
        do j = 1, num_bands
          if (lwindow(j, loop_kpt)) then
            counter = counter + 1
            eigval_opt(counter, loop_kpt) = eigval(j, loop_kpt)
          end if
        end do
      end do

      ! rotate eigval into the optimal subspace
      ! in general eigval would be a matrix at each kpoints
      ! but we choose u_matrix_opt such that the Hamiltonian is
      ! diagonal at each kpoint. (I guess we should check it here)

      if (.not. lsitesymmetry) then                                                                             !YN:
        do loop_kpt = 1, num_kpts
          do j = 1, num_wann
            do m = 1, ndimwin(loop_kpt)
              eigval2(j, loop_kpt) = eigval2(j, loop_kpt) + eigval_opt(m, loop_kpt)* &
                                     real(conjg(u_matrix_opt(m, j, loop_kpt))*u_matrix_opt(m, j, loop_kpt), dp)
            enddo
          enddo
        enddo
      else                                                                                                     !YN:
        ! u_matrix_opt are not the eigenvectors of the Hamiltonian any more                                   !RS:
        ! so we have to calculate ham_k in the following way                                                  !RS:
        do loop_kpt = 1, num_kpts                                                                                !RS:
          utmp(1:ndimwin(loop_kpt), :) = &                                                                     !RS:
            matmul(u_matrix_opt(1:ndimwin(loop_kpt), :, loop_kpt), u_matrix(:, :, loop_kpt))                    !RS:
          do j = 1, num_wann                                                                                    !RS:
            do i = 1, j                                                                                        !RS:
              do m = 1, ndimwin(loop_kpt)                                                                     !RS:
                ham_k(i, j, loop_kpt) = ham_k(i, j, loop_kpt) + eigval_opt(m, loop_kpt)*conjg(utmp(m, i))*utmp(m, j) !RS:
              enddo                                                                                        !RS:
              if (i .lt. j) ham_k(j, i, loop_kpt) = conjg(ham_k(i, j, loop_kpt))                                   !RS:
            enddo                                                                                           !RS:
          enddo                                                                                              !RS:
        enddo                                                                                                 !RS:
      endif                                                                                                    !YN:

    else
      eigval2(1:num_wann, :) = eigval(1:num_wann, :)
    end if

    ! At this point eigval2 contains num_wann values which belong to the wannier subspace.

    ! Rotate Hamiltonian into the basis of smooth bloch states
    !          H(k)=U^{dagger}(k).H_0(k).U(k)
    ! Note: we enforce hermiticity here

    if (.not. lsitesymmetry .or. .not. have_disentangled) then !YN:
      do loop_kpt = 1, num_kpts
        do j = 1, num_wann
          do i = 1, j
            do m = 1, num_wann
              ham_k(i, j, loop_kpt) = ham_k(i, j, loop_kpt) + eigval2(m, loop_kpt)* &
                                      conjg(u_matrix(m, i, loop_kpt))*u_matrix(m, j, loop_kpt)
            enddo
            if (i .lt. j) ham_k(j, i, loop_kpt) = conjg(ham_k(i, j, loop_kpt))
          enddo
        enddo
      enddo
    endif                                                  !YN:

    have_ham_k = .true.

100 continue

    ! Fourier transform rotated hamiltonian into WF basis
    ! H_ij(k) --> H_ij(R) = (1/N_kpts) sum_k e^{-ikR} H_ij(k)
!~    if (.not.allocated(ham_r)) then
!~      allocate(ham_r(num_wann,num_wann,nrpts),stat=ierr)
!~      if (ierr/=0) call io_error('Error in allocating ham_r in hamiltonian_get_hr')
!~    end if

    ham_r = cmplx_0

    if (.not. use_translation) then

      do irpt = 1, nrpts
        do loop_kpt = 1, num_kpts
          rdotk = twopi*dot_product(kpt_latt(:, loop_kpt), real(irvec(:, irpt), dp))
          fac = exp(-cmplx_i*rdotk)/real(num_kpts, dp)
          ham_r(:, :, irpt) = ham_r(:, :, irpt) + fac*ham_k(:, :, loop_kpt)
        enddo
      enddo

      have_translated = .false.

    else

      allocate (shift_vec(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating shift_vec in hamiltonian_get_hr')
      call internal_translate_centres()

      do irpt = 1, nrpts
        do loop_kpt = 1, num_kpts
          do i = 1, num_wann
            do j = 1, num_wann
              ! ham_r(j,i,irpt)
              ! interaction btw j at 0 and i at irvec(:,irpt)
              irvec_tmp(:) = irvec(:, irpt) + shift_vec(:, i) - shift_vec(:, j)
              rdotk = twopi*dot_product(kpt_latt(:, loop_kpt), real(irvec_tmp(:), dp))
              fac = exp(-cmplx_i*rdotk)/real(num_kpts, dp)
              ham_r(j, i, irpt) = ham_r(j, i, irpt) + fac*ham_k(j, i, loop_kpt)
            end do
          end do
        enddo
      enddo

      have_translated = .true.

    end if

    ! [lp] if required, compute the minimum diistances
!     if (use_ws_distance) then
!         allocate(irdist_ws(3,ndegenx,num_wann,num_wann,nrpts),stat=ierr)
!         if (ierr/=0) call io_error('Error in allocating irdist_ws in hamiltonian_get_hr')
!         allocate(wdist_ndeg(num_wann,num_wann,nrpts),stat=ierr)
!         if (ierr/=0) call io_error('Error in allocating wcenter_ndeg in hamiltonian_get_hr')
    !
!         call ws_translate_dist(nrpts, irvec)
!     endif

    have_ham_r = .true.

200 continue

    if (allocated(shift_vec)) then
      deallocate (shift_vec, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating shift_vec in hamiltonian_get_hr')
    end if

    if (timing_level > 1) call io_stopwatch('hamiltonian: get_hr', 2)

    return

  contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !====================================================!
    subroutine internal_translate_centres()
      !! Translate the centres of the WF into the home cell
      !====================================================!

      use w90_parameters, only: num_wann, real_lattice, recip_lattice, wannier_centres, &
        num_atoms, atoms_pos_cart, translation_centre_frac, &
        automatic_translation, num_species, atoms_species_num, lenconfac
      use w90_io, only: stdout, io_error
      use w90_utility, only: utility_cart_to_frac, utility_frac_to_cart

      implicit none

      ! <<<local variables>>>
      integer :: iw, ierr, nat, nsp, ind
      real(kind=dp), allocatable :: r_home(:, :), r_frac(:, :)
      real(kind=dp) :: c_pos_cart(3), c_pos_frac(3)
      real(kind=dp) :: r_frac_min(3)

!~      if (.not.allocated(wannier_centres_translated)) then
!~         allocate(wannier_centres_translated(3,num_wann),stat=ierr)
!~         if (ierr/=0) call io_error('Error in allocating wannier_centres_translated &
!~              &in internal_translate_wannier_centres')
!~      end if

      allocate (r_home(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating r_home in internal_translate_centres')
      allocate (r_frac(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating r_frac in internal_translate_centres')
      r_home = 0.0_dp; r_frac = 0.0_dp

      if (automatic_translation) then
        ! Calculate centre of atomic positions
        c_pos_cart = 0.0_dp; c_pos_frac = 0.0_dp
        do nsp = 1, num_species
          do nat = 1, atoms_species_num(nsp)
            c_pos_cart(:) = c_pos_cart(:) + atoms_pos_cart(:, nat, nsp)
          enddo
        enddo
        c_pos_cart = c_pos_cart/num_atoms
        ! Cartesian --> fractional
        call utility_cart_to_frac(c_pos_cart, translation_centre_frac, recip_lattice)
      end if
      ! Wannier function centres will be in [c_pos_frac-0.5,c_pos_frac+0.5]
      r_frac_min(:) = translation_centre_frac(:) - 0.5_dp

      ! Cartesian --> fractional
      do iw = 1, num_wann
        call utility_cart_to_frac(wannier_centres(:, iw), r_frac(:, iw), recip_lattice)
        ! Rationalise r_frac - r_frac_min to interval [0,1]
        !  by applying shift of -floor(r_frac - r_frac_min)
        shift_vec(:, iw) = -floor(r_frac(:, iw) - r_frac_min(:))
        r_frac(:, iw) = r_frac(:, iw) + real(shift_vec(:, iw), dp)
        ! Fractional --> Cartesian
        call utility_frac_to_cart(r_frac(:, iw), r_home(:, iw), real_lattice)
      end do

      ! NEVER overwrite wannier_centres
      !wannier_centres = r_home

      if (on_root) then
        write (stdout, '(1x,a)') 'Translated centres'
        write (stdout, '(4x,a,3f10.6)') 'translation centre in fractional coordinate:', translation_centre_frac(:)
        do iw = 1, num_wann
          write (stdout, 888) iw, (r_home(ind, iw)*lenconfac, ind=1, 3)
        end do
        write (stdout, '(1x,a78)') repeat('-', 78)
        write (stdout, *)
      endif
      wannier_centres_translated = r_home

      deallocate (r_frac, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating r_frac in internal_translate_centres')
      deallocate (r_home, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating r_home in internal_translate_centres')

      return

888   format(2x, 'WF centre ', i5, 2x, '(', f10.6, ',', f10.6, ',', f10.6, ' )')

    end subroutine internal_translate_centres

  end subroutine hamiltonian_get_hr

  !============================================!
  subroutine hamiltonian_write_hr()
    !============================================!
    !!  Write the Hamiltonian in the WF basis
! jmlihm: write as binary file, not text file
    !============================================!

    use w90_io, only: io_error, io_stopwatch, io_file_unit, &
      seedname, io_date
    use w90_parameters, only: num_wann, timing_level, wannier_centres, real_lattice, &
      do_write_bin, do_write_text

    integer            :: i, j, irpt, file_unit
    integer :: reclen
    character(len=33) :: header
    character(len=9)  :: cdate, ctime
    integer :: reclen

    if (hr_written) return

    if (timing_level > 1) call io_stopwatch('hamiltonian: write_hr', 1)

    ! write the  whole matrix with all the indices

    ! write text file
    if (do_write_text) then
      file_unit = io_file_unit()
      open (file_unit, file=trim(seedname)//'_hr.dat', form='formatted', &
            status='unknown', err=101)

      call io_date(cdate, ctime)
      header = 'written on '//cdate//' at '//ctime

      write (file_unit, *) header ! Date and time
      write (file_unit, *) num_wann
      write (file_unit, *) nrpts
      write (file_unit, '(15I5)') (ndegen(i), i=1, nrpts)
      do irpt = 1, nrpts
        do i = 1, num_wann
          do j = 1, num_wann
            write (file_unit, '(5I5,2F12.6)') irvec(:, irpt), j, i, &
              ham_r(j, i, irpt)
          end do
        end do
      end do
      close (file_unit)
    end if

    ! write binary file
    if (do_write_bin) then
      file_unit = io_file_unit()
      inquire (iolength=reclen) ndegen(:)
      open (file_unit, file=trim(seedname)//'_ndegen.bin', form='unformatted', &
            status='unknown', access='direct', recl=reclen, err=101)
      write (file_unit, rec=1) ndegen
      close (file_unit)

      inquire (iolength=reclen) irvec(:, :)
      open (file_unit, file=trim(seedname)//'_irvec.bin', form='unformatted', &
            status='unknown', access='direct', recl=reclen, err=101)
      write (file_unit, rec=1) irvec
      close (file_unit)

      inquire (iolength=reclen) ham_r(:, :, :)
      open (file_unit, file=trim(seedname)//'_hr.bin', form='unformatted', &
            status='unknown', access='direct', recl=reclen, err=101)
      write (file_unit, rec=1) ham_r
      close (file_unit)

      inquire (iolength=reclen) wannier_centres(:, :)
      open (file_unit, file=trim(seedname)//'_wcenter.bin', form='unformatted', &
            status='unknown', access='direct', recl=reclen, err=101)
      write (file_unit, rec=1) wannier_centres
      close (file_unit)

      inquire (iolength=reclen) real_lattice(:, :)
      open (file_unit, file=trim(seedname)//'_reallatt.bin', form='unformatted', &
            status='unknown', access='direct', recl=reclen, err=101)
      write (file_unit, rec=1) real_lattice
      close (file_unit)
    end if

    hr_written = .true.

    if (timing_level > 1) call io_stopwatch('hamiltonian: write_hr', 2)

    return

101 call io_error('Error: hamiltonian_write_hr: problem opening file '//trim(seedname)//'_hr.dat')

  end subroutine hamiltonian_write_hr

  !================================================================================!
  subroutine hamiltonian_wigner_seitz(count_pts)
    !================================================================================!
    !! Calculates a grid of points that fall inside of (and eventually on the
    !! surface of) the Wigner-Seitz supercell centered on the origin of the B
    !! lattice with primitive translations nmonkh(1)*a_1+nmonkh(2)*a_2+nmonkh(3)*a_3
    !================================================================================!

    use w90_constants, only: eps7, eps8
    use w90_io, only: io_error, io_stopwatch, stdout
    use w90_parameters, only: iprint, mp_grid, real_metric, timing_level

    ! irvec(i,irpt)     The irpt-th Wigner-Seitz grid point has components
    !                   irvec(1:3,irpt) in the basis of the lattice vectors
    ! ndegen(irpt)      Weight of the irpt-th point is 1/ndegen(irpt)
    ! nrpts             number of Wigner-Seitz grid points

    implicit none

    logical, intent(in) :: count_pts
    !! Only count points and return

    integer       :: ndiff(3)
    real(kind=dp) :: dist(125), tot, dist_min
    integer       :: n1, n2, n3, i1, i2, i3, icnt, i, j

    if (timing_level > 1) call io_stopwatch('hamiltonian: wigner_seitz', 1)

    ! The Wannier functions live in a supercell of the real space unit cell
    ! this supercell is mp_grid unit cells long in each direction
    !
    ! We loop over grid points r on a unit cell that is 8 times larger than this
    ! primitive supercell.
    !
    ! One of these points is in the W-S cell if it is closer to R=0 than any of the
    ! other points, R (where R are the translation vectors of the supercell)

    ! In the end nrpts contains the total number of grid
    ! points that have been found in the Wigner-Seitz cell

    nrpts = 0
    do n1 = -mp_grid(1), mp_grid(1)
      do n2 = -mp_grid(2), mp_grid(2)
        do n3 = -mp_grid(3), mp_grid(3)
          ! Loop over the 125 points R. R=0 corresponds to
          ! i1=i2=i3=0, or icnt=63
          icnt = 0
          do i1 = -2, 2
            do i2 = -2, 2
              do i3 = -2, 2
                icnt = icnt + 1
                ! Calculate distance squared |r-R|^2
                ndiff(1) = n1 - i1*mp_grid(1)
                ndiff(2) = n2 - i2*mp_grid(2)
                ndiff(3) = n3 - i3*mp_grid(3)
                dist(icnt) = 0.0_dp
                do i = 1, 3
                  do j = 1, 3
                    dist(icnt) = dist(icnt) + real(ndiff(i), dp)*real_metric(i, j) &
                                 *real(ndiff(j), dp)
                  enddo
                enddo
              enddo
            enddo

            ! AAM: On first pass, we reference unallocated variables (ndegen,irvec)

          enddo
          dist_min = minval(dist)
          if (abs(dist(63) - dist_min) .lt. eps7) then
            nrpts = nrpts + 1
            if (.not. count_pts) then
              ndegen(nrpts) = 0
              do i = 1, 125
                if (abs(dist(i) - dist_min) .lt. eps7) ndegen(nrpts) = ndegen(nrpts) + 1
              end do
              irvec(1, nrpts) = n1
              irvec(2, nrpts) = n2
              irvec(3, nrpts) = n3
              !
              ! Record index of r=0
              if (n1 == 0 .and. n2 == 0 .and. n3 == 0) rpt_origin = nrpts
            endif
          end if

          !n3
        enddo
        !n2
      enddo
      !n1
    enddo
    !
    if (count_pts) return

    if (iprint >= 3 .and. on_root) then
      write (stdout, '(1x,i4,a,/)') nrpts, ' lattice points in Wigner-Seitz supercell:'
      do i = 1, nrpts
        write (stdout, '(4x,a,3(i3,1x),a,i2)') '  vector ', irvec(1, i), irvec(2, i), &
          irvec(3, i), '  degeneracy: ', ndegen(i)
      enddo
    endif
    ! Check the "sum rule"
    tot = 0.0_dp
    do i = 1, nrpts
      tot = tot + 1.0_dp/real(ndegen(i), dp)
    enddo
    if (abs(tot - real(mp_grid(1)*mp_grid(2)*mp_grid(3), dp)) > eps8) then
      call io_error('ERROR in hamiltonian_wigner_seitz: error in finding Wigner-Seitz points')
    endif

    if (timing_level > 1) call io_stopwatch('hamiltonian: wigner_seitz', 2)

    return

  end subroutine hamiltonian_wigner_seitz

  !============================================!
  subroutine hamiltonian_write_rmn()
    !! Write out the matrix elements of r
    !============================================!
    use w90_parameters, only: m_matrix, wb, bk, num_wann, num_kpts, kpt_latt, &
      nntot, write_bvec, do_write_bin, do_write_text
    use w90_constants, only: twopi, cmplx_i
    use w90_io, only: io_error, io_file_unit, seedname, io_date

    implicit none

    complex(kind=dp)     :: fac
    real(kind=dp)        :: rdotk
    integer              :: loop_rpt, m, n, nkp, ind, nn, file_unit
    complex(kind=dp)     :: position(3)
    complex(kind=dp), allocatable :: position_save(:, :, :, :)
    character(len=33) :: header
    character(len=9)  :: cdate, ctime
    integer :: reclen

    ! save positions for later use
    allocate (position_save(num_wann, num_wann, nrpts, 3))

    ! write text file
    if (do_write_text) then
      file_unit = io_file_unit()
      open (file_unit, file=trim(seedname)//'_r.dat', form='formatted', status='unknown', err=101)
      call io_date(cdate, ctime)

      header = 'written on '//cdate//' at '//ctime
      write (file_unit, *) header ! Date and time
      write (file_unit, *) num_wann
      write (file_unit, *) nrpts
    end if

    do loop_rpt = 1, nrpts
      do m = 1, num_wann
        do n = 1, num_wann
          position(:) = 0._dp
          do nkp = 1, num_kpts
            rdotk = twopi*dot_product(kpt_latt(:, nkp), real(irvec(:, loop_rpt), dp))
            fac = exp(-cmplx_i*rdotk)/real(num_kpts, dp)
            do ind = 1, 3
              do nn = 1, nntot
                if (m .eq. n) then
                  ! For loop_rpt==rpt_origin, this reduces to
                  ! Eq.(32) of Marzari and Vanderbilt PRB 56,
                  ! 12847 (1997). Otherwise, is is Eq.(44)
                  ! Wang, Yates, Souza and Vanderbilt PRB 74,
                  ! 195118 (2006), modified according to
                  ! Eqs.(27,29) of Marzari and Vanderbilt
                  position(ind) = position(ind) - &
                                  wb(nn)*bk(ind, nn, nkp)*aimag(log(m_matrix(n, m, nn, nkp)))*fac
                else
                  ! Eq.(44) Wang, Yates, Souza and Vanderbilt PRB 74, 195118 (2006)
                  position(ind) = position(ind) + &
                                  cmplx_i*wb(nn)*bk(ind, nn, nkp)*m_matrix(n, m, nn, nkp)*fac
                endif
              end do
            end do
          end do
          if (do_write_text) write (file_unit, '(5I5,6F12.6)') irvec(:, loop_rpt), n, m, position(:)
          position_save(n, m, loop_rpt, :) = position(:)
        end do
      end do
    end do
    if (do_write_text) close (file_unit)
!    end if

    if (do_write_bin) then
      file_unit = io_file_unit()
      inquire (iolength=reclen) position_save(:, :, :, :)
      open (file_unit, file=trim(seedname)//'_r.bin', form='unformatted', &
            status='unknown', access='direct', recl=reclen, err=101)
      write (file_unit, rec=1) position_save
      close (file_unit)
    end if

    if (allocated(position_save)) deallocate (position_save)

    return

101 call io_error('Error: hamiltonian_write_rmn: problem opening file '//trim(seedname)//'_r')

  end subroutine hamiltonian_write_rmn

  !============================================!
  subroutine hamiltonian_write_tb()
    !============================================!
    !! Write in a single file all the information
    !! that is needed to set up a Wannier-based
    !! tight-binding model:
    !! * lattice vectors
    !! * <0n|H|Rn>
    !! * <0n|r|Rn>
    !============================================!

    use w90_io, only: io_error, io_stopwatch, io_file_unit, &
      seedname, io_date
    use w90_parameters, only: real_lattice, num_wann, timing_level, &
      m_matrix, wb, bk, num_kpts, kpt_latt, nntot
    use w90_constants, only: twopi, cmplx_i

    integer            :: i, j, irpt, ik, nn, idir, file_unit
    character(len=33) :: header
    character(len=9)  :: cdate, ctime
    complex(kind=dp)   :: fac, pos_r(3)
    real(kind=dp)      :: rdotk

    if (tb_written) return

    if (timing_level > 1) call io_stopwatch('hamiltonian: write_tb', 1)

    file_unit = io_file_unit()
    open (file_unit, file=trim(seedname)//'_tb.dat', form='formatted', &
          status='unknown', err=101)

    call io_date(cdate, ctime)
    header = 'written on '//cdate//' at '//ctime

    write (file_unit, *) header ! Date and time
    !
    ! lattice vectors
    !
    write (file_unit, *) real_lattice(1, :) !a_1
    write (file_unit, *) real_lattice(2, :) !a_2
    write (file_unit, *) real_lattice(3, :) !a_3
    !
    write (file_unit, *) num_wann
    write (file_unit, *) nrpts
    write (file_unit, '(15I5)') (ndegen(i), i=1, nrpts)
    !
    ! <0n|H|Rm>
    !
    do irpt = 1, nrpts
      write (file_unit, '(/,3I5)') irvec(:, irpt)
      do i = 1, num_wann
        do j = 1, num_wann
          write (file_unit, '(2I5,3x,2(E15.8,1x))') j, i, ham_r(j, i, irpt)
        end do
      end do
    end do
    !
    ! <0n|r|Rm>
    !
    do irpt = 1, nrpts
      write (file_unit, '(/,3I5)') irvec(:, irpt)
      do i = 1, num_wann
        do j = 1, num_wann
          pos_r(:) = 0._dp
          do ik = 1, num_kpts
            rdotk = twopi*dot_product(kpt_latt(:, ik), real(irvec(:, irpt), dp))
            fac = exp(-cmplx_i*rdotk)/real(num_kpts, dp)
            do idir = 1, 3
              do nn = 1, nntot
                if (i == j) then
                  ! For irpt==rpt_origin, this reduces to
                  ! Eq.(32) of Marzari and Vanderbilt PRB 56,
                  ! 12847 (1997). Otherwise, is is Eq.(44)
                  ! Wang, Yates, Souza and Vanderbilt PRB 74,
                  ! 195118 (2006), modified according to
                  ! Eqs.(27,29) of Marzari and Vanderbilt
                  pos_r(idir) = pos_r(idir) - &
                                wb(nn)*bk(idir, nn, ik)*aimag(log(m_matrix(i, i, nn, ik)))*fac
                else
                  ! Eq.(44) Wang, Yates, Souza and Vanderbilt PRB 74, 195118 (2006)
                  pos_r(idir) = pos_r(idir) + &
                                cmplx_i*wb(nn)*bk(idir, nn, ik)*m_matrix(j, i, nn, ik)*fac
                endif
              end do
            end do
          end do
          write (file_unit, '(2I5,3x,6(E15.8,1x))') j, i, pos_r(:)
        end do
      end do
    end do

    close (file_unit)

    tb_written = .true.

    if (timing_level > 1) call io_stopwatch('hamiltonian: write_tb', 2)

    return

101 call io_error('Error: hamiltonian_write_tb: problem opening file ' &
                  //trim(seedname)//'_tb.dat')

  end subroutine hamiltonian_write_tb

  !============================================!
  subroutine hamiltonian_get_spnr()  ! jmlihm
    !============================================!
    !                                            !
    !!  Calculate the Pauli matrices in the WF basis
    !                                            !
    !============================================!
    use w90_constants, only: dp, pi, cmplx_0
    use w90_parameters, only: num_wann, ndimwin, num_kpts, num_bands, &
      timing_level, have_disentangled, spn_formatted, lsitesymmetry
    use w90_io, only: io_error, io_stopwatch, stdout, seedname, &
      io_file_unit

    implicit none

    complex(kind=dp), allocatable :: spn_o(:, :, :, :), spn_q(:, :, :, :), spn_temp(:, :)
    real(kind=dp)                 :: s_real, s_img
    integer, allocatable          :: num_states(:)
    integer                       :: i, j, ii, jj, m, n, spn_in, ik, is, &
                                     winmin, nb_tmp, nkp_tmp, ierr, s, counter
    character(len=60)             :: header

    if (lsitesymmetry) then
      if (on_root) write (stdout, '(1x,a)') &
        'hamiltonian_get_spnr not implemented for symmetry-adapted WF'
      return
    end if

    if (timing_level > 1) call io_stopwatch('hamiltonian: get_spnr', 1)

    if (.not. allocated(spn_r)) then
      allocate (spn_r(num_wann, num_wann, nrpts, 3))
    else
      return ! been here before
    end if

    allocate (spn_o(num_bands, num_bands, num_kpts, 3))
    allocate (spn_q(num_wann, num_wann, num_kpts, 3))

    allocate (num_states(num_kpts))
    do ik = 1, num_kpts
      if (have_disentangled) then
        num_states(ik) = ndimwin(ik)
      else
        num_states(ik) = num_wann
      endif
    enddo

    ! Read from .spn file the original spin matrices <psi_nk|sigma_i|psi_mk>
    ! (sigma_i = Pauli matrix) between ab initio eigenstates
    !
    spn_in = io_file_unit()
    if (spn_formatted) then
      open (unit=spn_in, file=trim(seedname)//'.spn', form='formatted', &
            status='old', err=109)
      write (stdout, '(/a)', advance='no') &
        ' Reading spin matrices from '//trim(seedname)//'.spn in hamiltonian_get_spnr : '
      read (spn_in, *, err=110, end=110) header
      write (stdout, '(a)') trim(header)
      read (spn_in, *, err=110, end=110) nb_tmp, nkp_tmp
    else
      open (unit=spn_in, file=trim(seedname)//'.spn', form='unformatted', &
            status='old', err=109)
      write (stdout, '(/a)', advance='no') &
        ' Reading spin matrices from '//trim(seedname)//'.spn in hamiltonian_get_spnr : '
      read (spn_in, err=110, end=110) header
      write (stdout, '(a)') trim(header)
      read (spn_in, err=110, end=110) nb_tmp, nkp_tmp
    endif
    if (nb_tmp .ne. num_bands) &
      call io_error(trim(seedname)//'.spn has wrong number of bands')
    if (nkp_tmp .ne. num_kpts) &
      call io_error(trim(seedname)//'.spn has wrong number of k-points')
    if (spn_formatted) then
      do ik = 1, num_kpts
        do m = 1, num_bands
          do n = 1, m
            read (spn_in, *, err=110, end=110) s_real, s_img
            spn_o(n, m, ik, 1) = cmplx(s_real, s_img, dp)
            read (spn_in, *, err=110, end=110) s_real, s_img
            spn_o(n, m, ik, 2) = cmplx(s_real, s_img, dp)
            read (spn_in, *, err=110, end=110) s_real, s_img
            spn_o(n, m, ik, 3) = cmplx(s_real, s_img, dp)
            ! Read upper-triangular part, now build the rest
            spn_o(m, n, ik, 1) = conjg(spn_o(n, m, ik, 1))
            spn_o(m, n, ik, 2) = conjg(spn_o(n, m, ik, 2))
            spn_o(m, n, ik, 3) = conjg(spn_o(n, m, ik, 3))
          end do
        end do
      enddo
    else
      allocate (spn_temp(3, (num_bands*(num_bands + 1))/2), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating spm_temp in hamiltonian_get_spnr')
      do ik = 1, num_kpts
        read (spn_in) ((spn_temp(s, m), s=1, 3), m=1, (num_bands*(num_bands + 1))/2)
        counter = 0
        do m = 1, num_bands
          do n = 1, m
            counter = counter + 1
            spn_o(n, m, ik, 1) = spn_temp(1, counter)
            spn_o(m, n, ik, 1) = conjg(spn_temp(1, counter))
            spn_o(n, m, ik, 2) = spn_temp(2, counter)
            spn_o(m, n, ik, 2) = conjg(spn_temp(2, counter))
            spn_o(n, m, ik, 3) = spn_temp(3, counter)
            spn_o(m, n, ik, 3) = conjg(spn_temp(3, counter))
          end do
        end do
      end do
      deallocate (spn_temp, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating spm_temp in hamiltonian_get_spnr')
    endif

    close (spn_in)

    ! Transform to projected subspace, Wannier gauge
    !
    spn_q(:, :, :, :) = cmplx_0
    do ik = 1, num_kpts
      do is = 1, 3
        call get_gauge_overlap_matrix( &
          ik, num_states(ik), &
          ik, num_states(ik), &
          spn_o(:, :, ik, is), spn_q(:, :, ik, is))
      enddo !is
    enddo !ik

    call fourier_q_to_R(spn_q(:, :, :, 1), spn_r(:, :, :, 1))
    call fourier_q_to_R(spn_q(:, :, :, 2), spn_r(:, :, :, 2))
    call fourier_q_to_R(spn_q(:, :, :, 3), spn_r(:, :, :, 3))

    call hamiltonian_write_spnr()

    if (timing_level > 1 .and. on_root) call io_stopwatch('hamiltonian: get_spnr', 2)
    return

109 call io_error &
      ('Error: Problem opening input file '//trim(seedname)//'.spn')
110 call io_error &
      ('Error: Problem reading input file '//trim(seedname)//'.spn')
  contains

    !==========================================================
    subroutine get_gauge_overlap_matrix(ik_a, ns_a, ik_b, ns_b, S_o, S, H)
      !==========================================================
      !
      ! Wannier-gauge overlap matrix S in the projected subspace
      !
      ! TODO: Update this documentation of this routine and
      ! possibliy give it a better name. The routine has been
      ! generalized multiple times.
      !
      !==========================================================

      use w90_constants, only: dp, cmplx_0
      use w90_parameters, only: num_wann, eigval, u_matrix_opt, u_matrix
      use w90_utility, only: utility_zgemmm
      complex(kind=dp), allocatable :: v_matrix(:, :, :)

      integer, intent(in) :: ik_a, ns_a, ik_b, ns_b

      complex(kind=dp), dimension(:, :), intent(in)            :: S_o
      complex(kind=dp), dimension(:, :), intent(out), optional :: S, H

      integer :: wm_a, wm_b
      integer :: loop_kpt, j, m, i

      allocate (v_matrix(num_bands, num_wann, num_kpts), stat=ierr)
      if (ierr /= 0) &
        call io_error('Error allocating v_matrix in hamiltonian_get_spnr private subroutine')
      ! u_matrix and u_matrix_opt are stored on root only
      if (.not. have_disentangled) then
        v_matrix = u_matrix
      else
        v_matrix = cmplx_0
        do loop_kpt = 1, num_kpts
          do j = 1, num_wann
            do m = 1, ndimwin(loop_kpt)
              do i = 1, num_wann
                v_matrix(m, j, loop_kpt) = v_matrix(m, j, loop_kpt) &
                                           + u_matrix_opt(m, i, loop_kpt)*u_matrix(i, j, loop_kpt)
              enddo
            enddo
          enddo
        enddo
      endif

      call get_win_min(ik_a, wm_a)
      call get_win_min(ik_b, wm_b)

      call utility_zgemmm(v_matrix(1:ns_a, 1:num_wann, ik_a), 'C', &
                          S_o(wm_a:wm_a + ns_a - 1, wm_b:wm_b + ns_b - 1), 'N', &
                          v_matrix(1:ns_b, 1:num_wann, ik_b), 'N', &
                          S, eigval(wm_a:wm_a + ns_a - 1, ik_a), H)
    end subroutine get_gauge_overlap_matrix

    !=========================================================!
    subroutine fourier_q_to_R(op_q, op_R)
      !==========================================================
      !
      !! Fourier transforms Wannier-gauge representation
      !! of a given operator O from q-space to R-space:
      !!
      !! O_ij(q) --> O_ij(R) = (1/N_kpts) sum_q e^{-iqR} O_ij(q)
      !
      !==========================================================

      use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
      use w90_parameters, only: num_kpts, kpt_latt

      implicit none

      ! Arguments
      !
      complex(kind=dp), dimension(:, :, :), intent(in)  :: op_q
      !! Operator in q-space
      complex(kind=dp), dimension(:, :, :), intent(out) :: op_R
      !! Operator in R-space

      integer          :: ir, ik
      real(kind=dp)    :: rdotq
      complex(kind=dp) :: phase_fac

      op_R = cmplx_0
      do ir = 1, nrpts
        do ik = 1, num_kpts
          rdotq = twopi*dot_product(kpt_latt(:, ik), irvec(:, ir))
          phase_fac = exp(-cmplx_i*rdotq)
          op_R(:, :, ir) = op_R(:, :, ir) + phase_fac*op_q(:, :, ik)
        enddo
      enddo
      op_R = op_R/real(num_kpts, dp)

    end subroutine fourier_q_to_R

    !===============================================
    subroutine get_win_min(ik, win_min)
      !===============================================
      !
      !! Find the lower bound (band index) of the
      !! outer energy window at the specified k-point
      !
      !===============================================

      use w90_constants, only: dp
      use w90_parameters, only: num_bands, lwindow, have_disentangled

      implicit none

      ! Arguments
      !
      integer, intent(in)  :: ik
      !! Index of the required k-point
      integer, intent(out) :: win_min
      !! Index of the lower band of the outer energy window

      integer :: j

      if (.not. have_disentangled) then
        win_min = 1
        return
      endif

      do j = 1, num_bands
        if (lwindow(j, ik)) then
          win_min = j
          exit
        end if
      end do

    end subroutine get_win_min
  end subroutine hamiltonian_get_spnr

  !============================================!
  subroutine hamiltonian_write_spnr() ! jmlihm
    !============================================!
    !!  Write the Pauli matrices in the WF basis
    !============================================!

    use w90_io, only: io_error, io_stopwatch, io_file_unit, &
      seedname, io_date
    use w90_parameters, only: num_wann, timing_level, do_write_bin, do_write_text
    integer            :: i, j, irpt, file_unit
    character(len=33) :: header
    character(len=9)  :: cdate, ctime
    integer :: reclen

    if (spnr_written) return

    if (timing_level > 1) call io_stopwatch('hamiltonian: write_spnr', 1)

    ! write the whole matrix with all the indices

    if (do_write_bin) then
      file_unit = io_file_unit()
      inquire (iolength=reclen) spn_r(:, :, :, :)
      open (file_unit, file=trim(seedname)//'_spnr.bin', form='unformatted', &
            status='unknown', access='direct', recl=reclen, err=111)
      write (file_unit, rec=1) spn_r
      close (file_unit)
    end if

    if (do_write_text) then
      file_unit = io_file_unit()
      open (file_unit, file=trim(seedname)//'_spnr.dat', form='formatted', &
            status='unknown', err=112)
      call io_date(cdate, ctime)
      header = 'written on '//cdate//' at '//ctime
      write (file_unit, *) header ! Date and time
      write (file_unit, *) num_wann
      write (file_unit, *) nrpts
      write (file_unit, '(15I5)') (ndegen(i), i=1, nrpts)
      do irpt = 1, nrpts
        do i = 1, num_wann
          do j = 1, num_wann
            write (file_unit, '(5I5,6F12.6)') irvec(:, irpt), j, i, &
              spn_r(j, i, irpt, 1), spn_r(j, i, irpt, 2), spn_r(j, i, irpt, 3)
          end do
        end do
      end do
      close (file_unit)
    end if

    spnr_written = .true.
    if (timing_level > 1) call io_stopwatch('hamiltonian: write_spnr', 2)

    return

111 call io_error('Error: hamiltonian_write_spnr: problem opening file '//trim(seedname)//'_spnr.bin')
112 call io_error('Error: hamiltonian_write_spnr: problem opening file '//trim(seedname)//'_spnr.dat')

  end subroutine hamiltonian_write_spnr

end module w90_hamiltonian
