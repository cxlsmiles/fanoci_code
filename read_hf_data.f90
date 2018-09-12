module read_hf_data
    use globals
    use strings
    implicit none
    contains


    subroutine read_n_mos()
        integer          :: io
        character(len=5) :: temp1
        double precision :: temp2

        open(unit=118, file='HForbenergy.txt')
        do
          read(118,*,iostat=io)n_mo,temp1,temp2
          if (io/=0) exit
        end do
        rewind(118)
        close(118)

        write(*,*)'Total number of orbitals from HForbenergy.txt = ', n_mo

    end subroutine read_n_mos


    subroutine read_input()
        integer :: funit
        integer :: stat
        integer :: k
        integer, parameter :: StrMax=40, Nmax = 80
        character(len=1) :: delims
        character(len=StrMax), dimension(Nmax) :: h_args
        character(len=StrMax), dimension(Nmax) :: p_args


        namelist /prop/ prop_type
        namelist /init/ hs_in, ps_in
        namelist /final/ h_f_min, h_f_max, p_f_min, p_f_max
        namelist /control/ n_occ, prnt, gam, moint, eig_thresh

        funit = 139

        open (funit, file='input.txt',status='old')
        read (funit,nml=prop)
        read (funit,nml=init)
        read (funit,nml=final)
        read (funit,nml=control)
        close(funit)

        delims = ' '
        call parse(hs_in, delims, h_args, n_hs_in)
        call parse(ps_in, delims, p_args, n_ps_in)
    
        allocate(holes_arr(n_hs_in), particles_arr(n_ps_in))
        
        do k = 1, n_hs_in
           read(h_args(k),*,iostat=stat)holes_arr(k)
        end do

        do k = 1, n_ps_in
           read(p_args(k),*,iostat=stat)particles_arr(k)
        end do
        
    end subroutine



    subroutine read_mo_energies ()
        integer*8 :: i, temp
        character(3) :: sym


        open(unit=117, file='HForbenergy.txt')
        do i = 1, n_mo
            read(117,*)temp,sym,emo(i)
        end do
        close(117)

    end subroutine read_mo_energies




    subroutine read_2e_integrals (n_, int2e_)
        integer*8, intent(in) :: n_
        integer*8 :: i
        integer*8 :: j
        double precision :: temp
        double precision, dimension(n_) :: int2e_
        int2e(:) = 0.d0

        open(unit=119, file='moint.txt') !, form='unformatted')
        do i = 1, n_
            read(119,*,end=100)j,temp
            int2e_(j) = temp
        end do
        100 int2e_(j) = temp
        close(119)

    end subroutine read_2e_integrals



end module read_hf_data
