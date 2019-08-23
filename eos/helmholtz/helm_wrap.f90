      subroutine helm_read_table_c(tab_file_name_c) bind(c)
        use iso_c_binding
        implicit none
        character(kind=c_char), dimension(1), intent(in) :: tab_file_name_c
        character(1024) :: tab_file_name
        integer :: i, n

        n = 1024
        i = 1
        do while((i.le.n-1).and.(tab_file_name_c(i).ne.c_null_char))
          tab_file_name(i:i) = tab_file_name_c(i)
          i = i + 1
        enddo
        tab_file_name(i:i) = c_null_char
        i = i + 1
        do while(i.le.n)
          tab_file_name(i:i) = ' '
          i = i + 1
        enddo

        call read_helm_table(tab_file_name)

        return
      end

      subroutine helm_range_rho_ye_c(rho_ye_min_c, rho_ye_max_c) bind(c)
        use iso_c_binding
        implicit none
        include 'helm_table_storage.dek'
        real (c_double), intent(out) :: rho_ye_min_c, rho_ye_max_c

        rho_ye_min_c = d(1)
        rho_ye_max_c = d(imax)

        return
      end

      subroutine helm_range_temp_c(temp_min_c, temp_max_c) bind(c)
        use iso_c_binding
        implicit none
        include 'helm_table_storage.dek'
        real (c_double), intent(out) :: temp_min_c, temp_max_c

        temp_min_c = t(1)
        temp_max_c = t(jmax)

        return
      end

      integer function offset(rank)
        include 'helm_vector_eos.dek'
        integer rank
        offset = 4*rank
        if(2*offset + 2 .gt. nrowmax) then
            write(*,*) 'Too many threads!'
            stop
        endif
        return
      end

      subroutine helm_range_eps_c(rank, rho_c, abar_c, ye_c, &
                                  eps_min_c, eps_max_c, &
                                  eosfail_c) bind(c)
        use iso_c_binding
        implicit none
        include 'helm_table_storage.dek'
        include 'helm_vector_eos.dek'
        integer (c_int), intent(in)   :: rank
        real (c_double), intent(in)   :: rho_c, abar_c, ye_c
        real (c_double), intent(out)  :: eps_min_c, eps_max_c
        integer (c_int), intent(out) :: eosfail_c
        integer oft, offset

        oft               = offset(rank)
        jlo_eos           = 2*oft + 1
        jhi_eos           = 2*oft + 2
        den_row(jlo_eos)  = rho_c
        den_row(jhi_eos)  = rho_c
        temp_row(jlo_eos) = t(1)
        temp_row(jhi_eos) = t(jmax)
        abar_row(jlo_eos) = abar_c
        abar_row(jhi_eos) = abar_c
        zbar_row(jlo_eos) = ye_c*abar_c
        zbar_row(jhi_eos) = ye_c*abar_c

        call helmeos(jlo_eos, jhi_eos, eosfail)

        eps_min_c = etot_row(jlo_eos)
        eps_max_c = etot_row(jhi_eos)

        eosfail_c = eosfail

        return
      end

      subroutine helm_eos_t_c(rank, rho_c, temp_c, abar_c, ye_c, &
                            press_c, eps_c, entr_c, cs_c, dedt_c, &
                            eosfail_c) bind(c)
        use iso_c_binding
        implicit none
        include 'helm_vector_eos.dek'
        integer (c_int), intent(in)   :: rank
        real (c_double), intent(in)   :: rho_c, temp_c, abar_c, ye_c
        real (c_double), intent(out)  :: press_c, eps_c, entr_c, &
                                         cs_c, dedt_c
        integer (c_int), intent(out) :: eosfail_c
        integer oft, offset

        oft               = offset(rank)
        jlo_eos           = oft + 1
        jhi_eos           = oft + 1
        den_row(jlo_eos)  = rho_c
        temp_row(jlo_eos) = temp_c
        abar_row(jlo_eos) = abar_c
        zbar_row(jlo_eos) = ye_c*abar_c

        call helmeos(jlo_eos, jhi_eos, eosfail)

        press_c = ptot_row(jlo_eos)
        eps_c   = etot_row(jlo_eos)
        entr_c  = stot_row(jlo_eos)
        cs_c    = cs_row(jlo_eos)
        dedt_c  = cv_row(jlo_eos)

        eosfail_c = eosfail

        return
      end

      subroutine helm_eos_e_c(rank, rho_c, eps_c, abar_c, ye_c, &
                             temp_c, press_c, entr_c, cs_c, dedt_c, &
                             eosfail_c) bind(c)
        use iso_c_binding
        implicit none
        include 'helm_vector_eos.dek'
        include 'helm_table_storage.dek'
        integer (c_int), intent(in)    :: rank
        real (c_double), intent(in)    :: rho_c, eps_c, abar_c, ye_c
        real (c_double), intent(inout) :: temp_c
        real (c_double), intent(out)   :: press_c, entr_c, cs_c, dedt_c
        integer (c_int), intent(out)   :: eosfail_c

        integer :: it
        double precision :: delta, dt

        double precision, parameter :: prec = 1.0d-10
        integer, parameter :: maxit = 100
        integer oft, offset

        oft               = offset(rank)
        jlo_eos           = oft + 1
        jhi_eos           = oft + 1
        den_row(jlo_eos)  = rho_c
        temp_row(jlo_eos) = temp_c
        abar_row(jlo_eos) = abar_c
        zbar_row(jlo_eos) = ye_c*abar_c

        call helmeos(jlo_eos, jhi_eos, eosfail)
        if(eosfail) then
          eosfail_c = .true.
          return
        endif

        delta = dabs(etot_row(jlo_eos) - eps_c)
        if(delta .lt. prec*eps_c) then
          temp_c    = temp_row(jlo_eos)
          press_c   = ptot_row(jlo_eos)
          entr_c    = stot_row(jlo_eos)
          cs_c      = cs_row(jlo_eos)
          dedt_c    = cv_row(jlo_eos)
          eosfail_c = .false.
          return
        endif

        it = 1
        do while((it .lt. maxit) .and. (delta .gt. prec*eps_c))
          dt = - (etot_row(jlo_eos) - eps_c)/cv_row(jlo_eos)
          temp_row(jlo_eos) = temp_row(jlo_eos) + dt
          if(temp_row(jlo_eos) .lt. t(jlo_eos)) then
              temp_row(jlo_eos) = t(jlo_eos)
          endif
          if(temp_row(jlo_eos) .gt. t(jmax)) then
              temp_row(jlo_eos) = t(jmax)
          endif

          call helmeos(jlo_eos, jhi_eos, eosfail)
          if(eosfail) then
            eosfail_c = .true.
            return
          endif

          delta = dabs(etot_row(jlo_eos) - eps_c)
          it = it + 1
        enddo

        if(it.eq.maxit) then
          write(6,*) 'temperature root did not converge'
          write(6,*) 'it = ', it, 'maxit = ', maxit
          write(6,*) 'rel err = ', delta/eps_c, 'required = ', prec
          eosfail_c = .true.
          return
        endif

        eosfail_c = .false.

        temp_c  = temp_row(jlo_eos)
        press_c = ptot_row(jlo_eos)
        entr_c  = stot_row(jlo_eos)
        cs_c    = cs_row(jlo_eos)
        dedt_c  = cv_row(jlo_eos)

        return
      end
