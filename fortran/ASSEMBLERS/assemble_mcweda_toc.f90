


! ============================= FROM    assemble_mcweda ()

if ((itheory .eq. 1 .or. itheory .eq. 2) .and. Kscf .eq. 1) then
   call get_ewald ! (nprocs, my_proc, kforce, icluster, itheory, iordern)
end if

if(itheory_xc .eq. 2 ) then
   call assemble_olsxc_1c   ! no-interpolation
endif

if (Kscf .eq. 1) then
   call assemble_sVNL()     ! doscentrosPP()
   call assemble_2c()       ! doscentros() 
   call assemble_2c_PP()    ! no-interpolation
end if

if ( (itheory_xc .eq. 1) .or. (itheory_xc .eq. 2) ) then
   !write (*,*) ' Assemble SN-xc exchange-correlation interactions. '
   if (itheory .eq. 1) then
      call average_ca_rho()  !  doscentros(), doscentrosS(),trescentros(),trescentrosS()  ... really complicated function
   else
      call average_rho()     !  doscentros(), doscentrosS(),trescentros(),trescentrosS()  ... really complicated function
   endif
end if

if      (itheory_xc .eq. 1) then
   call assemble_snxc_on ()   ! no-interpolation
   call assemble_snxc_off()   ! no-interpolation
else if (itheory_xc .eq. 2 ) then
   call assemble_olsxc_on()  ! no-interpolation
   call assemble_olsxc_off() ! doscentros()
end if ! if (itheory_xc = 2)



if (itheory .eq. 1) then
   call assemble_ca_2c()  ! doscentros()
endif
if (Kscf .eq. 1) then
   call assemble_3c()     ! trescentros()
   call assemble_3c_PP()  ! no-interpolation
end if
if (itheory .eq. 1) then
   call assemble_ca_3c()  ! trescentros()
   call assemble_lr   ()  ! no-interpolation
else




!! ============ FOR ( itheory .eq. 1) and  (itheory_xc .eq. 2 )  =============

    if ( Kscf .eq. 1) then
        call get_ewald         
        call assemble_sVNL()     ! doscentrosPP()
        call assemble_2c()       ! doscentros() 
        call assemble_2c_PP()    ! no-interpolation
    end if
    call average_ca_rho()        ! doscentros(), doscentrosS(),trescentros(),trescentrosS()  ... really complicated function
    call assemble_olsxc_on()     ! no-interpolation
    call assemble_olsxc_off()    ! doscentros()
    call assemble_ca_2c()        ! doscentros()
    if (Kscf .eq. 1) then
        call assemble_3c()       ! trescentros()
        call assemble_3c_PP()    ! no-interpolation
    end if
    call assemble_ca_3c()        ! trescentros()
    call assemble_lr   ()        ! no-interpolation












