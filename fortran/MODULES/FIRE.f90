module FIRE
    implicit none

	real FIRE_finc     
	real FIRE_fdec     
	real FIRE_acoef
	real FIRE_Ftot
	real FIRE_acoef0   
	real FIRE_falpha   
	integer FIRE_Nmin  ! currently not used
	real FIRE_mass     
	real FIRE_Fclamp 
	real FIRE_dt  
	real FIRE_dtmax
	real FIRE_dtmin

    contains

subroutine init_FIRE( )
	!use optimization
	!use MD
	!use fragments
	use configuration
	use loops, only : dt
	logical file_exists
	inquire(file="FIRE.optional", exist=file_exists ) 
	! set FIRE parameters
	FIRE_dt = dt
	if ( file_exists ) then
		write (*,*) " loading parameters from FIRE.optional "
		open (unit= 33, file="FIRE.optional", status='unknown')
		read (33,*)  FIRE_finc 
		read (33,*)  FIRE_fdec 
		read (33,*)  FIRE_acoef0  
		read (33,*)  FIRE_falpha 
		read (33,*)  FIRE_Nmin 
		read (33,*)  FIRE_mass  
		read (33,*)  FIRE_Fclamp  
		read (33,*)  FIRE_dtmax
		read (33,*)  FIRE_dtmin
		close (33)   
	else
		write (*,*) " No FIRE.optional => setting default parameters "
		FIRE_finc     = 1.1D0
		FIRE_fdec     = 0.5D0
		FIRE_acoef0   = 0.1D0
		FIRE_falpha   = 0.99D0
		FIRE_Nmin     = 5        ! currently not used
		FIRE_mass     = 4.0D0
		FIRE_Fclamp   = 10.0D0   ! too big force
		FIRE_dtmax    = FIRE_dt
		FIRE_dtmin    = FIRE_dt*0.1
	end if
	! set FIRE varaibles
	FIRE_dt       = FIRE_dtmax * 0.25
	FIRE_acoef    = FIRE_acoef0

	vatom(:,:) = 0

	!if ( .not. allocated( fragxyz ) ) then
	!	allocate( fragxyz(3,natoms) ) 
	!	fragxyz( :,:) = 0
	!end if

end subroutine init_FIRE

! ======================================
! =============== move_ions_FIRE( )
! ======================================

subroutine move_ions_FIRE( istep, iwrtxyz )
	!use outputs, only: iwrtxyz
	!use optimization
	use configuration
	use forces
	!use fragments
	use energy, only: deltaFmax

 implicit none
! == arguments 
	integer, intent(in) :: istep
	integer, intent(in) :: iwrtxyz
! == Local variables
	integer iatom, k
	real ff,vv,vf, cF, cV, dtv
 	character (200 ) xyz_headline
! == Procedure
	
! projection of force to velocity <v|f>
	ff = 0.D0
	vv = 0.D0
	vf = 0.D0
    deltaFmax = 0.0d0
	do iatom = 1, natoms
		do k=1,3
			!if ( fragxyz(k,iatom) .eq. 0 ) then 
				deltaFmax = max(deltaFmax, abs(ftot(k,iatom)))
				ff = ff + ftot(k,iatom)**2
				vv = vv + vatom(k,iatom)**2
				vf = vf + vatom(k,iatom) * ftot(k,iatom)
			!end if ! fragxyz(k,iatom)
		end do ! k
	end do ! iatom
	FIRE_Ftot = sqrt( ff )
	

! FIRE update depending of projection <v|f>
	if ( vf .lt. 0 ) then
		vatom(:,:)   = 0 
		!FIRE_dt      = FIRE_dt * FIRE_fdec
		FIRE_dt      = max( FIRE_dt * FIRE_fdec, FIRE_dtmin )
		FIRE_acoef   = FIRE_acoef0
    else
		cF           =     FIRE_acoef * sqrt(vv/ff)
		cV           = 1 - FIRE_acoef
		do iatom = 1, natoms
			do k=1,3
				!if ( fragxyz(k,iatom) .eq. 0 ) then 
					vatom(k,iatom) = cV * vatom(k,iatom)    +    cF *  ftot(k,iatom)
				!end if ! fragxyz(k,iatom)
			end do ! k
		end do ! iatom
		FIRE_dt     = min( FIRE_dt * FIRE_finc, FIRE_dtmax ) 
		FIRE_acoef  = FIRE_acoef   * FIRE_falpha
    end if

!  normal MD step using leap-frog 
	dtv  = FIRE_dt / FIRE_mass 
	if( deltaFmax .gt. FIRE_Fclamp ) then  ! if force is too big
		dtv        = dtv * FIRE_Fclamp / deltaFmax
		vatom(:,:) = vatom(:,:) * 0.5 
	end if
	do iatom = 1, natoms
		do k=1,3
			!if ( fragxyz(k,iatom) .eq. 0 ) then 
				vatom(k,iatom) = vatom(k,iatom)   +  dtv     * ftot (k,iatom)
				ratom(k,iatom) = ratom(k,iatom)   +  FIRE_dt * vatom(k,iatom)
			!end if
		enddo
	enddo

	write ( xyz_headline, '(A, i6, 6f16.8)' ) " #### FIRE: i,Fmax,|F|,v,<v|f>,dt: ", istep, deltaFmax, FIRE_Ftot, vv, sqrt(vv), vf, FIRE_dt  
	write ( *, '(A)' )  xyz_headline
	!write(*,*) "DEBGU move_ions_FIRE() : iwrtxyz ", iwrtxyz
	if ( iwrtxyz .eq. 1 ) call write_to_xyz( xyz_headline, istep )

end subroutine move_ions_FIRE


! ======================================
!        OUTPUTS ROUTINES
! ======================================

! =============== write_bas(  )
subroutine write_bas(  )
	use configuration, only: natoms,ratom
	use interactions, only: imass
	use charges, only: nzx
 implicit none
! == variables 
	 integer iatom, in1
! == body
	open (unit = 17, file = 'answer.bas', status = 'unknown')
	write (17,*) natoms
	do iatom = 1, natoms
		in1 = imass(iatom)
		write (17, '(i5, 3f16.8)'  ) nzx(in1), ratom(:,iatom)
	enddo
	close (unit = 17)
end subroutine write_bas


! ===============  write_to_xyz(  )
subroutine write_to_xyz( headline, istep )
	use configuration, only: natoms,ratom,symbol
	use interactions, only: imass
 implicit none
! == arguments 
	character (*), intent(in) :: headline
	integer      , intent(in) :: istep
! == variables 
	 integer iatom, in1
! == body
	if ( istep .eq. 1 ) then
		open (unit = 18, file = 'answer.xyz', status = 'unknown')
	else
		open (unit = 18, file = 'answer.xyz', status = 'unknown', position = 'append')
	endif
	write (18,*) natoms
	write (18,'(A)') headline
	do iatom = 1, natoms
		in1 = imass(iatom)
		write (18, '(A, 3f16.8)' ) symbol(iatom), ratom(:,iatom)
	enddo
	close (unit = 18)
end subroutine write_to_xyz

end module FIRE









