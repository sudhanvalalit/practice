program capacitor
  use mod_cap

  implicit none
  include "mpif.h"
!!$   Datei:  capacitor.c
!!$   Author: BCM
!!$   Datum:  20.06.2009 
!!$   Datum:  25.06.2012 
!!$
!!$   Berechnung des elektrostatischen Potentials Phi(z,r) eines Plattenkondensators
!!$   in Zylinderkoordinaten mit scheibenfoermigen Platten im Abstand [2*cposition] 
!!$   mit Radius [radius] in einem geerdeten Zylinder der Laenge [2*n_z] und Radius [n_r]

! Vorwaertsdeklaration der Funktionen:

!!$void GetUserParam( int, char *[] ) ! Eingabe ueber die Kommandozeile: 
!!$void InitPhi( double[] )           ! Initialisierung des Potentials 
!!$double Functional( double[] )      ! Funktional 
!!$double Er( int, int, double[] )    ! E-Feld in r-Richtung 
!!$double Ez( int, int, double[] )    ! E-Feld in z-Richtung 
!!$long idx( int, int )               ! Indizierung des Potentialfeldes 

! Globale Variablen: 

  integer :: iter, i, k, ierr, myid, size
  double precision :: delta, deltamx, Phi_old, U, resulting
  double precision, allocatable,dimension(:) :: Phi

  !integer :: idx
  !double precision :: Er,Ez,Functional
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
!  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, size)

!  if (size .ne. 4) MPI_ABORT(MPI_COMM_WORLD, 1)

  call Start()
  
if(myid==0) then
  !call GetUserParam() ! Parameter ueber die Kommandozeile 
  
  eps = eps*V ! Absolute Genauigkeit des Potentials 
 
  dim = ( n_r + 1 ) * ( n_z + 1 ) ! # der Feldvariable Phi 
  allocate(Phi(0:dim)) ! Speicherplatzreservierung Phi
  call InitPhi(Phi) ! Initialisiere Phi und setze Randwerte 
  
!  Phi(0) = 1.d0
  
  iter = 0 ! Iterationszaehler 
  
  ! Ausgabe der Parameter: 
  write(6,*) "# Kondensator:"
  write(6,*) "# =========== "
  write(6,*) "# Gitter: ",n_z," x ",n_r," Punkte"
  write(6,*) "# Kondensatorposition: z=",cposition," Radius: r=", radius
  write(6,*) "# Kondensatorspannung: V=", V
  write(6,*) "# Anfangspotential: Phi0=", phi0
  write(6,*) "# Relaxationsparameter: Omega=", w
  write(6,*) "# Maximale Abweichung: eps=",eps
  write(6,*)
  
  
  if (o==0) then ! Ausgabe fuer maximale Abweichung ( und Energie ): 
     !                1         2         3         4         5         6 
     !       123456789012345678901234567890123456789012345678901234567890 
     write(6,*) "#Iteration      maximale Abweichung                  Energie"
     write(6,*) "#_________ ________________________ ________________________"
  end if
  
  if (o==2) then ! Ausgabe direkt nach gnuplot 
     write(6,*) "#set hidden" 
     write(6,*) "#set contour base"  ! Hoehelinien 
     write(6,*) "#set cntrpar levels incr 0,",V/10,",",V
     write(6,*) "set view map"  ! Hoehelinien 
     write(6,*) "set xr[",0,":",n_z,"]\n"  ! Plotbereich 
     write(6,*) "set yr[",0,":",n_r,"]\n" 
     write(6,*) "set zr[",0,":",V,"]"
     write(6,*) "set xla 'z'"  ! Beschriftung 
     write(6,*) "set yla 'r'" 
     write(6,*) "set zla 'Phi'" 
     write(6,*) "set title 'Kondensator: Omega=",w,"   Iteration=",iter,"'" 
     write(6,*) "set label 1 'Iteration:",iter,", E:",Functional(Phi),"' at ",0,",",n_r/2,",",V," left"
     write(6,*) "splot '-' u 1:2:3 t 'Phi' w pm3d lt 1" ! Data von stdin 
     do k=0,n_r ! Schleife ueber Gitter 
        do i=0, n_z
           write(6,*) i, k, Phi(idx(i,k)), Ez(i,k,Phi), Er(i,k,Phi)
        end do
        write(6,*) 
     end do
     write(6,*) "e"  ! 'e' beendet Data 
     write(6,*) "pause mouse"  ! weiter nach Mausklick 
  end if
 
!     do k=0, n_r
!        do i=0, n_z
!           write(6,*) i, k, Phi(idx(i,k)), Ez(i,k,Phi), Er(i,k,Phi)
!        end do
!  write(6,*) phi(idx(5,k))
!        write(6,*)
!     end do
 
  deltamx = 100000000.d0
end if
  do while(deltamx > eps) ! Relaxationsschritte: 
!	if (myid==0) then     
	deltamx = 0.0
     ! Schleife ueber Gitterpunkte in radialer Richtung: 
     do k=0, n_r-1
        ! Schleife ueber Gitterpunkte in axialer Richtung: 
        do i=1, n_z-1
           ! Potential auf dem Kondensator ist konstant V 
           if (.not.( (i==cposition) .and. (k<=radius) )) then 
           if (k==0) then ! Zentrale Achse ist speziell: 
              U = ( Phi(idx(i+1,k)) + Phi(idx(i-1,k)) + 4 * Phi(idx(i,k+1)) ) / 6.0
           else  ! Punkte im Innern sind standard: 
              U = (   Phi(idx(i,k+1)) + Phi(idx(i,k-1)) + Phi(idx(i+1,k)) + Phi(idx(i-1,k)) ) / 4.0 &
                   & + ( Phi(idx(i,k+1)) - Phi(idx(i,k-1)) ) / ( 8.0 * k )
           end if
           Phi_old = Phi(idx(i,k))
           Phi(idx(i,k)) = Phi_old + w * ( U - Phi_old )
           ! Berechne maximale Abweichung: 
           delta = abs(Phi(idx(i,k))-Phi_old)
           if ( delta > deltamx ) then
              deltamx = delta
           end if
	end if
        end do
     end do
     
     iter = iter + 1 ! Zaehle Iterationen 
!	end if
!     print*, "iter=", iter
     if (o==0) then ! Ausgabe fuer maximale Abweichung ( und Energie ): 
	resulting = Functional(Phi)
	if(myid==0) then
        	write(6,*) iter, deltamx, resulting
	end if
     end if
     
     if (o==2) then ! o=2: Ausgabe nach gnuplot fuer Animation 
        write(6,*) "set title 'Kondensator: Omega=%g   Iteration=%d'", w,iter
        write(6,*) "set label 1 'Iteration:%i, E:%g' at %i,%i,%f left", iter, Functional(Phi), 0, n_r/2, V
        write(6,*) "splot '-' u 1:2:3 w pm3d lt 1" ! Data von stdin 
        do k=0, n_r
           do i=0, n_z
              write(6,*) i, k, Phi(idx(i,k)), Ez(i,k,Phi), Er(i,k,Phi)
              write(6,*)
           end do
           write(6,*) "e" ! 'e' beendet Data 
        end do
     end if
  end do ! Ende Relaxationsschritte 
  
  
  if (o==2) then ! Ausgabe direkt nach gnuplot 
     write(6,*) "pause mouse"  ! weiter nach Mausklick 
  end if
  
  ! o=1: Ausgabe von Phi, Ez, Er: 
  if (o==1) then
     !                1         2         3         4         5 
     !       12345678901234567890123456789012345678901234567890 
     write(6,*) "#        z         r                      Phi" 
     write(6,*) "                       Ez                       Er"
     write(6,*) "#_________ _________ ________________________" 
     write(6,*) " ________________________ ________________________" 
 end if
     
!     do k=0, n_r
!        do i=0, n_z
!           write(6,*) i, k, Phi(idx(i,k)), Ez(i,k,Phi), Er(i,k,Phi)
!        end do
!  write(6,*) phi(idx(5,k))
!        write(6,*)
!     end do
!  end if


  CALL MPI_FINALIZE(ierr)
end program capacitor
