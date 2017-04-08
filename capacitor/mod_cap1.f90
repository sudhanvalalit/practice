module mod_cap
  implicit none
  
  integer :: cposition, radius, n_r, n_z, o, dim
  double precision :: eps, V, Phi0, w
  
contains
  subroutine Start()
    implicit none
    
    cposition = 25      ! Platten Lage '-F' 
    radius    = 25      ! Platten Radius '-R' 
    n_r  = 500      ! # Gitterpunkte r '-r' 
    n_z  = 500      ! # Gitterpunkte z '-z' 
    eps  = 0.10e-2! relative Genauigkeit Potential '-e' 
    V    = 1000.0    ! Kondensatorspannung '-V'  
    Phi0 = 0.0    ! Anfangspotential '-U' 
    w    = 1.0    ! Relaxationsparameter '-w' 
    o    = 0    ! Ausgabe-option '-o' 
  end subroutine Start

! =========================== IDX ==================================== 
integer function idx(i,k)
  !use mod_cap
  implicit none
  integer :: i,k
  ! indizierung einer Feldvariable der Dimension (n_z + 1) * (n_r + 1) 
  idx = i + (( n_z + 1 ) * k)
  return
end function idx
! =========================== INIT ==================================== 
subroutine InitPhi(Phi)
  !use mod_cap
  implicit none
  
  integer :: i, k ! Gitterindizes 
  !integer :: idx
  double precision, dimension(0:dim) :: Phi
  
  ! Erst mal Alle Gitterpunkte mit Phi0 vorbesetzen 
  Phi = 0.d0
  do k=0,radius
     ! Potential auf dem Kondensator ist konstant V 
     Phi(idx(cposition,k)) = V
  end do
end subroutine InitPhi
! =========================== Functional ============================== 
double precision function Functional(Phi)
  !use mod_cap
  implicit none

  INCLUDE 'mpif.h'

  integer :: ierr,npe,myid
  
  ! Berechne das "Energie"-Funktional: 
  
!!$                r_mx z_mx    /       2        2 \
!!$             1 /    /        | /dPhi\   /dPhi\  |
!!$    E[Phi] = - | dr | dz   r | |----| + |----|  | 
!!$             2 /    /        | \ dr /   \ dz /  |
!!$              0    0         \                  /
  
  
  integer :: i, k ! Gitterindizes
  !integer :: idx
  integer :: status(MPI_STATUS_SIZE)
  double precision :: sum, term, dPhi, res
  double precision,dimension(0:dim) :: Phi
  
!  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,npe,ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  
  sum = 0.0
 
  do k=myid+1,n_r,npe
!	print*, myid
  !do k=1,n_r
     do i=1,n_z
        dPhi = ( Phi(idx(i,k)) - Phi(idx(i-1,k)) ) ! dPhi / dz 
        term = dPhi * dPhi
        dPhi = ( Phi(idx(i,k)) - Phi(idx(i,k-1)) ) ! dPhi / dr 
        term = term + dPhi * dPhi
        sum = sum + term * dble(k) ! Zylinderkoordinaten: " * r" 
     end do
  end do

  if(myid==0) then
     res = sum
     do i=1,npe-1
        call MPI_RECV(sum,1,MPI_REAL8,i,0, MPI_COMM_WORLD, status, ierr)
        res = res + sum
     end do	
     Functional = res * 0.5 ! eigentlich wird E/h berechnet  
  return
  else
  
!  CALL MPI_REDUCE(sum,res,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  CALL MPI_SEND(sum,1,MPI_REAL8,0,0,MPI_COMM_WORLD,ierr)

  end if

!  CALL MPI_FINALIZE(ierr)
  

end function Functional

! =========================== Er ============================== 
double precision function Er(i,k,Phi)
  !use mod_cap
  implicit none
  
  ! Berechne die r-komponente des elektrischen Feldes bei z_i, r_k 
  
  integer :: i, k
  !integer :: idx
  double precision :: res
  double precision, dimension(0:dim) :: Phi
  
  if (k==0) then
     res = 0.0
  elseif (k<n_r) then
     res = - 0.5 * ( Phi(idx(i,k+1)) - Phi(idx(i,k-1)) )
  else
     res = 0.0
  end if
  Er = res
  return
end function Er

! =========================== Ez ============================== 
double precision function Ez(i,k,Phi)
  !use mod_cap
  implicit none
  
  ! Berechne die z-komponente des elektrischen Feldes bei z_i, r_k 
  
  integer :: i, k
  !integer :: idx
  double precision :: res
  double precision, dimension(0:dim) :: Phi
    
  if (i==0) then
     res = - Phi(idx(i,k))
  elseif (i<n_z) then
     res = - 0.5 * ( Phi(idx(i+1,k)) - Phi(idx(i-1,k)) )
  else
     res = 0.0
  end if
  Ez = res
  return
end function Ez

! =========================== GETUSERPARAM ==================================== 

subroutine GetUserParam()
  !use mod_cap
  implicit none

! Variablen / Konstanten: 
  integer :: argc, i
  character*40 :: argv, endptr
  character*200 :: usage
  character*40 ::error_message
 
  !usage = "capacitor [-F <Plattenposition> -R <Plattenradius> -U <Anfangspotential>  -V <Spannung> -e <relative Genauigkeit>  -r <# r> -z <# z> -w <Relaxationsparameter> -o <Ausgabeoption>]"
  
  !error_message = "# FEHLER(GetuserParam): falsche Option: "

  write(6,*) 'F:'
  read(5,*) cposition
  write(6,*) 'R:'
  read(5,*) radius
  write(6,*) 'U:'
  read(5,*) Phi0
  write(6,*) 'V:'
  read(5,*) V
  write(6,*) 'e:'
  read(5,*) eps
  write(6,*) 'r:'
  read(5,*) n_r
  write(6,*) 'z:'
  read(5,*) n_z
  write(6,*) 'w:'
  read(5,*) w
  

!!$  if (argc>1) then ! falls es ueberhaupt Parameter gibt ... 
!!$     ! Es gibt immer mindestens 1 Parameter: argv[0] = Kommandoname 
!!$     do i=1, argc-1
!!$        ! Parameter sind 2 Charaktere lang und sollten mit '-' anfaengen ... 
!!$        if ( (strlen(argv(i))==2) .and. (argv(i)(0) == '-') ) then 
!!$           if(argv[i][1]=='w') then
!!$              w = strtod( argv(i+1), &endptr) ! String --> double 
!!$			! ... sollte mit Lehrzeichen beendet werden ... 
!!$                        if ( (.not.isspace(*endptr) && (*endptr) != 0) ) { 
!!$			    printf(error_message)
!!$			    printf(usage)
!!$			    printf("\n")
!!$			    exit(1)
!!$                        }
!!$                        break
!!$                    case 'U':
!!$                        Phi0 = strtod( argv[++i], &endptr) 
!!$                        if ( (!isspace(*endptr) && (*endptr) != 0) ) { 
!!$			    printf(error_message)
!!$			    printf(usage)
!!$			    printf("\n")
!!$			    exit(1)
!!$                        }
!!$                        break
!!$                    case 'V':
!!$                        V = strtod( argv[++i], &endptr) 
!!$                        if ( (!isspace(*endptr) && (*endptr) != 0) ) { 
!!$			    printf(error_message)
!!$			    printf(usage)
!!$			    printf("\n")
!!$			    exit(1)
!!$                        }
!!$                        break
!!$                    case 'e':
!!$                        eps = strtod( argv[++i], &endptr)
!!$                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
!!$			    printf(error_message)
!!$			    printf(usage)
!!$			    printf("\n")
!!$			    exit(1)
!!$                        }
!!$                        break
!!$                    case 'F':
!!$			! String --> long --> int ( 10: Dezimalzahl ) 
!!$                        cposition = (int) strtol( argv[++i], &endptr, 10)
!!$                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
!!$			    printf(error_message)
!!$			    printf(usage)
!!$			    printf("\n")
!!$			    exit(1)
!!$                        }
!!$                        break
!!$                    case 'R':
!!$                        radius = (int) strtol( argv[++i], &endptr, 10)
!!$                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
!!$			    printf(error_message)
!!$			    printf(usage)
!!$			    printf("\n")
!!$			    exit(1)
!!$                        }
!!$                        break
!!$                    case 'r':
!!$                        n_r = (int) strtol( argv[++i], &endptr, 10)
!!$                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
!!$			    printf(error_message)
!!$			    printf(usage)
!!$			    printf("\n")
!!$			    exit(1)
!!$                        }
!!$                        break
!!$                    case 'z':
!!$                        n_z = (int) strtol( argv[++i], &endptr, 10)
!!$                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
!!$			    printf(error_message)
!!$			    printf(usage)
!!$			    printf("\n")
!!$			    exit(1)
!!$                        }
!!$                        break
!!$                    case 'o':
!!$                        o = (int) strtol( argv[++i], &endptr, 10)
!!$                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
!!$			    printf(error_message)
!!$			    printf(usage)
!!$			    printf("\n")
!!$			    exit(1)
!!$                        }
!!$                        break
!!$                    default:
!!$			printf(error_message)
!!$			printf(usage)
!!$			printf("\n")
!!$			exit(1)
!!$                }
!!$            } else {
!!$		printf(error_message)
!!$		printf(usage)
!!$		printf("\n")
!!$		exit(1)
!!$            } ! end of: if-then-else 
!!$        end do ! end-of: for  
!!$    end if ! end-of: if 

end subroutine GetUserParam
!!$
!!$Kompilieren:
!!$gcc capacitor.c -l m -o capacitor
!!$Ausfuehren in Verbindung mit gnuplot z.B.:
!!$./capacitor -F 6 -R 6 -r 30 -z 60 -V 500 -o 2 -e 0.5e-2 -U 0 -w 1.7| gnuplot
!!$oder als Bildschirmausgabe:
!!$./capacitor
!!$# Kondensator:
!!$# =========== 
!!$# Gitter: 25 x 25 Punkte
!!$# Kondensatorposition: z=5 Radius: r=5
!!$# Kondensatorspannung: V=1000.000000
!!$# Anfangspotential: Phi0=0.000000
!!$# Relaxationsparameter: Omega=1
!!$# Maximale Abweichung: eps=10
!!$
!!$#Iteration      maximale Abweichung                  Energie
!!$#_________ ________________________ ________________________
!!$         1    3.216388702392579e+02    1.057529752577351e+07
!!$         2    1.851851851851852e+02    8.563781479521886e+06
!!$         3    1.156724215534979e+02    7.606258041563846e+06
!!$         4    8.094986443329907e+01    7.027946880043454e+06
!!$         5    6.926442166505618e+01    6.636039406175652e+06
!!$         6    5.711388006216117e+01    6.355297511441418e+06
!!$         7    4.703664059749542e+01    6.148120609742193e+06
!!$         8    3.824758957342215e+01    5.991826313770101e+06
!!$         9    3.180690612096416e+01    5.871424161616748e+06
!!$        10    2.664496281509363e+01    5.776687317526440e+06
!!$        11    2.196427346622710e+01    5.700568349378514e+06
!!$        12    1.817195055240620e+01    5.638180863156057e+06
!!$        13    1.612524443203483e+01    5.586106697652645e+06
!!$        14    1.438133266063471e+01    5.541926048868635e+06
!!$        15    1.289062975161625e+01    5.503901260247231e+06
!!$        16    1.161085742222792e+01    5.470764679563234e+06
!!$        17    1.063083031516152e+01    5.441576057238697e+06
!!$        18    9.797194535102562e+00    5.415626137300598e+06
!!$


end module mod_cap
