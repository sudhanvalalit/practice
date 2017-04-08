/* Datei:  capacitor.c
   Author: BCM
   Datum:  20.06.2009 
   Datum:  25.06.2012 */

/* Berechnung des elektrostatischen Potentials Phi(z,r) eines Plattenkondensators
   in Zylinderkoordinaten mit scheibenfoermigen Platten im Abstand [2*cposition] 
   mit Radius [radius] in einem geerdeten Zylinder der Laenge [2*n_z] und Radius [n_r]
*/

/* Zunaechst Definitionen der Standardbibliotheken einbinden mit #include */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>       /* mathematische Funktionen */ 

/* Vorwaertsdeklaration der Funktionen: */

void GetUserParam( int, char *[] ); /* Eingabe ueber die Kommandozeile: */
void InitPhi( double[] );           /* Initialisierung des Potentials */
double Functional( double[] );      /* Funktional */
double Er( int, int, double[] );    /* E-Feld in r-Richtung */
double Ez( int, int, double[] );    /* E-Feld in z-Richtung */
long idx( int, int );               /* Indizierung des Potentialfeldes */

/* Globale Variablen: */

int cposition =    5;      /* Platten Lage '-F' */
int radius    =    5;      /* Platten Radius '-R' */
int n_r       =   25;      /* # Gitterpunkte r '-r' */
int n_z       =   25;      /* # Gitterpunkte z '-z' */
double eps    =    1.0e-2; /* relative Genauigkeit Potential '-e' */
double V      = 1000.0;    /* Kondensatorspannung '-V' */ 
double Phi0   =    0.0;    /* Anfangspotential '-U' */
double w      =    1.0;    /* Relaxationsparameter '-w' */
int o         =      0;    /* Ausgabe-option '-o' */

/* ======================= MAIN ===========================*/

int main( int argc, char *argv[] ) {
 
    int iter;              /* Iterationsschrittzaehler */
    long dim;              /* zur Speichplatzreservierung Phi */
    double delta, deltamx; /* (maximale) Abweichung in einer Iteration */
    double *Phi;           /* Potentialfeld */ 
    double Phi_old;        /* alter Potentialwert */
    double U;              /* neuer Potentialwert */
    int i, k;              /* Indizes Gitterpunkte */

    GetUserParam( argc, argv ); /* Parameter ueber die Kommandozeile */

    eps *= V; /* Absolute Genauigkeit des Potentials */

    dim = ( (long) n_r + 1 ) * ( (long) n_z + 1 ); /* # der Feldvariable Phi */
    Phi = calloc( dim, sizeof(double) ); /* Speicherplatzreservierung Phi */
    InitPhi( Phi ); /* Initialisiere Phi und setze Randwerte */

    iter = 0; /* Iterationszaehler */

    /* Ausgabe der Parameter: */
    printf( "# Kondensator:\n" );
    printf( "# =========== \n" );
    printf( "# Gitter: %i x %i Punkte\n", n_z, n_r );
    printf( "# Kondensatorposition: z=%i; Radius: r=%i\n", cposition, radius );
    printf( "# Kondensatorspannung: V=%f\n", V );
    printf( "# Anfangspotential: Phi0=%f\n", Phi0 );
    printf( "# Relaxationsparameter: Omega=%g\n", w );
    printf( "# Maximale Abweichung: eps=%g\n\n", eps );

    
    if (o==0) { /* Ausgabe fuer maximale Abweichung ( und Energie ): */
	/*                1         2         3         4         5         6 */
	/*       123456789012345678901234567890123456789012345678901234567890 */
	printf( "#Iteration      maximale Abweichung                  Energie\n" );
	printf( "#_________ ________________________ ________________________\n" );
    }
    
    if (o==2) { /* Ausgabe direkt nach gnuplot */
	printf( "#set hidden\n" );
	printf( "#set contour base\n" ); /* Hoehelinien */
	printf( "#set cntrpar levels incr 0,%f,%f\n", V/10, V );
	printf( "set view map \n" ); /* Hoehelinien */
	printf( "set xr[%i:%i]\n", 0, n_z ); /* Plotbereich */
	printf( "set yr[%i:%i]\n", 0, n_r );
	printf( "set zr[%i:%f]\n", 0, V );
	printf( "set xla 'z'\n" ); /* Beschriftung */
	printf( "set yla 'r'\n" );
	printf( "set zla 'Phi'\n" );
	printf( "set title 'Kondensator: Omega=%g   Iteration=%d'\n", w,iter );
	printf( "set label 1 'Iteration:%i, E:%g' at %i,%i,%f left\n", 
		iter, Functional(Phi), 0, n_r/2, V );
	printf( "splot '-' u 1:2:3 t 'Phi' w pm3d lt 1\n"); /* Data von stdin */
	for (k=0; k<=n_r; k++) { /* Schleife ueber Gitter */
	    for (i=0; i<=n_z; i++) {
		printf( "%10i%10i%25.15e%25.15e%25.15e\n", 
			i, k, Phi[idx(i,k)], Ez(i,k,Phi), Er(i,k,Phi) );
	    }
	    printf( "\n" );
	}
	printf( "e\n" ); /* 'e' beendet Data */
	fflush(stdout); /* flush */
	printf( "pause mouse\n" ); /* weiter nach Mausklick */
	fflush(stdout);
    }

    do { /* Relaxationsschritte: */
	deltamx = 0.0;
	/* Schleife ueber Gitterpunkte in radialer Richtung: */
	for (k=0; k<n_r; k++) {
	    /* Schleife ueber Gitterpunkte in axialer Richtung: */
	    for (i=1; i<n_z; i++) {
		/* Potential auf dem Kondensator ist konstant V */
		if ( (i==cposition) && (k<=radius) ) continue; 
		if (k==0) { /* Zentrale Achse ist speziell: */
		    U = ( Phi[idx(i+1,k)] + Phi[idx(i-1,k)] 
			  + 4 * Phi[idx(i,k+1)] ) / 6.0;
		} else { /* Punkte im Innern sind standard: */
		    U = (   Phi[idx(i,k+1)] + Phi[idx(i,k-1)] 
			  + Phi[idx(i+1,k)] + Phi[idx(i-1,k)] ) / 4.0
			+ 
			( Phi[idx(i,k+1)] - Phi[idx(i,k-1)] ) / ( 8.0 * k );
		}
		Phi_old = Phi[idx(i,k)];
		Phi[idx(i,k)] = Phi_old + w * ( U - Phi_old );
		/* Berechne maximale Abweichung: */
		delta = fabs(Phi[idx(i,k)]-Phi_old);
		if ( delta > deltamx ) {
		    deltamx = delta;
		}
	    }
	}		   
	
	iter++; /* Zaehle Iterationen */

	if (o==0) { /* Ausgabe fuer maximale Abweichung ( und Energie ): */
	    printf( "%10i%25.15e%25.15e\n", iter, deltamx, Functional(Phi) );
	}

	if (o==2) { /* o=2: Ausgabe nach gnuplot fuer Animation */
            printf( "set title 'Kondensator: Omega=%g   Iteration=%d'\n", w,iter );
	    printf( "set label 1 'Iteration:%i, E:%g' at %i,%i,%f left\n", 
		    iter, Functional(Phi), 0, n_r/2, V );
	    printf( "splot '-' u 1:2:3 w pm3d lt 1\n"); /* Data von stdin */
	    for (k=0; k<=n_r; k++) {
		for (i=0; i<=n_z; i++) {
		    printf( "%10i%10i%25.15e%25.15e%25.15e\n", 
			    i, k, Phi[idx(i,k)], Ez(i,k,Phi), Er(i,k,Phi) );
		}
		printf( "\n" );
	    }
	    printf( "e\n" ); /* 'e' beendet Data */
	    fflush(stdout); /* flush */
	}

    } while (deltamx > eps); /* Ende Relaxationsschritte */


    if (o==2) { /* Ausgabe direkt nach gnuplot */
	printf( "pause mouse\n" ); /* weiter nach Mausklick */
	fflush(stdout);
    }

    /* o=1: Ausgabe von Phi, Ez, Er: */
    if (o==1) {
	/*                1         2         3         4         5 */
	/*       12345678901234567890123456789012345678901234567890 */
	printf( "#        z         r                      Phi" );
	printf( "                       Ez                       Er\n");
	printf( "#_________ _________ ________________________" );
	printf( " ________________________ ________________________\n" );

	for (k=0; k<=n_r; k++) {
	    for (i=0; i<=n_z; i++) {
		printf( "%10i%10i%25.15e%25.15e%25.15e\n", 
			i, k, Phi[idx(i,k)], Ez(i,k,Phi), Er(i,k,Phi) );
	    }
	    printf( "\n" );
	}
    }

    return 0;
}
/* =========================== IDX ==================================== */
long idx( int i, int k ) {
    /* indizierung einer Feldvariable der Dimension (n_z + 1) * (n_r + 1) */
    return (long) i + ( (long) n_z + 1 ) * (long) k;
}
/* =========================== INIT ==================================== */
void InitPhi( double Phi[] ) {

    int i, k; /* Gitterindizes */

    /* Erst mal Alle Gitterpunkte mit Phi0 vorbesetzen */
    for (k=0; k<=n_r; k++) {
	for (i=0; i<=n_z; i++) {
	    Phi[idx(i,k)] = Phi0;
	}
    }
    for (k=0; k<=n_r; k++) {
	/* Potential in der Mitte (z=0) ist immer Null */
	Phi[idx(0,k)] = 0.0; 
    }
    for (k=0; k<=radius; k++) {
	/* Potential auf dem Kondensator ist konstant V */
	Phi[idx(cposition,k)] = V;
    }
    for (i=0; i<=n_z; i++) {
	/* Zylindermantel ist geerdet */
	Phi[idx(i,n_r)] = 0;
    }
    for (k=0; k<=n_r; k++) {
	/* Zylinderkappe ist geerdet */
	Phi[idx(n_z,k)] = 0;
    }
}
/* =========================== Functional ============================== */
double Functional( double Phi[] ) {

    /* Berechne das "Energie"-Funktional: 

                r_mx z_mx    /       2        2 \
             1 /    /        | /dPhi\   /dPhi\  |
    E[Phi] = - | dr | dz   r | |----| + |----|  | 
             2 /    /        | \ dr /   \ dz /  |
              0    0         \                  /
    */

    int i, k; /* Gitterindizes */
    double sum, term, dPhi; 

    sum = 0.0;
    for (k=1; k<=n_r; k++) {
	for (i=1; i<=n_z; i++) {
	    dPhi = ( Phi[idx(i,k)] - Phi[idx(i-1,k)] ); /* dPhi / dz */
	    term = dPhi * dPhi;
	    dPhi = ( Phi[idx(i,k)] - Phi[idx(i,k-1)] ); /* dPhi / dr */
	    term += dPhi * dPhi;
	    sum += term * (double) k; /* Zylinderkoordinaten: " * r" */
	}
    }
    sum *= 0.5; /* eigentlich wird E/h berechnet */ 

    return sum;
}

/* =========================== Er ============================== */
double Er( int i, int k, double Phi[] ) {

    /* Berechne die r-komponente des elektrischen Feldes bei z_i, r_k */

    double res;
    
    if (k==0) {
	res = 0.0;
    } else if (k<n_r) {
	res = - 0.5 * ( Phi[idx(i,k+1)] - Phi[idx(i,k-1)] );
    } else {
	res = 0.0;
    }

    return res;
}

/* =========================== Ez ============================== */
double Ez( int i, int k, double Phi[] ) {

    /* Berechne die z-komponente des elektrischen Feldes bei z_i, r_k */

    double res;
    
    if (i==0) {
	res = - Phi[idx(i,k)];
    } else if (i<n_z) {
	res = - 0.5 * ( Phi[idx(i+1,k)] - Phi[idx(i-1,k)] );
    } else {
	res = 0.0;
    }

    return res;
}

/* =========================== GETUSERPARAM ==================================== */

void GetUserParam( int argc, char *argv[] ){

/* Variablen / Konstanten: */
    int i;
    char* endptr;
    const char* usage = 
        "capacitor [-F <Plattenposition> -R <Plattenradius> -U <Anfangspotential>\
 -V <Spannung> -e <relative Genauigkeit>\
 -r <# r> -z <# z> -w <Relaxationsparameter> -o <Ausgabeoption>]";
    const char* error_message =
        "# FEHLER(GetuserParam): falsche Option: ";

    if (argc>1) { /* falls es ueberhaupt Parameter gibt ... */
	/* Es gibt immer mindestens 1 Parameter: argv[0] = Kommandoname */
        for (i=1; i<argc; i++){
            /* Parameter sind 2 Charaktere lang und sollten mit '-' anfaengen ... */
            if ( (strlen(argv[i])==2) && (argv[i][0] == '-') ) { 
                switch (argv[i][1]) { 
                    case 'w':
                        w = strtod( argv[++i], &endptr); /* String --> double */
			/* ... sollte mit Lehrzeichen beendet werden ... */
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) { 
			    printf(error_message);
			    printf(usage);
			    printf("\n");
			    exit(1);
                        }
                        break;
                    case 'U':
                        Phi0 = strtod( argv[++i], &endptr); 
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) { 
			    printf(error_message);
			    printf(usage);
			    printf("\n");
			    exit(1);
                        }
                        break;
                    case 'V':
                        V = strtod( argv[++i], &endptr); 
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) { 
			    printf(error_message);
			    printf(usage);
			    printf("\n");
			    exit(1);
                        }
                        break;
                    case 'e':
                        eps = strtod( argv[++i], &endptr);
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
			    printf(error_message);
			    printf(usage);
			    printf("\n");
			    exit(1);
                        }
                        break;
                    case 'F':
			/* String --> long --> int ( 10: Dezimalzahl ) */
                        cposition = (int) strtol( argv[++i], &endptr, 10);
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
			    printf(error_message);
			    printf(usage);
			    printf("\n");
			    exit(1);
                        }
                        break;
                    case 'R':
                        radius = (int) strtol( argv[++i], &endptr, 10);
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
			    printf(error_message);
			    printf(usage);
			    printf("\n");
			    exit(1);
                        }
                        break;
                    case 'r':
                        n_r = (int) strtol( argv[++i], &endptr, 10);
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
			    printf(error_message);
			    printf(usage);
			    printf("\n");
			    exit(1);
                        }
                        break;
                    case 'z':
                        n_z = (int) strtol( argv[++i], &endptr, 10);
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
			    printf(error_message);
			    printf(usage);
			    printf("\n");
			    exit(1);
                        }
                        break;
                    case 'o':
                        o = (int) strtol( argv[++i], &endptr, 10);
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
			    printf(error_message);
			    printf(usage);
			    printf("\n");
			    exit(1);
                        }
                        break;
                    default:
			printf(error_message);
			printf(usage);
			printf("\n");
			exit(1);
                }
            } else {
		printf(error_message);
		printf(usage);
		printf("\n");
		exit(1);
            } /* end of: if-then-else */
        } /* end-of: for */ 
    } /* end-of: if */

}
/*
Kompilieren:
gcc capacitor.c -l m -o capacitor
Ausfuehren in Verbindung mit gnuplot z.B.:
./capacitor -F 6 -R 6 -r 30 -z 60 -V 500 -o 2 -e 0.5e-2 -U 0 -w 1.7| gnuplot
oder als Bildschirmausgabe:
./capacitor
# Kondensator:
# =========== 
# Gitter: 25 x 25 Punkte
# Kondensatorposition: z=5; Radius: r=5
# Kondensatorspannung: V=1000.000000
# Anfangspotential: Phi0=0.000000
# Relaxationsparameter: Omega=1
# Maximale Abweichung: eps=10

#Iteration      maximale Abweichung                  Energie
#_________ ________________________ ________________________
         1    3.216388702392579e+02    1.057529752577351e+07
         2    1.851851851851852e+02    8.563781479521886e+06
         3    1.156724215534979e+02    7.606258041563846e+06
         4    8.094986443329907e+01    7.027946880043454e+06
         5    6.926442166505618e+01    6.636039406175652e+06
         6    5.711388006216117e+01    6.355297511441418e+06
         7    4.703664059749542e+01    6.148120609742193e+06
         8    3.824758957342215e+01    5.991826313770101e+06
         9    3.180690612096416e+01    5.871424161616748e+06
        10    2.664496281509363e+01    5.776687317526440e+06
        11    2.196427346622710e+01    5.700568349378514e+06
        12    1.817195055240620e+01    5.638180863156057e+06
        13    1.612524443203483e+01    5.586106697652645e+06
        14    1.438133266063471e+01    5.541926048868635e+06
        15    1.289062975161625e+01    5.503901260247231e+06
        16    1.161085742222792e+01    5.470764679563234e+06
        17    1.063083031516152e+01    5.441576057238697e+06
        18    9.797194535102562e+00    5.415626137300598e+06
*/
