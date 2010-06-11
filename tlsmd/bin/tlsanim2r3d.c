/*
 * tlsanim2r3d.c
 *
 * Create a modular Raster3D input file from a list of CA coordinates
 * with associated chain, segment, and model information.
 * We assume that a suitable header file is being created elsewhere.
 *
 * LAST UPDATE: 2010-06-10
 *
 * TODO:
 *   - command-line options for trace radius, etc
 *   - color half-bonds at segment break
 */

#include <stdio.h>
#include <math.h>

typedef struct CA {
    int model;		/* Which chain-trace within the animation */
    char chain;		/* 1 char CHAINID */
    int segment;	/* TLS segment index */
    int libration;	/* 1/2/3  which libration component */
    double x,y,z;	/* The coordinates of this CA atom */
} CA;

typedef struct color {
    double R,G,B;
} color;

/* This color table should match the one in Colors.py used by TLSMD */
static color segcolor[] = {
   {0.750, 0.750, 0.750},	/* Grey */
   {0.000, 0.000, 1.000},	/* Blue */
   {0.000, 1.000, 0.000},	/* Green */
   {1.000, 0.000, 1.000}, 	/* Magenta */
   {1.000, 0.000, 0.000},	/* Red */
   {0.000, 1.000, 1.000},	/* Cyan */
   {1.000, 1.000, 0.000},	/* Yellow */
   {1.000, 0.500, 1.000},	/* Violet */
   {0.500, 0.000, 1.000},	/* PurpleBlue */
   {1.000, 0.600, 0.500},	/* Salmon */
   {0.500, 1.000, 0.500},	/* Lime */
   {0.500, 0.500, 1.000},	/* Slate */
   {0.000, 1.000, 0.500},	/* BlueGreen */
   {1.000, 0.000, 0.500},	/* HotPink */
   {1.000, 0.500, 0.000},	/* Orange */
   {0.500, 1.000, 0.000},	/* YellowGreen */
   {0.500, 0.000, 1.000}, 	/* BlueViolet */
   {0.000, 0.500, 1.000},	/* Marine */
   {0.750, 0.750, 0.000}, 	/* Olive */
   {0.750, 0.000, 0.750},	/* Purple */
   {0.000, 0.750, 0.750}, 	/* Teal */
   {0.500, 0.100, 0.100},	/* Ruby */
   {0.100, 0.500, 0.100},	/* Forest */
   {0.100, 0.100, 0.500},	/* Deep */
   {0,0,0}
};

static double dist(CA *a, CA *b) {
    return (a->x - b->x)*(a->x - b->x)
	 + (a->y - b->y)*(a->y - b->y)
	 + (a->z - b->z)*(a->z - b->z);
}

int main(int argc, char *argv[]) {
    /* Internal constants - may add command line options later */
    double radius = 0.2;	/* Radius of trace in Angstroms */
    double gap = 5.0*5.0;	/* Square of longest allowable CA-CA distance */
    int nca = 0;		/* Number of CAs read so far for this segment */
    int maxcolor = sizeof(segcolor) / sizeof(color);
    int color = 0;
    int desired_libration_group = 0;

    /* Bookkeeping */
    CA previous, current;
    double thickness = radius;
    int ierr;

    /* Informational only */
    fprintf(stderr,"tlsanim2r3d version 0.2\n");
    fprintf(stderr,"\tmaxcolors = %d\n",maxcolor);

    while (8 == (ierr = scanf("%d %d %c %d %d %lf %lf %lf",
			&desired_libration_group,
			&current.model, &current.chain, &current.segment,
			&current.libration, &current.x, &current.y, &current.z))) {
	if (0) printf("%2d %2d %1c %d %d %g %g %g\n",
			desired_libration_group,
			current.model, current.chain, current.segment,
			current.libration, current.x, current.y, current.z);

	/* Ignore unwanted libration groups.
	 * So far there is no way to specify which one we want,
	 * so arbitrarily take the first
	 */
	if (!desired_libration_group) {
	    continue;
	}

	/* Start new segment */
	if (nca == 0) {
	    previous = current;
	    nca = 1;
	    continue;
	}

	/* Check for end of model or chain */
	if ((current.model != previous.model)
	||  (current.chain != previous.chain)) {
	    previous = current;
	    nca = 1;
	    continue;
	}

	/* Check for gap in chain */
	if (dist(&current,&previous) > gap) {
	    previous = current;
	    nca = 1;
	    continue;
	}

	/* OK, we have both ends of a chain-trace segment */
	color = current.segment % maxcolor;
	thickness = (current.model == 0) ? radius : radius * 1.1;
	printf("3\n %.3f %.3f %.3f  %4.2f  %.3f %.3f %.3f  %4.2f",
		previous.x, previous.y, previous.z, thickness,
		current.x, current.y, current.z, thickness);
	printf("   %5.3f %5.3f %5.3f\n",
		segcolor[color].R, segcolor[color].G, segcolor[color].B); 
	previous = current;
	nca++;
    }

    return 0;
}
