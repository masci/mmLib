
/* incomplete gamma functions */

/* Returns the value ln[gamma(XX)] for XX>0.  Full accuracy is obtained for XX>1.
 * For 0<XX<1, the reflection formula can be used first.
 */
static void
gammaln(double xx, double *retval) {
	int j;
	double x, ser, tmp;
	
	double cof[6] = {76.18009173, -86.50532033, 24.01409822, 
			-1.231739516, 0.120858003E-2, -0.536382E-5};
	double stp  = 2.50662827465;
	double half = 0.5;
	double one  = 1.0;
	double fpf  = 5.0;

	x = xx - one;
	tmp = x + fpf;
	tmp = (x + half) * log(tmp) - tmp;
	ser = one;

	for (j = 0; j < 6; j++) {
		x += one;
		ser += cof[j] / x;
	}

	*retval = tmp + log(stp*ser);
}

static void
gser(double a, double x, double *gamser, double *gln) {
	int i;
	double ap, sum, del;

	int itmax = 100;
	double eps = 3.0E-7;
 
	if (x <= 0.0) {
		if (x < 0.0) return;
		*gamser = 0.0;
		return;
	}

	ap = a;
	sum = 1.0 / a;
	del = sum;

	for (i = 0; i < itmax; i++) {
		ap += 1.0;
		del = del * x/ap;
		sum += del;

		if (fabs(del) < fabs(sum)*eps) break;
	}

	gammaln(a, gln);
	*gamser = sum * exp(-x + a*log(x) - *gln);
}

static void
gcf(double a, double x, double *gamser, double *gln) {
	int n;
	double gold, a0, a1, b0, b1, fac, an, ana, g, anf;
	
	int itmax = 100;
	double eps = 3.0E-7;

	g    = 0.0;
	gold = 0.0;
	a0   = 1.0;
	a1   = x;
	b0   = 0.0;
	b1   = 1.0;
	fac  = 1.0;

	for (n = 1; n <= itmax; n++) {
		an = n;
		ana = an - a;
		a0 = (a1 + a0*ana) * fac;
		b0 = (b1 + b0*ana) * fac;
		anf = an * fac;
		a1 = x*a0 + anf*a1;
		b1 = x*b0 + anf*b1;

		if (a1 != 0.0) {
			fac = 1.0 / a1;
			g = b1*fac;
			if (fabs((g-gold)/g) < eps) break;
			gold = g;
		}
	}

	gammaln(a, gln);
	*gamser = exp(-x + a*log(x) - *gln) * g;
}

static void 
gammaq(double a, double x, double *retval) {
	double gamser, gammcf, gln;
	
	if (x < 0.0 || a < 0.0) return;

	if (x < a+1.0) {
		gser(a, x, &gamser, &gln);
		*retval = 1.0 - gamser;
	} else {
		gcf(a, x, &gammcf, &gln);
		*retval = gammcf;
	}
}
