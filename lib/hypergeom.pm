package hypergeom;

use warnings;
use strict;

use Inline (
	    Config =>
	    DIRECTORY => 'lib/.Inline',
	   );

use Inline Config    => (
			 'FORCE_BUILD'       => 0,
			 'CLEAN_AFTER_BUILD' => 0,
			);

use Inline C => << 'END_C_CODE';
/** Begin of the C Code **/

#define MIN(x,y) ((x<y)?x:y)
#define MAX(x,y) ((x>y)?x:y)

double gammaln(double x) {
    double y,tmp,ser,
	cof[6]={76.18009173,-86.50532033,24.01409822,
		-1.231739516,0.120858003e-2,-0.536382e-5};
    int    i;
    
    if (x==0.0 || x==1.0)
	return 0.0;
    
    y=x-1.0;
    tmp=y+5.5;
    tmp-=(y+0.5)*log(tmp);
    ser=1.0;
    for (i=0;i<=5;i++) {
	y+=1.0;
	ser+=cof[i]/y;
    }
    
    return log(x) +(log(2.50662827465*ser)-tmp);
}

double choose(int A, int B){

    double choose_res = gammaln(A) - (gammaln(B)+gammaln(A-B));

//    printf ("(%d %d) = %.5f\n",A,B,choose_res);
    return choose_res;
}

double hypergeom_c(int N, int M, int k, int x)
{
    double prob = 0;  /* Final probability */
    double lprob; /* ln(prob) */
    double pvalue = 0;
    double i;     /* Counter */


    lprob = choose(M,x)+choose((N-M),(k-x))-choose(N,k);
    prob = exp(lprob);

//    printf ("\nprob = %.10f\n",prob);
    
    for (i=x; i<=MIN(k,M); i++){
	lprob = choose(M,i)+choose((N-M),(k-i))-choose(N,k);
	pvalue += exp(lprob);
    }
//    printf ("pvalue = %.10f\n",pvalue);

    return pvalue;
}

/** End of C Code **/
END_C_CODE

#######################################################
############## Continue Perl Code #####################
#######################################################

sub new
{
    my ($class) = @_;
    bless {
	_N => $_[1],
	_M => $_[2],
	_k => $_[3],
	_x => $_[4],
    }, $class;
}

sub over_represented_pval
{
    my ($self) = @_;
    return 1 if (
		! $self->{_N} ||
		! $self->{_M} ||
		! $self->{_k} ||
		! $self->{_x});
    my $pv = hypergeom_c ($self->{_N},$self->{_M},$self->{_k},$self->{_x});
    return $pv;
}


1
