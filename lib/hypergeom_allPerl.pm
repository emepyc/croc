package hypergeom_allPerl;

use warnings;
use strict;
use List::Util qw/min/;

sub gammaln
{
    my $x = shift @_;
    my ($y, $tmp, $ser);

    my @cof = qw(76.18009173 -86.50532033 24.01409822 -1.231739516 0.120858003e-2 -0.536382e-5);
    my $i;

    return 0 if ($x == 0 || $x == 1);

    $y = $x-1;
    $tmp = $y + 5.5;
    $tmp -= ($y+0.5)*log ($tmp);
    $ser = 1;
    for ($i=0; $i<5; $i++){
	$y+=1;
	$ser+=$cof[$i]/$y;
    }

    return log($x) + (log(2.50662827465*$ser)-$tmp);
}

sub choose
{
    my ($A, $B) = @_;

    my $choose_res = &gammaln($A) - (&gammaln($B)+&gammaln($A-$B));
    return $choose_res;
}


sub hypergeom_c
{
    my ($N, $M, $k, $x) = @_;

    my $prob = 0;
    my $lprob;
    my $pvalue = 0;

    $lprob = &choose($M,$x) + &choose(($N-$M),($k-$x)) - &choose($N,$k);

    for (my $i=$x; $i<=min($k,$M); $i++){
	$lprob = &choose($M,$i) + &choose(($N-$M),($k-$i)) - &choose($N,$k);
	$pvalue += exp($lprob);
    }
    return $pvalue;

}

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
