package multtest;

use warnings;
use strict;
use Data::Dumper;

use Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw/&bonferroni &benjamini_hochberg &benjamini_liu &holms/;

###################################################
##                                               ##
## For the thory behind these procedures see:    ##
##                                               ##
## Benjamini et al. Behavioural Brain Research   ##
## (2001) 125:279-284                            ##
##                                               ##
###################################################

sub bonferroni
{
	my ($plim,$prefs) = @_;
	if (@$prefs) {
		return $plim / @$prefs;
	} else { return 0 }
}

sub rholms
{
	my ($plim,$prefs) = @_;
	my $b = bonferroni ($plim,$prefs);
	if ( $b < $prefs->[0]){
		return $b;
	} else {
		shift @$prefs;
		rholms ($plim,$prefs)
	}
}

sub holms
{
	my ($plim,$prefs) = @_;
	my @sprefs = sort {$a <=> $b} @$prefs;
	return rholms ($plim,\@sprefs);
}

sub benjamini_hochberg
{
	my ($plim,$prefs) = @_;
	my @sprefs = sort {$b <=> $a} @$prefs;
	my $l = @sprefs;
	for my $pos (0..$#sprefs){
		my $c = $plim*($pos+1)/$l;
		return $c if ($sprefs[$pos] < $c);
	}
	return 0;
}

sub benjamini_liu
{
	my ($plim,$prefs) = @_;
	my @sprefs = sort {$a <=> $b} @$prefs;
	my $l = @sprefs;
	for my $pos (0..$#sprefs){
		my $c = $l+1-($pos+1);
		my $c1 = $plim * $l /($c*$c);
		my $co = $c1 < 0.05 ? $c1 : 0.05;
		return $co if ($sprefs[$pos] > $co);
	}
	return 0;
}

1;
