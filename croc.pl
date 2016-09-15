#!/usr/bin/perl -w

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

# User modules:
use lib qw(./lib);
use multtest qw/bonferroni benjamini_hochberg benjamini_liu holms/;
## hypergeom is loaded dynamically

our $VERSION = 0.01;

###########################################
## Options
my %chr_ref;
my $regulated_file;
my $window=30000;
my $gwindow;
my $offset=10000;
my $pval_lim = 0.05;
my $min_genes_in_cluster = 3;
my $help = 0;
my $version = 0;
my $ref_file = "";
my $fdr_method = "BH";
my $allPerl = 0;
my $genomeRef = 0;
my $options = GetOptions ("r|reg=s"         => \$regulated_file,
						  "f|ref:s"         => \$ref_file,
						  "w|window:i"      => \$window,
						  "g|gwindow:i"     => \$gwindow,
						  "o|offset:i"      => \$offset,
						  "p|pval:f"        => \$pval_lim,
						  "m|min_genes:i"   => \$min_genes_in_cluster,
						  "t|fdr:s"         => \$fdr_method,
						  "a|genomeref"     => \$genomeRef,
						  "v|version"       => \$version,
						  "h|help"          => \$help,
						  "allPerl"         => \$allPerl);

if ($allPerl){
	require hypergeom_allPerl;
	import hypergeom_allPerl;
} else {
	require hypergeom;
	import hypergeom;
}

my %opts = (
	    'window' => $window,
	    'offset' => $offset,
	    'gwindow' => $gwindow
	   );


if ( (!defined $regulated_file) || ($help) ){
	phelp();
	exit(0) if ($help);
	exit(1);
}

if ($version){
	pversion();
	exit (0);
}

if ($fdr_method !~ /^B[OHL]$/ and $fdr_method !~ /^[HN]O$/){
	print "$0 : fdr method $fdr_method is not recognised. Possibilities are: BO, HO, NO, BH and BL\n";
	phelp();
	exit (1);
}

# Read the reference
open my $p, "<", $ref_file or die $!;
my $chr_lengths = get_chrms_lenghts($p);


if (!defined $gwindow){
	print STDERR "### Options:\n# Reg_file:$regulated_file\nRef_file:$ref_file\n# window_size=$window\n# offset=$offset\n# pval=$pval_lim\n# min_genes=$min_genes_in_cluster\n\n";
} else {
	print STDERR "# Reg_file:$regulated_file\nRef_file:$ref_file\n# window_size=$gwindow genes\n# pval=$pval_lim\n# min_genes=$min_genes_in_cluster\n\n";
}

my $arrayRegGenes_ref = extract_genes ($regulated_file);
my ($chrs,$n_all_genes,$n_reg_genes, $unknown_genes) = get_all_genes_by_chrs2 ($p,$arrayRegGenes_ref);
my @clusters=();
for my $chr (sort keys %{$chrs}){
#  pr_report ($chr,$n_all_genes,$n_reg_genes);
  my $cluster = cluster_finder (${$chrs}{$chr},${$chr_lengths}{$chr},${$n_all_genes}{$chr},${$n_all_genes}{'all'},${$n_reg_genes}{$chr},${$n_reg_genes}{'all'},\%opts,$chr,$fdr_method) if ((${$n_reg_genes}{$chr}) && (${$n_reg_genes}{$chr} >= $min_genes_in_cluster));
  push @clusters,$cluster if (defined $cluster);
}
print STDERR "\n";
print STDERR "Unknown genes: ";
print STDERR join (" ", keys %$unknown_genes);
print STDERR "\n\n";

#Dump \@clusters;

for my $cl (@clusters){
  print_clusters ($cl);
}

########## SUBROUTINES ##########

sub pr_report
{
  my ($chr,$n_all_genes,$n_reg_genes) = @_;
  my $porc_ch = sprintf "%.1f", 100 * (${$n_all_genes}{$chr} || 0) / ${$n_all_genes}{'all'};
  print "$chr:" . (${$n_all_genes}{$chr} || 0) .":$n_all_genes->{all}:$porc_ch:";
  if (defined ${$n_reg_genes}{$chr}){
    my $porc_reg_ch = sprintf "%.1f", 100 * (${$n_reg_genes}{$chr} / ${$n_reg_genes}{'all'});
    print STDERR "${$n_reg_genes}{$chr}:$n_reg_genes->{all}:$porc_reg_ch:\n";
  } else {
    print STDERR "0\n";
  }
}

sub get_all_genes_by_chrs2
  {
    my ($f_all,$arr_reg) = @_;
    my %seen_refseqs = ();
    my %reg_genes = map {$_ => 1} @{$arr_reg}; ## Pueden ser names or refseqs
    my %chrs;
    my %unknown_genes = %reg_genes;
    my %last_gene = ();
    my %n_reg_genes;
    my %n_all_genes;
    while (my $refgene_line = <$f_all> ){
      next if ($refgene_line =~ /^#/);
      chomp $refgene_line;
      my %gene;
      my ($gene, $chrm, $begin, $end, $name) = split /\t/,$refgene_line;
      if (! $seen_refseqs{$gene}){
	$seen_refseqs{$gene} = 1;
	@gene{qw/gene chrm begin end name/} = ($gene, $chrm, $begin, $end, $name);
	$gene{is_reg} = $reg_genes{$gene} || $reg_genes{$name} || 0;
	if (! %last_gene){ %last_gene = %gene; delete $unknown_genes{$gene}; next; }
	delete $unknown_genes{$gene} if (defined $unknown_genes{$gene});
	delete $unknown_genes{$name} if (defined $unknown_genes{$name});
	if (defined $last_gene{name} and $gene{name} eq $last_gene{name}){
	  %last_gene = join_genes(\%last_gene, \%gene);
	} else {
	  push @{$chrs{$last_gene{chrm}}},{%last_gene};
	  $n_all_genes{$chrm}++;
	  $n_all_genes{'all'}++;
	  if ($last_gene{is_reg}){
	    $n_reg_genes{$last_gene{chrm}}++;
	    $n_reg_genes{'all'}++;
	  }
	  %last_gene = %gene;
	}
      }
    }
    push @{$chrs{$last_gene{chrm}}},{%last_gene};
    $n_all_genes{$last_gene{chrm}}++;
    $n_all_genes{'all'}++;
    $n_reg_genes{$last_gene{chrm}}++ if ($last_gene{'is_reg'});
    $n_reg_genes{'all'}++ if ($last_gene{'is_reg'});
    close ($f_all);
    return (\%chrs,\%n_all_genes,\%n_reg_genes, \%unknown_genes);
  }

sub join_genes
  {
    my ($ref1,$ref2) = @_;
    my %gene;
    @gene{qw/name strand chrm cdsStart cdsEnd exonStarts exonEnds/} = @{$ref2}{qw/name strand chrm cdsStart cdsEnd exonStarts exonEnds/};
    $gene{begin} = $ref1->{begin} < $ref2->{begin} ? $ref1->{begin} : $ref2->{begin};
    $gene{end}   = $ref1->{end}   > $ref2->{end}   ? $ref1->{end}   : $ref2->{end};
    $gene{gene}  = $ref1->{gene} .";".$ref2->{gene};
    $gene{is_reg} = $ref1->{is_reg} || $ref2->{is_reg} || 0;
    return %gene;
  }

sub extract_genes
{
    my $file = shift @_;
    open (F, $file) or die "Cannot find file $file\n\n";

    my %gene_refseqs = ();
    while (my $refseqs = <F>){
      next if $refseqs=~/^$/;
      chomp $refseqs;
#      $loadedGenes++;
      $gene_refseqs{$refseqs} = 1;
    }
    close (F);
    return [keys %gene_refseqs];
}

sub cluster_finder
  {
    my ($ref_ch, $length, $N, $Nall, $k, $kall, $opts_ref, $chr_name, $fdr_meth) = @_;
    do {$N=$Nall; $k=$kall} if ($genomeRef);
    my $n_clusters = 0;
    my @clusters;
    my @pvals;

    my $next_genes;
    if (defined ${$opts_ref}{'gwindow'}){
      $next_genes = get_next_genes_gwindow ($ref_ch,${$opts_ref}{'gwindow'});
    } else {
      $next_genes = get_next_genes ($ref_ch, ${$opts_ref}{'window'}, ${$opts_ref}{'offset'});
    }

    while (my ($M_genes, $x_genes) = $next_genes->()){
      last unless defined $M_genes;
      my $M = scalar (@{$M_genes});
      my $x = scalar (@{$x_genes});
## Perform hypergeometric distribution test if there are enough genes
      my $hpobj = $allPerl ? hypergeom_allPerl -> new($N,$M,$k,$x) : hypergeom -> new ($N, $M, $k, $x);
      my $pvals = $hpobj -> over_represented_pval ();
      push @pvals, $pvals;

      if ( $pvals < $pval_lim and $x>=$min_genes_in_cluster){
		$n_clusters++;
		my %cluster;
		$cluster{"chr"}    = $chr_name;
		$cluster{"begin"}  = ${$x_genes}[0]{'begin'};
		$cluster{"end"}    = ${$x_genes}[-1]{'end'};
		$cluster{"genes"}  = [map {$_->{name}} @{$x_genes}];
		$cluster{"pvalue"} = $pvals;
		push @clusters, {%cluster};
      }
    }
my %multtests = (
	'BO' => \&bonferroni,
	'BH' => \&benjamini_hochberg,
	'BL' => \&benjamini_liu,
	'HO' => \&holms,
	'NO' => sub {$pval_lim}
	);
## Performs multtest correction        
	my $qlim = $multtests{$fdr_meth}->($pval_lim,\@pvals);
	print STDERR "$chr_name : QLIM => $qlim\n";
	my @goodClusters = grep {$_->{"pvalue"} <= $qlim} @clusters;
    return \@goodClusters;
  }

sub get_next_genes_gwindow
  {
    my ($arr_ref, $gwindow) = @_;
    my ($from,$to) = (0,$gwindow-1);

    my $getnext = sub
      {
	return (undef,undef) if ($from+$to >= scalar @{$arr_ref});
	my @genes = @{$arr_ref}[$from..$from+$to];
	my @genes_reg = grep {${$_}{'is_reg'} == 1} @genes;
	$from++;
	return ([@genes],[@genes_reg]);
      };
    return $getnext;

  }

sub get_next_genes
  {
    my ($arr_ref, $window, $offset) = @_;
    my $pos2beg = 0;
    my $p = 0;

    my $getnext = sub
      {
	return (undef,undef) if ($pos2beg==scalar @{$arr_ref});
	$p+=$offset;
	my @genes = ();
	my @genes_reg = ();
	for (my $i = $pos2beg; $i < scalar @{$arr_ref}; $i++){
	  $pos2beg++ if (${$arr_ref}[$i]{'begin'} < ($p - $window + $offset));
	  return (\@genes, \@genes_reg) if (${$arr_ref}[$i]{'begin'} > ($p));
	  push @genes, ${$arr_ref}[$i];
	  push @genes_reg, ${$arr_ref}[$i] if (${$arr_ref}[$i]{'is_reg'} == 1);
	}
	return ([@genes], [@genes_reg])
      };

    return $getnext;
  }

#	for (my $nround = 0; $nround < scalar (keys %{$hash_ref})-$gwindow+1; $nround++){
#	    for (my $ng = 0; $ng < $gwindow; $ng++){

sub print_clusters
  {
    my ($cluster_aref) = @_;
    filter_overlaping_clusters ($cluster_aref);
    my $count = 1;
    foreach my $cluster_ref (@{$cluster_aref}){
    	if (defined $cluster_ref){
    		my $chr = ${$cluster_ref}{"chr"};
			my $begin = ${$cluster_ref}{"begin"};
			my $end = ${$cluster_ref}{"end"};
			my $score = ${$cluster_ref}{"pvalue"};
			$score = $score if ($score =~ /\d+/);
			my @genes = @{${$cluster_ref}{"genes"}};
			if ($score < $pval_lim){
				printf "%s\t%s\t%s\t%d\t%d\t%g\t.\t.\t",
				$chr,
				 "CROC",
				  "cluster",
				  	$begin,
				  	 $end,
				  	  $score;
	  			for (my $i=0; $i<scalar @genes; $i++){
	    			print $genes[$i]." ";
	  			}
	  			print "\n";
	  			$count++;
			}
      	}
    }
  }

sub filter_overlaping_clusters  ### -f --filter_overlaping 
{
    my $good_clusters_aref = shift @_;
    my %non_overlaping;

    for (my $n=0; $n<(scalar (@{$good_clusters_aref})-1); $n++){
	if (defined (${$good_clusters_aref}[$n])){
	    if ( (${${$good_clusters_aref}[$n]}{"chr"}) eq ${${$good_clusters_aref}[($n+1)]}{"chr"}  && 
		 (${${$good_clusters_aref}[$n]}{"end"} > ${${$good_clusters_aref}[($n+1)]}{"begin"}) ){
		${${$good_clusters_aref}[($n+1)]}{"begin"} = ${${$good_clusters_aref}[$n]}{"begin"};
		${${$good_clusters_aref}[($n+1)]}{"pvalue"} = 
		    ${${$good_clusters_aref}[($n+1)]}{"pvalue"}>${${$good_clusters_aref}[($n)]}{"pvalue"}? ${${$good_clusters_aref}[$n]}{"pvalue"}:${${$good_clusters_aref}[($n+1)]}{"pvalue"};

		my %ref_2_genes_in_cluster = map { $_->{"name"} => $_ } @{${${$good_clusters_aref}[$n]}{"ref2genes"}},@{${${$good_clusters_aref}[($n+1)]}{"ref2genes"}};
		@{${${$good_clusters_aref}[($n+1)]}{"ref2genes"}} = sort {$a->{"begin"} <=> $b->{"begin"}} values %ref_2_genes_in_cluster;
		my %genes_in_cluster = map {$_ => 1} @{${${$good_clusters_aref}[$n]}{"genes"}},@{${${$good_clusters_aref}[($n+1)]}{"genes"}};
		@{${${$good_clusters_aref}[($n+1)]}{"genes"}} = keys %genes_in_cluster;

#		foreach my $name (@{${${$good_clusters_aref}[$n]}{"genes"}}){
#		    unshift @{${${$good_clusters_aref}[($n+1)]}{"genes"}}, $name if !&isin_array (${${$good_clusters_aref}[($n+1)]}{"genes"}, $name);
#		}
		undef (${$good_clusters_aref}[$n]);
	    }
	}
    }

#     for (my $n=0; $n<scalar @{$all_clusters_aref}; $n++){
# 	if (defined (${$all_clusters_aref}[$n])){
# 	    if (${${$all_clusters_aref}[$n]}{"end"} < ${${$all_clusters_aref}[$n+1]}{"begin"}){
# 		${${$all_clusters_aref}[$n]}{"end"} = ${${$all_clusters_aref}[$n+1]}{"end"};
# 		${${$all_clusters_aref}[$n]}{"pvalue"} = (${${$all_clusters_aref}[$n]}{"pvalue"} + ${${$all_clusters_aref}[$n+1]}{"pvalue"}) / 2;
# 		foreach my $name (@{${${$all_clusters_aref}[$n+1]}{"genes"}}){
# 		    push @{${${$all_clusters_aref}[$n]}{"genes"}}, $name if !&isin_array (${${$all_clusters_aref}[$n]}{"genes"}, $name);
# 		}
# 		undef (${$all_clusters_aref}[$n+1]);
# 	    }
# 	}
#     }

}

sub get_chrms_lenghts
  {
    my ($fh) = @_;

    local $/ = "\n\n";
    my $chunk = <$fh>;

    my $maxLength = 0;
    my @lines = split /\n/,$chunk;
    shift @lines if ($lines[0] =~ /^#/);
    my %chrs_length = map {my ($chr,$length) = split /\s+/; $maxLength=$length if ($length>$maxLength);$chr => $length} @lines;
    $chrs_length{max} = $maxLength;
    return \%chrs_length;
  }
  
sub phelp
{
print <<EOF;
croc.pl [options] --reg <file> --ref <file> [-options]
Options:
    -w | --window [integer]    : Width of the length to calculate the clusters (Default: 30000)
    -o | --offset [integer]    : Offset between windows (Default: 10000). Only has effect with -w
    -g | --gwindow [integer]   : Defines the window by number of genes instead of by number of nucleotides (and sets offset to 1 gene)
    -p | --pval [float]        : Max allowed pvalue per cluster before multiple testing correction (Default: 0.05)
    -m | --min_genes [integer] : Min number of genes that can for a cluster
    -a | --genomeref           : Tells the program to use the entire genome as reference in the statistical calculations (by default, the references are each chromosome)
    -t | --fdr [string]        : Method to apply to correct multiple testing. Possibilities are:
                                  BO => Bonferroni
                                  HO => Holms
                                  BH => Benjamini & Hochberg (default)
                                  BL => Benjamini & Liu
                                  NO => No correction
    --allPerl                  : Use pure Perl code to perform the hypergeometric distribution tests
                                 Using this option, the script will run slower, but will run in systems where you can't use Inline::C
    -h | --help                : Displays this document
EOF
}

sub pversion
{
	print "croc.pl v$VERSION\n";
}