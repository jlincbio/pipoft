use 5.016;
use strict;
use warnings;

use constant NA         => 'NA';
use constant EPS        => 2.220446e-16;
use constant EPS_LOG    => -36.04365;
use constant EPS_LOGMIN => -999;
use constant EPS_LOGMAX => 999;
use constant EPS_ZERO   => 0.00001;
use constant EPS_ITER   => 0.00000001;
use constant INFINITE   => 'Inf';
use constant FMAX_SAMP  => 100000;
use constant FMIN       => 1/65535;
use constant FMAX       => 65535;
use constant REF_KEGGDB => 'KEGG_DATASET_PREFIX';

sub check {
  my $m = shift;
  eval "use $m";
  if (($@) && ($@ =~ m/^Can\'t\ locate/)) {
	  print $@;
	  return 1;
  }
  return 0;
}

my @reqs = (
	'Cwd',
	'File::Fetch',
	'File::Temp',
	'Storable',
	'Getopt::Long',
	'Parallel::ForkManager',
	'Scalar::Util',
	'List::Util',
	'Statistics::TTest',
	'File::Type',
	'PerlIO::gzip',
#	'blahblahblahyadayadayada',
	'XML::Simple'
);

my @missing; my $err = 0;
foreach my $i (@reqs) {
	my $s = check $i;
	if ($s) {
		push @missing, $i;
		$err++;
	}
}

print "\n";

if (scalar @missing) {
	print "----------------------------------------\n";
	print "      Missing Perl Packages Found:      \n\n";
	print "      $_\n" for @missing;
	print "\n----------------------------------------\n\n";
}

exit $err;