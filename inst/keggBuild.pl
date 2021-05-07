use 5.016;
use strict;
use warnings;
use Parallel::ForkManager;
use File::Fetch;
use File::Temp 'tempdir';
use Storable;
use Getopt::Long;

# downloads KEGG pathways for a specified organism (human) to build
# a database of gene lists in plain text format and in binary

$File::Fetch::WARN = 0;

sub getCores {
	my $cpuCount = shift;
	my $smp = 1; my $fixed = 2;
	if (defined $cpuCount) {
		$cpuCount = int $cpuCount;
		$fixed = 1;
		$smp = 0 if ($cpuCount < 2);
	} else {
		if (-e '/proc/cpuinfo') {
			open my $CPUIN, '<', '/proc/cpuinfo';
			$cpuCount = 0;
			while (<$CPUIN>) {
				$cpuCount++ if ($_ =~ m/^processor/);
			}
			close $CPUIN;
		} elsif ($^O eq 'darwin') {
			$cpuCount = `sysctl -n hw.ncpu`;
		} elsif ($^O eq 'MSWin32') {
			$cpuCount = $ENV{NUMBER_OF_PROCESSORS} // 1;
		} else {
			$cpuCount = 0;
		}
		chomp $cpuCount;
		if ($cpuCount < 2) {
			$cpuCount = 0; $smp = 0;
		};
	}	
	return sprintf("%u", 0.6 * $cpuCount);
}
sub quickname {
	my ($prefix, $suffix, $fine) = @_;
	my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
	$year = ($year % 100);
	$mon = sprintf "%02u", $mon + 1;
	my $day = sprintf "%02u", $mday;
	my $strtime = join '', $mon, $day, $year;
	if (defined $fine) {
		# include hours and minutes
		$strtime = $strtime.$hour.$min;
	}
	my $name = $prefix.'_'.$strtime.'.'.$suffix;
	return $name;
}
my ($opt_prefix, $opt_org, $opt_core, $tmp_path) = @ARGV;

# probably need to get a list of organisms here for later
$opt_org = 'hsa' if (lc $opt_org ne 'hsa'); 
print "KEGG organism specified: $opt_org\n";

my $opCores = getCores($opt_core);
print "Number of threads to be used for data retrieval and processing: $opCores\n";
print "Retrieving a list of KEGG pathways...\n";
my $org_fetch = File::Fetch->new(uri => "http://rest.kegg.jp/list/pathway/$opt_org");
my $org_path = $org_fetch->fetch(to => $tmp_path) or die $org_fetch->error;
my %pathways; my %genes; my $i = 0;
open my $LIST_PATHWAYS, '<', $org_path or die "Can't read pathway information!\n";
while (<$LIST_PATHWAYS>) {
	chomp;
	my @s = split "\t", $_;
	$s[0] =~ s/^path\://;
	$pathways{$s[0]} = $s[1];
}
close $LIST_PATHWAYS;

my $pm = new Parallel::ForkManager($opCores); $i = 0;
$pm-> run_on_finish(sub{
	my @stats = @_;
	my @kegg = @{$stats[5]};
	if (@{$kegg[1]}) {
		# genes present
		$genes{$kegg[0]} = $kegg[1];
	} else {
		delete $pathways{$kegg[0]};
	}
});
print "Retrieving gene lists for KEGG pathways...\n";
foreach my $p (keys %pathways) {
	$pm->start and next;
	sleep 1;
	my @symbols; my $err_fetch = 0; my $err_path = 0;
	my $list_fetch = File::Fetch->new(uri => "http://rest.kegg.jp/get/$p");
	my $list_path = $list_fetch->fetch(to => $tmp_path) or do {
		printf "Error: %s\n", $list_fetch->error;
		$err_fetch = 1;
	};
	
	if ($err_fetch == 0) {
		# continue only if file is fetched
		open my $PATHWAY, '<', $list_path or do {
			print "Error: cannot read from retrieved file for $p!\n";
			$err_path++;
		};
		if ($err_path == 0) {
			my $l = 1; my $start;
			while (<$PATHWAY>) {
				chomp;
				$start = $l if (m/^GENE/);
				if (defined $start) {
					my @g = split ' ';
					if ($l == $start) {
						shift @g;
					}
					$l++;
	
					if ($g[0] =~ m/^\d/) {
						$g[1] =~ s/;//;
						push @symbols, $g[1];
					}
				}
			}
			close $PATHWAY;
		}
	}
	my @passback = ($p, \@symbols);
	my $err_total = $err_fetch + $err_path;
	print "Error encountered during retrieval of $p data!\n" if ($err_total > 0);
	$pm->finish($err_total, \@passback);
}
$pm->wait_all_children;

my $output_genes = quickname($opt_prefix, 'geneKegg');
my $output_paths = quickname($opt_prefix, 'pathKegg');
my $output_debug = quickname($opt_prefix, 'txt');
if (-e $output_genes) {
	# if file already exists, change to fine mode
	$output_genes = quickname($opt_prefix, 'geneKegg', 1);
	($output_paths = $output_genes) =~ s/geneKegg$/pathKegg/;
	($output_debug = $output_genes) =~ s/geneKegg$/txt/;
}

open my $OUTPUT_TEXT, '>', $output_debug or die "Can't write output to disk!\n";
foreach my $k (keys %genes) {
	print $OUTPUT_TEXT "\*\*$k\n";
	print $OUTPUT_TEXT "\#\#$pathways{$k}\n";
	print $OUTPUT_TEXT "$_\n" for @{$genes{$k}};
	print $OUTPUT_TEXT "\/\/\n";
}
close $OUTPUT_TEXT;

store \%genes, $output_genes;
store \%pathways, $output_paths;
sleep 1;

print "Output KEGG pathway data: $output_paths\n";
print "Output KEGG gene lists  : $output_genes\n";
print "Output KEGG summary file: $output_debug\n";
exit 0;