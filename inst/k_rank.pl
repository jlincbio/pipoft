use 5.016;
use strict;
use warnings;
use Storable;
use Getopt::Long;
use Statistics::TTest;
use File::Type;
use PerlIO::gzip;
use Cwd;
use Scalar::Util 'looks_like_number';

use constant EPS        => 2.220446e-16;
use constant EPS_LOG    => -36.04365;
use constant NA         => 'NA';

sub gzInput {
	my $file = shift;
	my $ftype = File::Type->new->checktype_filename($file);
	my $FHANDLE;
	if ($ftype =~ m/gzip/) {
		open $FHANDLE, '<:gzip', $file or die "Error opening $file\!\n";
	} else {
		open $FHANDLE, '<', $file or die "Error opening $file\!\n";
	}
	return $FHANDLE;
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
sub mean {
	my @array = @_;
	my $count = scalar @array; my $sum = 0;
	if ($count >= 1) {
		foreach my $i (@array) {
			$sum = $sum + $i;
		}
		return ($sum/$count);
	}
	return NA;
}
sub min {
	my @x = @_;
	my $m = $x[0];
	foreach my $i (@x) {
		if ($i < $m) {
			$m = $i;
		}
	}
	return $m;
}
sub log_odds {
	my ($x11, $x12, $x21, $x22) = @_;
	return log((($x11 + 0.01)/($x12 + 0.01))/(($x21 + 0.01)/($x22 + 0.01)));
}
sub log_limit {
	my $x = shift;
	if ($x < EPS) {
		return EPS_LOG;
	}
	return ((log $x)/(log 10));
}
sub significance {
	my ($alpha, @p) = @_;
	my $p = log_limit($alpha);
	my $n = 0;
	foreach my $m (@p) {
		$n++ if ($m < $p);
	}
	if ($n >= (scalar @p) * 0.6) {
		return 'Significant';
	} elsif ($n > 0) {
		return 'Marginal';
	}
	return 'N/S';
}
sub read_list {
	my ($file, $col) = @_; my (%genes, %empty);
	if (defined $file) {
		print STDERR "Reading from $file...\n";
		if (-e $file) {
			my $INPUT_LIST = gzInput($file);
			while (<$INPUT_LIST>) {
				chomp;
				next if ($_ =~ m/^\#/);
				if (defined $col) {
					# bed mode
					my @s = split "\t", $_;
					my $i = $col - 1;
					$i = 3 if ($col > scalar @s);
					next if ($s[$i] eq 'Intergenic');
					$s[$i] =~ s/Promoter\_//;
					$genes{$s[$i]}++;
				} else {
					$genes{$_}++;
				}
			}
			close $INPUT_LIST;
			return %genes;	 
		}
	}
	return %empty;
}
my $usage = join "\n", (
	'k_rank [options] --target <GENE1,GENE2,...> --kegg <DATABASE> --input <ARRAY_DATA>',
	'Parameters:',
	'   --kegg   KEGG gene database file (geneKegg)',
	'   --input  Processed 6-column microarray logFC table file',
	'   --target List of genes assumed to be on-target',
	'   --sites  PI polyamide binding data, Chem-seq or from motifBed (BED)',
	'   --prefix Prefix for the output data',
	'   --alpha  Significance level (default 0.05)'
);

die "Usage: $usage\n" unless (@ARGV);

my $opt_keggfs; # geneKegg file
my $opt_return;
my $opt_foreground;
my $opt_genes;
my $opt_bed;
my $opt_alpha;
my $opt_path;
my $opt_prefix = 'keggRank';

my (@headers, @entries);

GetOptions(
	'kegg=s' => \$opt_keggfs,
	'input=s' => \$opt_foreground,
	'target=s' => \$opt_genes,
	'sites=s' => \$opt_bed,
	'prefix=s' => \$opt_prefix,
	'path=s' => \$opt_path,
	'alpha=s' => \$opt_alpha,
	'rreturn' => \$opt_return);

unless ((defined $opt_genes) && (defined $opt_keggfs) && (defined $opt_foreground)) {
	print STDERR "Required options not specified!\n";
	exit 2;
}

$opt_path = cwd() unless (defined $opt_path);
$opt_alpha = 0.05 unless ((defined $opt_alpha) && (looks_like_number $opt_alpha));
push @headers, sprintf("Significance level set as log10 p < %.4f", log_limit($opt_alpha));

my $ref_genes; my $input_genes = $opt_keggfs;
print STDERR "Reading gene lists...\n";
if (-e $input_genes) {
	$ref_genes = retrieve $input_genes;
} else {
	print STDERR "Error: KEGG library pathway data cannot be found or is incomplete!\n";
	exit 1;
}

my (%key_fc, %key_pathways, %key_genes, %sites);
if (defined $opt_bed) {
	if (-e $opt_bed) {
		my $INPUT_BED = gzInput($opt_bed);
		while (<$INPUT_BED>) {
			chomp;
			my @s = split '\t', $_;
			(my $sym = $s[3]) =~ s/^[^_]*_//;
			next if ($sym eq 'Intergenic');
			$sites{$sym}++;
		}
		close $INPUT_BED;
	} else {
		die "Can't locate $opt_bed\! Check file path!\n";
	}
}

foreach my $gene (keys %sites) {
	# maps sites into %key_genes, but mark with 1 to denote data from bed input
	$key_genes{$gene} = 1;
}

my @primary = split '\,', uc $opt_genes;
print STDERR "Reference genes: ";
foreach my $rgene (@primary) {
	print STDERR "$rgene ";
	$key_genes{$rgene} = 2;
}
print STDERR "\n";

my %pathways; my %kegg = %{$ref_genes};
foreach my $p (keys %kegg) {
	my @genes = @{$kegg{$p}};
	foreach my $g (@genes) {
		# compiling a index of pathways per gene
		push @{$pathways{$g}}, $p;
		if ((exists $key_genes{$g}) && ($key_genes{$g} == 2)) {
			$key_pathways{$p}++;
		}
	}
}

my (%exp, %means); my $n = 0;
my $INPUT_FG = gzInput($opt_foreground);
#Chrom	Start	End	Symbol	RefSeq	logFC
while (<$INPUT_FG>) {
	chomp; next if (($n == 0) && ($_ =~ m/^Chrom\tStart/));
	my @s = split '\t', $_;
	my $sym = $s[3]; my $fc = pop @s;
	next unless (looks_like_number $fc);
	if (exists $key_genes{$sym}) {
		$key_fc{$sym} = $fc;
	}
	if (exists $pathways{$sym}) {
		my @k = @{$pathways{$sym}};
		foreach my $p (@k) {
			push @{$exp{$p}}, $fc;
			
		}
	}
}
close $INPUT_FG;

my (@fgPath, @bgPath); my $nfg = 0; my $nbg = 0;
foreach my $p (keys %exp) {
	$means{$p} = mean(@{$exp{$p}});
	if ((exists $key_pathways{$p}) && ($key_pathways{$p} > 0)) {
		push @fgPath, $means{$p};
		$nfg++;
	} else {
		push @bgPath, $means{$p};
		$nbg++;
	}
}
push @headers, "Enrichment significance: $nfg foreground pathways (background = $nbg)";
my $ttest = new Statistics::TTest;  
$ttest->set_significance(99);
$ttest->load_data(\@fgPath,\@bgPath);  

push @headers, (
	sprintf("Mean expression change (FG) = %.4f", mean(@fgPath)),
	sprintf("H0: mean(FG) = mean(BG) -> t = %.4f", $ttest->t_statistic),
	sprintf("H0: mean(FG) = mean(BG) -> log10 p = %.4f", log_limit($ttest->{t_prob})),
	sprintf("H0: var(FG) = var(BG) -> F = %.4f", $ttest->f_statistic)
);

my %refPvals; my %refMeans; my %refRanks;
foreach my $rgene (keys %key_fc) {
	my %ranks;
	foreach my $p (@{$pathways{$rgene}}) {
		my @fc = sort {$a <=> $b} @{$exp{$p}};
		if (scalar @fc > 2) {
			my %index;
			@index{@fc} = (1..scalar @fc);
			my $i;
			if (exists $key_fc{$rgene}) {
				my $i = ($index{$key_fc{$rgene}} / (scalar @fc));
				$ranks{$p} = $i if (defined $i);
			}
		}
	}
	if (scalar keys %ranks > 2) {
		my $ttest = new Statistics::TTest;  
		$ttest->set_significance(99);
		my @fgRef = values %ranks;
		my @bgRef = (0.5) x (scalar @fgRef); 
		$ttest->load_data(\@fgRef, \@bgRef);
		$refPvals{$rgene} = sprintf("%.4f", log_limit($ttest->{t_prob}));
		$refMeans{$rgene} = sprintf("%.4f", mean(@fgRef));
		$refRanks{$rgene} = \@fgRef;
	} 
}

my %pathOut;
my @pctTarget; my %pctBound; my %smed;
foreach my $r (keys %refPvals) {
	if ($key_genes{$r} == 2) {
		# primary target as specified by input
		my $naRanks = join "\t", ('NA') x (scalar @primary);
		push @entries, "$r\tTarget\t$refMeans{$r}\t$refPvals{$r}\t$naRanks\tNA";
		push @pctTarget, $refMeans{$r};
	} else {
		$pctBound{$r} = $refMeans{$r};
		if ($refPvals{$r} < log_limit($opt_alpha)) {
			$smed{$r}++;
		}
	}
}

my $pctTarget_min = min(@pctTarget);
if (scalar keys %pctBound > 0) {
	# other bound genes present
	push @headers, 
		sprintf("Rel. rank of other bound genes = %.4f (%u\/%u deviates from the median)",
			mean(values %pctBound), scalar (keys %smed), scalar (values %pctBound));
	my %pathIn;
	foreach my $r (sort keys %pctBound) {
		if (($pctBound{$r} < $pctTarget_min) && (exists $smed{$r})) {
			# lower than minimum cutoff and significant
			$pathOut{$r}++;
			foreach my $p (keys %key_pathways) {
				my %pg = map {$_ => 1} @{$kegg{$p}};
				if (exists $pg{$r}) {
					$pathIn{$r}++;
					delete $pathOut{$r};
				}
			}		
		}
	}
	push @headers, (
		sprintf("Genes with consistently lower relative ranks and in same pathways = %u", scalar keys %pathIn),
		sprintf("Genes with consistently lower relative ranks and in other pathways (n = %u)", scalar keys %pathOut)
	);
	
	foreach my $gene (sort keys %pathOut) {
		my @pPrimary;
		foreach my $prim (@primary) {
			my $ttest = new Statistics::TTest;
			$ttest->load_data($refRanks{$prim}, $refRanks{$gene});
			push @pPrimary, sprintf("%.4f", log_limit($ttest->{t_prob}));
			undef $ttest;
		}
		push @entries, sprintf("$gene\tOff-Target\t$pctBound{$gene}\t$refPvals{$gene}\t%s\t%s", join("\t", @pPrimary), significance($opt_alpha, @pPrimary));
	}
}

my @pHeader;
foreach my $i (@primary) {
	push @pHeader, "pDev\_$i";
}
push @headers, sprintf("Symbol\tClass\tmRelRank\tpDevMedian\t%s\tSignificance", join("\t", @pHeader));

my $f_output = $opt_path.'/'.quickname($opt_prefix, 'txt');
open my $OUTPUT_LOG, '>', $f_output or die "Can't write to device!\n";
print $OUTPUT_LOG "\# $_\n" for @headers;
print $OUTPUT_LOG "$_\n" for @entries;
close $OUTPUT_LOG;

print $f_output if (defined $opt_return);
exit;