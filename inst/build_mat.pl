use 5.016;
use strict;
use Encode;
use XML::Simple;
use PerlIO::gzip;
use File::Type;
use File::Path 'make_path';
use Getopt::Long;
use Cwd;

my $usage = join "\n", (
	'build_mat [options] --drugnames <drug_list> --effects <sider_effects>',
	'Parameters:',
	'   --drugnames SIDER drug ID table (TSV)',
	'   --effects   SIDER MedDRA side effects (TSV)',
	'   --genelist  Table of genes and uniprot ID (TSV) to limit search scope',
	'   --drugbank  DrugBank XML file for imposing keyword search restrictions',
	'   --uniprot   DrugBank drug-target table (TSV); required for --genelist',
	'   --filters   keywords to search through DrugBank records for limiting scope',
	'   --cutoff    minimum number of drugs for which a side effect is counted',
	'   --path      path to save file output (default: current directory)',
	'   --prefix    prefix for naming output files (default: pipoft)'
);

die "Usage: $usage\n" unless (@ARGV);

my ($opt_genelist, # gene list from Durante et al
	$opt_drugbank, # xml input
	$opt_drugbank_genes, # required if $opt_genelist
	$opt_keyterms, # filter terms for drugbank data
	$opt_sider_nametable, # drug names (required)
	$opt_sider_se, # side effect data (required);
	$opt_cutoff, # minimum number of side effects found in drug cohort
	$opt_path,
	$opt_prefix, # prefix for output
	$opt_return # return filenames in STDOUT
);

GetOptions(
	'drugnames=s' => \$opt_sider_nametable,
	'effects=s'   => \$opt_sider_se,
	'genelist=s'  => \$opt_genelist,
	'drugbank=s'  => \$opt_drugbank,
	'uniprot=s'   => \$opt_drugbank_genes,
	'filters=s'   => \$opt_keyterms,
	'cutoff=i'    => \$opt_cutoff,
	'path=s'      => \$opt_path,
	'prefix=s'    => \$opt_prefix,
	'rreturn'     => \$opt_return
) or die "Usage: $usage\n";

die "SIDER compound table is not specified!\n" unless (defined $opt_sider_nametable);
die "SIDER side effect table is not specified!\n" unless (defined $opt_sider_se);

if (defined $opt_genelist) {
	die "Gene restriction is on but DrugBank Drug/Target table is not specified!\n"
		unless (defined $opt_drugbank_genes);
}

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

my $tstring = time;
if (defined $opt_path) {
	my $m = make_path($opt_path);
	$opt_path = cwd() if ($m < 1);
} else {
	$opt_path = cwd();
}
$opt_prefix = 'pipoft' unless (defined $opt_prefix); 
my $output_inhibit = $opt_path.'/'.$opt_prefix.'_filterList-'.$tstring.'.txt';
# drugNames
my $output_drugs = $opt_path.'/'.$opt_prefix.'_drugList-'.$tstring.'.txt';
#sideEffects
my $output_se = $opt_path.'/'.$opt_prefix.'_sideEffects-'.$tstring.'.txt';
# drugMatrix
my $output_mat = $opt_path.'/'.$opt_prefix.'_drugMatrix-'.$tstring.'.csv';

### PARSING XML DATA
my %inhibitors;
if (defined $opt_drugbank) {
	print STDERR "Parsing DrugBank Data...\n";
	my @filters;
	if (defined $opt_keyterms) {
		@filters = split '\,', $opt_keyterms;
	} else {
		@filters = ('inhibit');
	}
	 
	my $sep_inhibitor = "\t";
	my $IN_XML = gzInput($opt_drugbank);
	my $xml = XMLin($IN_XML);

	my $sp;
	if (scalar @filters == 1) {
		$sp = ' ';
	} else {
		$sp = "\n\t";
	}
	printf STDERR "Filtering DrugBank data for drugs with descriptions matching:%s %s\n", $sp, join(' ', @filters);
	open my $OUT_INHIBIT, '>:gzip', $output_inhibit;
	foreach my $d (keys %{$xml -> {'drug'}}) {
	    my %cur = %{$xml->{'drug'}->{$d}};
	    my @rec;
	    foreach my $item (keys %cur) {
			my $p = 0;
			foreach my $f (lc @filters) {
				if (lc $cur{$item} =~ m/$f/) {
					$p++;
					last;
				}
			}
	        if ($p > 0) {
	        	# keyterm matches - process out weird white spaces
				# non-asciis and wides
	            $cur{$item} =~ s/\t//g;
	            $cur{$item} =~ s/\n//g;
	            $cur{$item} =~ s/^\s+|\s+$//g;
				$cur{$item} =~ s/[^[:ascii:]]//g;
				$cur{$item} =~ s/[^\x00-\x7f]//g;
	            push @rec, $cur{$item};
	        }
	    }
    
	    if (scalar @rec > 0) {
	       print $OUT_INHIBIT encode('utf-8', $d.$sep_inhibitor.$rec[0]."\n");
	    }
	    undef %cur;
	    undef @rec;
		$inhibitors{$d} = lc $d;
	}
	close $OUT_INHIBIT;
}

### IF LIST OF RESTRICTED GENES IS PRESENT
# precondition: list is line delimited
my (%listRestricted, %keepDrugs); # key = UniProt, value = Symbol
if (defined $opt_genelist) {
	print STDERR "Processing gene filters...\n";
	my $IN_GENELIST = gzInput($opt_genelist);
	while (<$IN_GENELIST>) {
	    chomp;
	    my @s = split '\t', $_;
		my $gene = $s[0]; my $uniprot = $s[1];
		$listRestricted{$uniprot} = $gene;
	}
	close $IN_GENELIST;
	
	my %keepGenes; # key = UniProt, value = array of drugs
	my $IN_TARGETS = gzInput($opt_drugbank_genes);
	my $h = 1; my ($col_name, $col_uniprot);
	while (<$IN_TARGETS>) {
		chomp; my @s = split '\,', $_;
		if ($h == 1) {
			foreach my $i (0..$#s) {
				if ($s[$i] eq 'Name') {
					$col_name = $i;
				}
				if ($s[$i] eq 'UniProt ID') {
					$col_uniprot = $i;
				}
				last if (defined $col_name && defined $col_uniprot);
			}
		}
		push @{$keepGenes{$s[$col_uniprot]}}, $s[$col_name] if ($h > 1);
		$h++;
	}
	close $IN_TARGETS;
	
	foreach my $uniprot (keys %listRestricted) {
		if (exists $keepGenes{$uniprot}) {
			foreach my $drug (@{$keepGenes{$uniprot}}) {
				my $k = lc $drug;
				$keepDrugs{$k}++;
			}
		}
	}
}

### PROCESSING DRUG TABLES TO INTEGRATE GENE TARGETS
print STDERR "Loading SIDER compound table...\n";
my (%listMatches, %listDrugs, %listEffects, @mat);
my $IN_DRUG = gzInput($opt_sider_nametable);
while (<$IN_DRUG>) {
    chomp;
    my @s = split '\t', $_;
    if (exists $listDrugs{$s[0]}) {
        print STDERR "Warning: duplicated $s[0] -> $s[1] pairing!\n";
    } else {
        $listDrugs{$s[0]} = $s[1];
    }
}
close $IN_DRUG;

if ((scalar (keys %keepDrugs)) > 0) {
	# list of restricted genes present -> keep only those drugs on list
	foreach my $drug (keys %listDrugs) {
		my $k = lc $listDrugs{$drug};
		delete $listDrugs{$drug} unless (exists $keepDrugs{$k});
	}
}

my $IN_EFFECTS = gzInput($opt_sider_se);
while (<$IN_EFFECTS>) {
    chomp;
    my @s = split '\t', $_;
    if (exists $listDrugs{$s[0]}) {
		my $f = pop @s;
		my $drug = $listDrugs{$s[0]};
		$listEffects{$f}++;
		$listMatches{$drug}{$f} = 1;
	}
}
close $IN_EFFECTS;

my $n_effects = scalar (keys %listEffects);
if (defined $opt_cutoff) {
	my $max_cutoff = scalar (keys %listDrugs);
	$opt_cutoff = $max_cutoff if ($opt_cutoff > $max_cutoff);
	# check side effect counts and delete those below cutoff
	foreach my $se (keys %listEffects) {
		delete $listEffects{$se} if ($listEffects{$se} < $opt_cutoff);
	}
}

my @drugList = sort (keys %listMatches);
my @sideEffects = sort (keys %listEffects);
foreach my $d (@drugList) {
    my @row;
    foreach my $f (@sideEffects) {
        if (exists $listMatches{$d}{$f}) {
            push @row, 1;
        } else {
            push @row, 0;
        }
    }
    push @mat, join(',', @row);
	undef @row;
}

print STDERR "Generating output...\n";
open my $OUT_DRUGNAMES, '>', $output_drugs;
print $OUT_DRUGNAMES "$_\n" for @drugList;
close $OUT_DRUGNAMES;

open my $OUT_SE, '>', $output_se;
print $OUT_SE "$_\n" for @sideEffects;
close $OUT_SE;

open my $OUT_MATRIX, '>', $output_mat;
print $OUT_MATRIX "$_\n" for @mat;
close $OUT_MATRIX;

my @summary = (
	"Drug effect matrix file: $output_mat",
	"Drug names (row names): $output_drugs",
	"Side effects (columns): $output_se",
	sprintf("Dimensions: %u (row) x %u (col)", scalar @drugList, scalar @sideEffects)
);
if (scalar (keys %listRestricted) > 0) {
	push @summary, sprintf("Results are limited to %u UniProt IDs", scalar keys %listRestricted);
}
if (defined $opt_cutoff) {
	if ($opt_cutoff > 0) {
		push @summary, sprintf("Side effects removed after imposing cutoff of %u = %u", 
			$opt_cutoff, $n_effects - (scalar @sideEffects));
	}
}

### SUMMARY OUTPUT AND SEND BACK TO R
print STDERR "$_\n" for @summary;
# use system(intern = TRUE to capture and port this back to R)
if (defined $opt_return) {
	print join(',', $output_mat, $output_drugs, $output_se);
}

exit;