use strict;
use warnings;
use Storable;

my ($ref_pathway, $ref_genes);
my ($input_pathway, $input_genes, $f_temp, $pretty) = @ARGV;
die "Error: input data files not specified!\n" unless ((defined $input_pathway) && (defined $input_genes));
print STDERR "Reading library files...\n";
if ((-e $input_pathway) && (-e $input_genes)) {
	$ref_pathway = retrieve $input_pathway;
	$ref_genes = retrieve $input_genes;
} else {
	print STDERR "Error: KEGG library pathway data cannot be found or is incomplete!\n";
	exit 1;
}
my %info = %{$ref_pathway};
my %kegg = %{$ref_genes};
my %acc = map {$_ => 1} (keys %info, keys %kegg);
my %r1;

my $err_input = 0;
foreach my $k (keys %acc) {
	if ((exists $info{$k}) && (exists $kegg{$k})) {
		# mismatch check
		if (ref($info{$k}) eq 'ARRAY') {
			print STDERR "Error: gene list detected instead of pathway information!\n";
			$err_input++;
			last;
		}
		if (ref($kegg{$k}) ne 'ARRAY') {
			print STDERR "Error: pathway information detected instead of gene list!\n";
			$err_input++;
			last;
		}		
		if ($pretty > 1) {
			$info{$k} =~ s/\s-\s[A-Z].*\s\(.*\)//;
		}
		my %genes = map {$_ => 1} @{$kegg{$k}};
		my $gs = join '::', sort (keys %genes);
		$r1{$k} = "$k\t$gs\t\'$info{$k}\'";		
	} else {
		print STDERR "Skipping $k due to missing records in either input\n";
	}
}
die "Error: possible input mismatch!\n" if ($err_input > 0);

open my $F_OUTPUT, '>', $f_temp or die "Error! Can't create temporary file at location $f_temp!\n";
print $F_OUTPUT "$_\n" for (sort values %r1);
close $F_OUTPUT;

exit 0;