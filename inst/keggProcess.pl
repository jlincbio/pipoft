use 5.016;
use strict;
use warnings;
use File::Fetch;

# generate a list of symbols found in a KEGG pathway from URL
# http://rest.kegg.jp/get/hsaXXXXX (XXXXX = KEGG pathway ID)
# modified for R compatibility

$File::Fetch::WARN = 0;

my $pathwayName = shift @ARGV or die "No pathway specified! Usage: $0 <pathway_code>\n";
$pathwayName = lc($pathwayName);
my $kegg_fetch = File::Fetch->new(uri => "http://rest.kegg.jp/get/$pathwayName");
my $kegg_file = $kegg_fetch->fetch() or die $kegg_fetch->error;

(my $kegg_code = $pathwayName) =~ s/\d//g;

open my $PATHWAY, '<', $kegg_file;
my $l = 1; my $geneSection;
my @geneNames;
my @geneCodes;

while (<$PATHWAY>) {
	chomp;
	if (m/^GENE/) {
		$geneSection = $l;
	} 
	
	next if (not defined $geneSection);
	my @lineItems = split ' ';
	if ($l == $geneSection) {
		shift @lineItems;
	}
	
	$l++;
	
	if ($lineItems[0] =~ m/^\d/) {
		$lineItems[1] =~ s/;//;
		push @geneNames, $lineItems[1];
		push @geneCodes, $lineItems[0];
	} else {
		last;
	}
}

open my $GENELIST, '>', "pathway_$pathwayName.txt" or die "General I/O Error!\n";
open my $COLORMAP, '>', "colorMap_$pathwayName.txt" or die "General I/O Error!\n";
open my $KEGGCODE, '>', "keggKey_$pathwayName.txt" or die "General I/O Error!\n";

print $COLORMAP "\#$kegg_code\tFG\n";
print $KEGGCODE "\#$kegg_code\tSymbol\n";
foreach my $i (0..$#geneNames) {
	print $COLORMAP "$geneCodes[$i]\t\#84AF84\,\#E0E0E0\n";
	print $KEGGCODE "$geneCodes[$i]\t$geneNames[$i]\n";
	print $GENELIST "$geneNames[$i]\n";
}

close $GENELIST;
close $COLORMAP;
close $KEGGCODE;

print "Gene list:  pathway_$pathwayName.txt\n";
print "Accession:  keggKey_$pathwayName.txt\n";
print "KEGGColor: colorMap_$pathwayName.txt\n";
unlink $pathwayName;