use 5.016;
use strict;
use warnings;
use File::Fetch;
use File::Temp 'tempdir';

$File::Fetch::WARN = 0;

# keggColorMaps.pl
# download the colored pathway maps from KEGG using KEGG Get API
# API manual: http://www.kegg.jp/kegg/docs/weblink.html)
# Working as of Aug 2018

my $pngHandle = 'img src'; # HTML tag for png link
my ($pathway, $input, $output) = @ARGV;
unless ((defined $pathway) && (defined $input)) {
	die "Error: pathway input not specified!\n";
}

my $base = 'http://www.kegg.jp/kegg-bin/show_pathway?'.$pathway;

my @cols; my $err = 0;
print "Generating color maps...\n";
open my $INPUT, '<', $input or die "Error: cannot read from $input!\n";
while (<$INPUT>) {
	chomp;
	next if ($_ =~ m/^\#/);
	my @entry = split '\t', $_;
	$entry[1] =~ s/\#/\%23/g;
	my $color = $entry[0].'%09'.$entry[1];
	push @cols, $color;
}

my $keggUrl = join '/', $base, @cols;

print "Connecting to KEGG...\n";
my $tempPath = tempdir(CLEANUP => 1);
my $keggFetch = File::Fetch->new(uri => $keggUrl);
my $keggLink = $keggFetch->fetch(to => $tempPath);
open my $KEGGTEMP, '<', $keggLink;
my ($pngLink, @pngRef);
while (<$KEGGTEMP>) {
	chomp;
	next unless ($_ =~ m/png/);
	push @pngRef, $_;
}
close $KEGGTEMP;
foreach my $i (@pngRef) {
	if (defined $pngLink) {
		last;
	} else {
		my @j = split '\"', $i;
		foreach my $k (0..$#j) {
			if ($j[$k] =~ m/$pngHandle/) {
				# next term is $pngLink
				$pngLink = $j[$k + 1];
				last;
			}
		}
	}
}

if (defined $pngLink) {
	print "Fetching pathway image...\n";
	my @pterms = split '\/', $pngLink;
	my $pname = pop @pterms;
	$pngLink = 'http://www.kegg.jp/'.$pngLink;
	my $pngFetch = File::Fetch->new(uri => $pngLink);
	$pngFetch->fetch() or die $pngFetch->error;
} else {
	print "Warning: cannot initiate image D/L. Paste this URL to your browser instead:\n\n";
	print "$keggUrl\n";
	$err = 2;
}

exit $err;