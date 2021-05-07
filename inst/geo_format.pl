use 5.016;
use strict;
use warnings;
use PerlIO::gzip;
use File::Type;
use List::Util 'uniq';

sub trunc {
	my $x = shift;
	$x =~ s/\/\/\//\/\//g;
	my @g = split " \/\/ ", $x; my @y;
	foreach my $i (@g) {
		next if ($i =~ m/^N[MR]\_/);
		next if ($i =~ m/^X[MR]\_/);
		next if ($i =~ m/^ENS[TG][0-9]/);
		next if ($i =~ m/^BC0[0-9]/);
		push @y, $i;
	}
	@y = uniq @y;
	return join('|', @y);
}

my ($f_input, $f_output, $gzip_output) = @ARGV;
die "No input file specified!\n" unless (defined $f_input);
unless (defined $f_output) {
	$f_output = $f_input.'.processed';
	print "Output not specified; defaults to:\n$f_output\n";
}

my ($fid, $fkey, $FINPUT, %res);
my $t_input = File::Type->new->checktype_filename($f_input);
if ($t_input =~ m/gzip/) {
	open $FINPUT, '<:gzip', $f_input or die "Error opening $f_input\!\n";
} else {
	open $FINPUT, '<', $f_input or die "Error opening $f_input\!\n";
}

while (<$FINPUT>) {
	next if ($_ =~ m/[\!\#\^]/);
	next if ($_ =~ m/^\s/);
	my $s = $_; chomp $s;
	my @l = split '\t', $s;
	if ((defined $fid) && (defined $fkey)) {
		my $tid = $l[$fid];
		my $tx = $l[$fkey];
		if (defined $tx) {
			if ($tx ne '') {
				if ($tx ne '---') {
					unless ($tx =~ m/\s/){
						if ($tx =~ m/\/\//) {
							my $y = trunc($tx);
							$res{$tid} = $y;
						} else {
							$res{$tid} = $tx;
						}
					}
				}
			}
		}
	} else {
		print "Multiple fields found:\n";
		foreach my $i (0..$#l) {
			print "$i -> $l[$i]\n";
		}
		printf "Select ID field [0 - %u]: ", $#l;
		$fid = <STDIN>; chomp $fid;
		printf "Select Symbol field [0 - %u]: ", $#l;
		$fkey = <STDIN>; chomp $fkey;
	}
}
close $FINPUT;

my $n = 0; my $FOUTPUT;
if ((defined $gzip_output) && ($gzip_output > 0)) {
open $FOUTPUT, '>:gzip', $f_output or die "Error writing to $f_output\!\n";
} else {
	open $FOUTPUT, '>', $f_output or die "Error writing to $f_output\!\n";
}
foreach my $i (keys %res) {
	# last if ($n > $N_DEBUG);
	if (defined ($res{$i})) {
		print $FOUTPUT "$i\t$res{$i}\n";
		#$n++;
	}
}
close $FOUTPUT and do {
	print "Converted expression data:\n$f_output\n";
	exit 0;
};

exit 1;