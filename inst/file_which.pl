use 5.016;
use strict;
use warnings;

my $cmd = shift @ARGV;
my $sep = ':'; my @exts;
if ($^O =~ m/MSWin/) {
	$sep = ';';
	@exts = split $sep, $ENV{PATHEXT};
}
my @paths = split $sep, $ENV{PATH};

# to account for executables with/without extensions for Win32
unshift @exts, ''; 

my $target = 'NULL';
if (defined $cmd) {
	my $f = 0;
	foreach my $i (@paths) {
		for my $j (@exts) {
			my $g = $i.'/'.$cmd.$j;
			if ((-e $g) && (-x $g)) {
				$target = $g;
				$f++;
				last;
			}
		}
		last if ($f > 0);
	}
}

$target =~ s/\\/\//g;
print $target;
exit;