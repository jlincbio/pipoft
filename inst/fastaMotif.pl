use 5.016;
use strict;
use warnings;
use Scalar::Util 'looks_like_number';
use PerlIO::gzip;
use File::Type;

# added $bed_output to allow R to designate temp files

## Version sorting via the alphanum algorithm
## Reference: http://www.DaveKoelle.com

## Note: this algorithm will give different results than GNU sort
##.      (invoked by sort -V) if the input begins with non-alphanumeric
##		 letters, e.g. literal ^, >, # characters, in which GNU sort
##		 does not sort those ahead of alphanumerics despite their
##		 ASCII codes being numerically smaller.

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

sub natsort {
	my @sorted;
	if (@_) {
		if ((scalar @_ == 1) && (ref $_[0] eq 'ARRAY')) {
			# array reference provided;
			my $array = shift;
			@sorted = sort {natsort_cmpfun($a,$b)} @{$array};
			return \@sorted;
		} else {
			# assume array passed to natsort
			@sorted = sort {natsort_cmpfun($a,$b)} @_;
		}
	}
	return @sorted;
}

sub natsort_cmpfun {
	# split strings into chunks
	my @a = natsort_split($_[0]);
	my @b = natsort_split($_[1]);
	
	while (@a && @b) {
		my $a_chunk = shift @a;
		my $b_chunk = shift @b;
		# comparison test: if $a and $b numeric, treat as so;
		# otherwise compare as strings, return result if unequal
		my $test =
			(($a_chunk =~ /\d/) && ($b_chunk =~ /\d/)) ? 
			$a_chunk <=> $b_chunk :
			$a_chunk cmp $b_chunk ;  
		return $test if ($test != 0);
	}
	return @a <=> @b; # return longer string.
}

# splitting function for numeric/non-numeric transitions
sub natsort_split {
	# split conditions:
	# zero width, digit preceded by non-digit or otherwise
	my @chunks = split m{
		(?=
		(?<=\D)\d |
		(?<=\d)\D)
	}x, $_[0];
	return @chunks;
}

my %alphabet_expansion = (
	'A' => 'A',
	'C' => 'C',
	'G' => 'G',
	'T' => 'T',
	'M' => '{A,C}',
	'R' => '{A,G}',
	'W' => '{A,T}',
	'S' => '{C,G}',
	'Y' => '{C,T}',
	'K' => '{G,T}',
	'V' => '{A,C,G}',
	'H' => '{A,C,T}',
	'D' => '{A,G,T}',
	'B' => '{C,G,T}',
	'N' => '{A,C,G,T}'
);

sub readSeq {
   my ($FH_SEQ, $fasta) = @_; my $file_not_empty = 0; 
   $fasta->{seq} = undef; # clear out previous sequence and put header in place
   $fasta->{header} = $fasta->{next_header} if ($fasta->{next_header});
   while (<$FH_SEQ>) {
      $file_not_empty = 1;
      next if /^\s*$/;  # skip blank lines
      chomp;    

      if (/^>/) { # fasta header line
	  	(my $h = $_) =~ s/^>//;
		if ($fasta->{header}) {
			$fasta->{next_header} = $h;
			return $fasta;
		} else { # first time through only
			$fasta->{header} = $h;
		}              
      } else {
		  s/\s+//;  # remove any white space
		  $fasta->{seq} .= $_;
      }         
   }    
   return $fasta if ($file_not_empty);
   $fasta->{header} = $fasta->{seq} = $fasta->{next_header} = undef; # cleanup
   return;       
}
sub letterExpand {
	my $motif = shift; my $output;
	my @m = split '', $motif;
	foreach my $i (@m) {
		$output .= $alphabet_expansion{$i};
	}
	return $output;
}
sub randomString {
	my $n = shift;
	my @ascii = (48..57, 65..90, 97..122); # ASCII bank
	my $string; $string .= chr($ascii[rand @ascii]) for 1..$n;
	return $string;
}
sub check_gzip {
	my ($file, $mode) = @_;
	my $test = File::Type->new->checktype_filename($file);
	return "$mode\:gzip" if ($test =~ m/gzip/);
	return $mode;
}

sub stringToMotif {
	my $hashRef = shift; my %bankHash = %{$hashRef};
	my %expMotif;
	foreach my $inMotif (keys %bankHash) {
		my $fe_motif = $bankHash{$inMotif};
		$inMotif = uc $inMotif; my $rcMotif = reverse $inMotif;
		$rcMotif =~ tr/ACGTRYKMBVDH/TGCAYRMKVBHD/;
		$inMotif = letterExpand($inMotif);
		$rcMotif = letterExpand($rcMotif);
		my @motifs = (glob($inMotif), glob($rcMotif));
		foreach my $key (@motifs) {
			$expMotif{$key} = $fe_motif;
		}
	}
	return \%expMotif;
}
sub fastaToMotifBed {
	my ($fasta_filename, $ref_motif, $output_motifBed, $allow_blockwrites, $only_counts) = @_;
	my (%seqbank, %motifBank, @array_motifKey);
	if (ref($ref_motif) eq 'ARRAY') {
		@array_motifKey = @{$ref_motif};
	} elsif (ref($ref_motif) eq 'HASH'){
		@array_motifKey = keys %{$ref_motif};
	} else {
		die "Error: improper motif entry detected!\n";
	}
	
	# $allow_blockwrites: integer indicate how many entries to save
	my $BED_TEMP; my $lineCount = 0;
	if (defined $allow_blockwrites) {
		my $output_tempBed = '.'.randomString(12).'.bed';
		open $BED_TEMP, '>', $output_tempBed;
	}
	
	# $only_counts: just add number without preparing bed
	my $mode_input = check_gzip($fasta_filename, '<');
	open my $FASTA_INPUT, $mode_input, $fasta_filename or die "Error: cannot open $fasta_filename!\n";
	while (readSeq($FASTA_INPUT, \%seqbank)) {
		foreach my $key (@array_motifKey) {
			if (defined $allow_blockwrites) {
				if ($lineCount >= $allow_blockwrites) {
					print $BED_TEMP "$_\n" for keys %motifBank;
					$lineCount = 0;
				}
			}
			(my $coordinate = $seqbank{header}) =~ s/\:/\-/g;
			my @header = split /\-/, $coordinate;
			# process header; proper format should be chrN:n1-n2, otherwise reset
			$header[0] = $seqbank{header} unless (scalar @header == 3);
			$header[1] = 0 unless ((defined $header[1]) && looks_like_number($header[1]));
			my $sequence = uc($seqbank{seq});
			my $offset = 0; my @keyBed;
			my $pos = index($sequence, $key, $offset);
			while($pos != -1) {
				my $l = sprintf "%s\t%u\t%u\t%s", $header[0], $header[1] + $pos, $header[1] + $pos + length($key), $key;
				$motifBank{$l}++; $lineCount++;
				$offset = $pos + 1; $pos = index($sequence, $key, $offset);
			}			
		}
	}
	close $FASTA_INPUT;
	
	my @motifBed = natsort(keys %motifBank);
	if (defined $output_motifBed) {
		open my $BED_OUTPUT, '>', $output_motifBed or die "Error: cannot open $output_motifBed\n";
		print $BED_OUTPUT "$_\n" for @motifBed;
		return 1;
	} else {
		return \@motifBed;
	}
}

my $usage = "perl fastaMotif.pl <fasta file> <motif> <output_file>";
my ($fasta_filename, $inMotif, $bed_output) = @ARGV;
die "Usage: $usage\n" unless (defined $fasta_filename && defined $inMotif);

# convert input to hash for backward compatibility with unreleased code
my %hashMotif; $hashMotif{$inMotif} = 1;
$inMotif = \%hashMotif;

my $arrayMotif = stringToMotif($inMotif);
my $outputBed = fastaToMotifBed($fasta_filename, $arrayMotif);

open my $OUTPUT, '>', $bed_output or die "Can't write to $bed_output!\n";
print $OUTPUT "$_\n" for @{$outputBed};
close $OUTPUT;
