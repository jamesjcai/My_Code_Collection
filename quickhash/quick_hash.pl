#!/usr/bin/perl
use File::Temp qw/ tempfile tempdir /;
use IO::String;

open FILE1, "<input1.txt" or die $!;
open FILE2, "<input2.txt" or die $!;
open FILE3, ">output.txt" or die $!;



$novalue=0;
my @LinesA = <FILE1>;
my @LinesB = <FILE2>;

	foreach $lineB (@LinesB) {
		my ($genename, @desc) = split ("\t", _trim($lineB));
		#if (!(@desc)) {$desc[0]=1; $novalue=1;}
		$gene_dic{$genename} = join("\t", @desc);
		#$gene_dic{$genename} = join(",",$gene_dic{$genename}, join("\t", @desc));
	}


	foreach $lineA (@LinesA) {
		my ($genename) = split ("\t", _trim($lineA));
		print FILE3 join("\t", $genename, $gene_dic{$genename}), "\n";

	}


close(FILE1);
close(FILE2);
close(FILE3);


#############
### SUBS  ###
#############

sub _trim {
	my($v) = @_;
	$v=~s/^\s+//;
	$v=~s/\s+$//;
	return $v;
}


sub _improveSecurity {
    $_ = shift;
    s/\-+(.*)/\1/g;
    s/(.*)[ \t]+\-(.*)/\1\2/g;
#    tr/\$\'\`\"\<\>\/\;\!\|/_/;	
    tr/\$\'\`\"\<\/\;\!\|/_/;	
    return($_);
}
