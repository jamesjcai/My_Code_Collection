#!/usr/bin/perl  

my $refseq =  $ARGV[0];  # gene annotation - UCSC
my $output = $ARGV[1]; # length of sequence read

open(REF, $refseq) || die("Could not open file!");
#open(SAM, $sam) || die("Could not open file!");
open (RRR, ">$output");
#my %genome = ();

my %knownGene = ();
my %infor = ();
my $line =  <REF>;
chomp($line);
while (<REF>) {
  chomp($_);
  my @transcript = split(/\t/);
  my $name  = $transcript[12];
  my $tran  = $transcript[1]; 
  if($knownGene{$name}{$tran} == NULL )
  {$knownGene{$name}{$tran} = 1; }
  else
  {
    $knownGene{$name}{$tran} =  $knownGene{$name}{$tran} +1; 
  }
   $infor{$name}{$tran}{"chrom"} = $transcript[2];
   $infor{$name}{$tran}{"strand"} = $transcript[3];
   $infor{$name}{$tran}{"txStart"} = $transcript[4];
   $infor{$name}{$tran}{"txEnd"} = $transcript[5];
   $infor{$name}{$tran}{"exonCount"} = $transcript[8];
   $infor{$name}{$tran}{"exonStarts"} = $transcript[9];
   $infor{$name}{$tran}{"exonEnds"} = $transcript[10];
   $infor{$name}{$tran}{"cds_s"} = $transcript[13];
   $infor{$name}{$tran}{"cds_e"} = $transcript[14];

}

    foreach my $name  (keys %knownGene )
    {   
        my $index = 0;
        foreach my $tran (keys %{$knownGene{$name}})
        {
            if(  $knownGene{$name}{$tran} > 1 ) 
            {$index= 1;}
        }
        if($index == 0) 
        {
            foreach my $tran (keys %{$knownGene{$name}})
           {

             my $chrom = $infor{$name}{$tran}{"chrom"};
             my $strand = $infor{$name}{$tran}{"strand"};
             my $ts = $infor{$name}{$tran}{"txStart"};
             my $te = $infor{$name}{$tran}{"txEnd"};
             my $ec = $infor{$name}{$tran}{"exonCount"};
             my $es = $infor{$name}{$tran}{"exonStarts"};
             my $ee = $infor{$name}{$tran}{"exonEnds"};
             print RRR "$name\t$tran\t$chrom\t$strand\t$ts\t$te\t$ec\t$es\t$ee\n";
           }
        }
    }
close(RRR);
close(REF);




