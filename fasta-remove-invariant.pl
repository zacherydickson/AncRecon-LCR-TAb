#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Bio::SeqIO;

if(@ARGV < 1){
    die "Usage: ".basename($0). "alnIn.fasta[.gz] >| alnOut.fasta\n";
}

my $file = shift @ARGV;

my $cmd = ($file =~ /\.gz$/) ? "zcat" : "cat";

open(my $fh, "$cmd $file |") or die "Could not open $file: $!\n";
my $SeqIO = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');

my @labels; #sequence ides corresponding to rows
my @columns; #array of array ref of characters 
    #@{$columns[i]} = a column in the alignment
    #map {$_->[j]} @columns = a row in the alignment
my @variants; #array of index keyed hash ref of uniq chars
my $length = undef;
my $N = 0;

while(my $seqObj = $SeqIO->next_seq){
    push(@labels,$seqObj->id);
    my @chars = split("",uc($seqObj->seq));
    $N++;
    if(defined $length){
        unless($seqObj->length == $length){
            die "All input sequences are not the same length: @{[$seqObj->id]} differs\n"
        }
    } else {
        @columns = map {[$_]} @chars;
        @variants = map {{$_ => 1}} @chars;
        $length = $seqObj->length;
        next;
    }
    for (my $i=0; $i< @columns; $i++){
        push(@{$columns[$i]},$chars[$i]);
        $variants[$i]->{$chars[$i]} = 1;
    }
}

close($fh);

my @variantSites;
for(my $i =0; $i < @columns;$i++){
    delete $variants[$i]->{"N"};
    push(@variantSites,$i) if(scalar(keys %{$variants[$i]})  > 1);
}

@columns = @columns[@variantSites];

my $Output = Bio::SeqIO->newFh(-format => 'fasta');
for(my $j = 0; $j < $N;$j++){
    my $seq = join("",map {$_->[$j]} @columns);
    my $seqObj = Bio::Seq->new(-id => $labels[$j], seq => $seq);
    print $Output $seqObj;
}

