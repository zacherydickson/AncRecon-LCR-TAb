#!/usr/bin/perl
use warnings;
use strict;
use Bio::AlignIO;
use Bio::SeqIO;
use Class::Struct;
use File::Basename;
use File::Temp;
use Zach::Logger;
use Zach::Util::File qw(OpenFileHandle);

sub LoadCDSInfo($);
sub TestExonCount(@);
sub OutputProt(@);
sub SecCheck(@);
sub AttemptRepair(@);

struct (INTERVAL => {s => '$', e => '$'});
struct (CDSINFO => {exon_frames => '@', exon_ivs => '@', has_ifs => '$', obj => '$'});

my $MAFFT = "mafft --auto --quiet";

if(@ARGV < 1){
    die "Usage: ".basename($0). " homologouscds.fna > homologousprot.fna\n";
}

my $Logger = Zach::Logger->new(level => "WARNING");

############################################################
#Main

sub main {
    my @cdsInfoList = LoadCDSInfo(shift);
    my (undef,$tID) = split(":",$cdsInfoList[0]->obj->id);
    unless(TestExonCount(@cdsInfoList)){ #Test if all cds have the same number of segments
        $Logger->Log("$tID has variable exon numbers and is irrepairable","ERROR");
    }
    #TODO: Determine the consensus frame for each exon
    if(!grep {!$_->has_ifs} @cdsInfoList){#All prots have an in frame stop codon
        unless(SecCheck(@cdsInfoList)){
            $Logger->Log("$tID has non-frameshift stops which are not Sec: irrepairable","ERROR");
        }
    } elsif(grep {$_->has_ifs} @cdsInfoList){ #Some prots have an in frame stop
        unless(AttemptRepair(@cdsInfoList)){
            $Logger->Log("$tID has irrepairable sequences","ERROR");
        }
    } #Else no repair is needeed, just translate
    OutputProt(@cdsInfoList);
} main (shift @ARGV);


############################################################
#Subroutines

sub LoadCDSInfo($){
    my $in = Bio::SeqIO->new(-fh => OpenFileHandle(shift,"seq","ERROR"));
    my @cdsInfoList; #protID keyed hash of CDSINFO objects
    while(my $seqObj = $in->next_seq){
        my (undef,$segInfoStr) = split(/ /, $seqObj->desc);
        my @ivStrList = split(/;/,substr($segInfoStr,2,-1));
        my @exonLenList = map {(-1*eval($_)+1)} @ivStrList;
        my @exonFrameList = map {$_ % 3} @exonLenList;
        my @exonIVList;
        my $offset = 0;
        foreach my $len (@exonLenList){
            #Intervals are 1 indexed to play nice with Bio::Seq::subseq
            push(@exonIVList,INTERVAL->new(s => $offset + 1, e => $offset + $len));
            $offset += $len;
        }
        my $bInFrameStop = ($seqObj->translate()->seq =~ m/\*./) ? 1 : 0;
        push(@cdsInfoList,CDSINFO->new(exon_frames => \@exonFrameList, has_ifs => $bInFrameStop, obj => $seqObj, exon_ivs => \@exonIVList));
    }
    return @cdsInfoList;
}

#Determine plurality # of exons from cds which do not have in frame stops
sub TestExonCount(@){
    my %count;
    my $cur = undef;
    foreach my $cdsInfo (@_){
        my $nExon = scalar(@{$cdsInfo->exon_frames});
        $cur = $nExon unless(defined $cur);
        return 0 if($cur != $nExon);
    }
    return 1;
}

#Given a set of CDS info, outputs the translated sequences to stdout
sub OutputProt(@){
    my $out = Bio::SeqIO->newFh(-format => 'fasta');
    foreach my $cdsInfo (@_){
        print $out $cdsInfo->obj->translate();
    }
}


#Given a set of CDSINFO objects, all of which have an in-frame stop,
#Determines if the common stops are all selenocystein residues
#If all exons are in the same frame, and the stop codon is tga, change the A
#   to a C to recode for cysteine
#Modifies the CDSINFO objects 
sub SecCheck(@){
    my @cdsInfoList = @_;
    #Take the first set of frames as the comparison to make
    my @frameList = @{$cdsInfoList[0]->exon_frames};
    #check if any of the frames in any of the other cds's don't match
    foreach my $cdsInfo (@cdsInfoList[1 .. $#cdsInfoList]){
        for(my $i = 0; $i < @frameList; $i++){
            return 0 if($cdsInfo->exon_frames->[$i] != $frameList[$i]);
        }
    }
    #All frames are consistent, therefore the in-frame stop isn't the result of a
    #frameshift
    ##Now locate all the stops and check if they are Opal
    foreach my $cdsInfo (@cdsInfoList){
        my @AAs = split(//,$cdsInfo->obj->translate()->seq);
        for(my $i = 0; $i < @AAs; $i++){
            if($AAs[$i] eq "*"){
                my $codon = substr($cdsInfo->obj->seq,($i*3),3);
                if($codon =~ m/[tuTU][gG][aA]/){
                    #Edit the sequence so the codon is now cystein
                    my $seq = $cdsInfo->obj->seq;
                    substr($seq,($i*3)+2,1,"C");
                    $cdsInfo->obj->seq($seq);
                } else {
                    return 0;
                }
            }
        }
    }
    return 1;
}

#Man this code kinda sucks... It just keeps looping...
#Given a set of CDSINFO objects, some, but not all, of which have an in-frame stop:
#Identifies the exon in which a frame-shift may have occured
sub AttemptRepair(@){
    my @cdsInfoList = @_;
    my %deviantExonDict = #cdsInfoList index keyed hash of exon indices
    #for each cds, identify the deviant exon(s), if any
    my @cFrameList;
    for(my $i = 0; $i < @{$cdsInfoList[0]->exon_frames}; $i++){
        #Assume all cds without stops have the same exon frame
        #   find the first cds without a stop, and set its frame for
        #   this exon as the cannonical frame
        my $j = 0;
        $j++ while($cdsInfoList[$j]->has_ifs);
        my $cFrame = $cdsInfoList[$j]->exon_frames->[$i];
        push(@cFrameList,$cFrame);
        #Iterate over the cds with stops and note any whose frame differs
        #   from cannon
        for($j = 0; $j < @cdsInfoList; $j++){
            next unless($cdsInfoList[$j]->has_ifs);
            if($cdsInfoList[$j]->exon_frames->[$i] != $cFrame){
                if(exists $deviantExonDict{$j}){
                    $Logger->Log("Multiple Frameshifts identified in @{[$cdsInfoList[$j]->obj->id]}","WARNING");
                    return 0;
                }
                $deviantExonDict{$j} = $i;
            }
        }
    }
    #At this point there is 0 or 1 deviant exon for each cds with a stop
    #   zero indicates a nonsense mutation, not frame-shift
    unless(keys %deviantExonDict){ #All non-sense
        $Logger->Log("Stop codons not the result of frame-shift","WARNING");
        return 0;
    }
    #Check for any potential irrepairable nonsense
    for(my $i = 0; $i < @cdsInfoList; $i++){
        my $cdsInfo = $cdsInfoList[$i];
        next unless($cdsInfo->has_ifs);
        unless(exists $deviantExonDict{$i}){
            $Logger->Log("Irrepairable nonsense mutation for @{[$cdsInfo->obj->id]}","WARNING");
            return 0;
        }
    }
    #All ifs now have one deviant exon
    #Reverse the deviant ExonDict
    my %cdsIdxByDeviant; #exon idx keyed has of arrays of cdsInfo idx
    while(my ($cdsIdx,$exonIdx) = each %deviantExonDict){
        $cdsIdxByDeviant{$exonIdx} = [] unless(exists $cdsIdxByDeviant{$exonIdx});
        push(@{$cdsIdxByDeviant{$exonIdx}},$cdsIdx);
    }
    foreach my $j (keys %cdsIdxByDeviant){
        my @exons = map {Bio::Seq->new(-id => $_->obj->id, -seq => $_->obj->subseq($_->exon_ivs->[$j]->s,$_->exon_ivs->[$j]->e))} @cdsInfoList;
        my @valid = 0 .. $#cdsInfoList;
        splice(@valid,$_,1) foreach (sort {$b <=> $a} @{$cdsIdxByDeviant{$j}});
        foreach my $i (@{$cdsIdxByDeviant{$j}}){
            my $alnObj = GetAlignment(@exons[($i,@valid)]);
            my $validAln = $alnObj->select(2,(@valid + 1));
            my ($alignedSeq) = $alnObj->each_seq;
            my $frameShiftSize = ($cdsInfoList[$i]->exon_frames->[$j] - $cFrameList[$j]) % 3;
            my $bChangeMade = 0;
            #If an insertion caused the frameshift there will be a shared gap
            #   in all non frameshifted sequences
            #   delete the nucleotide(s) at that position
            my $gapLine = $validAln->all_gap_line;
            while($gapLine =~ m/(-+)/g){
                my $pos = pos($gapLine) -1;
                my $len = length($1);
                next unless($len % 3 == $frameShiftSize % 3);
                my $seq = $cdsInfoList[$i]->obj->seq;
                my $offset = $pos + $cdsInfoList[$i]->exon_ivs->[$j]->s;
                substr($seq,$offset,$len,"");
                $cdsInfoList[$i]->obj->seq($seq);
                if($cdsInfoList[$i]->obj->translate()->seq =~ m/\*/){
                    $Logger->Log("Deleting inserted char(s) in @{[$cdsInfoList[$i]->obj->id]} did not repair sequence","WARNING");
                    return 0;
                }
                $bChangeMade = 1;
                last;
            }
            next if($bChangeMade);
            #If a deletion caused the frameshift there will be a gap which is only in
            #   frameshifted cds
            #   insert the consensus nucleotide at that position
            my $tmpSeq = $alignedSeq->seq;
            while($tmpSeq =~ m/(\-+)/g){
                my $pos = pos($tmpSeq) - 1;
                my $len = length($1);
                next unless(-$len %3 == $frameShiftSize %3);
                my $validGapAln = $validAln->slice($pos+1,$pos+$len,1);
                next if ($validGapAln->gap_line =~ m/-/);

                my $offset = $pos + $cdsInfoList[$i]->exon_ivs->[$j]->s;
                my $seq = $cdsInfoList[$i]->obj->seq;
                my $cons = $validGapAln->consensus_string;
                substr($seq,$offset,0,$cons);
                $cdsInfoList[$i]->obj->seq($seq);
                if($cdsInfoList[$i]->obj->translate()->seq =~ m/\*/){
                    $Logger->Log("Inserting consensus char(s) in @{[$cdsInfoList[$i]->obj->id]} did not repair sequence","WARNING");
                    return 0;
                }
                $bChangeMade=1;
                last;
            }
            unless($bChangeMade){
                $Logger->Log("Could not locate an approprate frameshift for @{[$cdsInfoList[$i]->obj->id]}","WARNING");
                return 0;
            }
        }
    }
    return 1;
}

sub GetAlignment(@){
    my $tmpFile = File::Temp->new();
    my $out = Bio::SeqIO->new(-fh => $tmpFile, -format => 'fasta');
    $out->write_seq($_) foreach(@_);
    my $mafftCmd = "$MAFFT --nuc ".$tmpFile->filename;
    my $fh = OpenFileHandle("$mafftCmd |","mafftcmd","ERROR");
    my $in = Bio::AlignIO->new(-fh => $fh, -format => 'fasta');
    return $in->next_aln;
}
