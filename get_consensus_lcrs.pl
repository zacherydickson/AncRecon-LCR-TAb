#!/usr/bin/perl

##Given a set of homologous protein sequences:
#   For each sequence identify the minimum entropy of a window in the protein, then run segA
#   using this as
#   the K1, and K2 some constant value higher than that, note the location of all LCR intervals
#   Generate an msa, using mafft, for the entire set of proteins
#   Generate a profile (including gaps) for the msa
#   Convert all sequence specific LCR intervals to alignment coordinates
#   Merge any overlapping intervals
#   For each protein calculate the length (# of non-gap positions) in the consensus LCR intervals
#   Calculate the entropy of the consensus sequence (excluding gaps) in the LCR intervals,
#   this can be used to decide on one interval per protein (i.e. the minimum entropy)
#(Don't remember why I thought I needed this) ####   The outputs would be ConsensusAln in CIGAR format (maybe in header)
#       Interval info (maybe in the header) Coords/ID, entropy
#       Protein info:  Sample, MinEntropy, CIGAR, Len1;Len2
#           The Lens would be in the same order as the interval info, which themselves are
#           ordered by increasing entropy
#           MinEntropy can be used to classify into LCR+/- based on a lowcut threshold
#   ##BEGIN/END ENTRY lines to separate information for each protein

use warnings;
use strict;
use Bio::SeqIO;
use Bio::AlignIO;
use Class::Struct;
use File::Basename;
use Getopt::Std;
use Zach::Logger;
use Zach::Util::Array qw(Unique Unzip);
use Zach::Util::File qw(OpenFileHandle);
use Zach::Util::Input qw(ProcessNumericOption);

############################################################
#Declarations

struct (INTERVAL => {start => '$', end => '$'});
struct (PROTINFO => {h => '$', regions => '$'});

sub CalcEntropy(@);
sub FindMinEntropy($);
sub RunSeg($$);
sub GetRegions($);
sub GetAln($);
sub GetCIGAR($);
sub GetConsensusRegions(\%\%);
sub CalcConsensusEntropy(\@$);
sub CalcConsensusRegionLen(\@$);


my @_AA = qw(A C D E F G H I K L M N P Q R S T V W Y);
my %_AAIDX = map {($_AA[$_] => [$_])} (0 .. $#_AA);
$_AAIDX{U} = $_AAIDX{C};
$_AAIDX{X} = [0 .. $#_AA];
my $_SEGA = "/home/zac/bin/segA";
my $_MAFFT = "mafft --auto --quiet --amino";
my %DEFAULT = (k => 0.3, v => 1, w => 15);

############################################################
#Handling Input

my %Opts;
getopts('i:k:w:v:h',\%Opts);

if(@ARGV < 1 or exists $Opts{h}){
    my $usage = basename($0). " [-options] Fasta[.gz] > output.creg";
    die "Usage: $usage\n\tUse -h for more info\n" unless(exists $Opts{h});
    die "===Description\n".
        "\tGiven a set of homologous protein sequences, identifies consensus LCRs:\n".
        "\t\tRegions within the consensus sequence which have minimum entropy in any\n".
        "\t\tof the input sequences\n".
        "===Usage\n\t$usage\n".
        "===Input\n".
        "\tFasta[.gz]\tA fasta formatted file containing proteins sequences\n".
        "===Output\n".
        "\tA irregular file format describing the consensus lcrs:\n".
        "\t\t#ID XXXX\n".
        #        "\t\t#CCIGAR XXXX\n".
        "\t\t#REG start end entropy\n".
        "\t\tsample\tmin_entropy\tCIGAR\tLen1;...\n".
        "===Options\n".
        "\t-i str\t[parsed from filename] An identifier for the protein\n". 
        "\t-k [0.0-∞)\t[$DEFAULT{k}] A value describing how much higher than\n".
        "\t\tthe minimum entropy The K2 threshold is\n".
        "\t-w [1-∞]\t[$DEFAULT{w}] An integer describing the length of sliding window to use\n".
        "\t-v [0-4]\t[$DEFAULT{v}] Verbosity: 0-Errors,1-Warnings,2-Info,3-Debug\n".
        "===Flags\n".
        "\t-h\tDisplay this message and exit\n";
} 

$Opts{k} = ProcessNumericOption($Opts{k},$DEFAULT{k},0,undef,0,"K1");
$Opts{w} = ProcessNumericOption($Opts{w},$DEFAULT{w},1,undef,1,"window");
$Opts{v} = ProcessNumericOption($Opts{v},$DEFAULT{v},0,4,1,"Verbosity");

my $Logger = Zach::Logger->new(level => $Zach::Logger::LEVEL_LOG{$Opts{v}});

############################################################
#Main

sub main {
    my $file = shift;
    #Parse the input file for a protein id if necessary
    $Opts{i} = basename($file,'.gz','.fna','.fa','.fasta') unless(exists $Opts{i});
    $Logger->title($Opts{i});
    my $alnObj = GetAln($file);
    my %regions = GetRegions($file);
    
    my %cigars = GetCIGAR($alnObj);
    my @cRegList = GetConsensusRegions(%regions,%cigars);
    my %cRegHListDict = CalcConsensusEntropy(@cRegList,$alnObj);
    #Order the consensus regions by ascending entropy
    my @order = sort {$cRegHListDict{_Consensus_}->[$a] <=> $cRegHListDict{_Consensus_}->[$b]} (0 .. $#cRegList);
    @cRegList = @cRegList[@order];
    foreach my $sID (keys %cRegHListDict){
        @{$cRegHListDict{$sID}} = @{$cRegHListDict{$sID}}[@order];
    }
    my %cRegLenListDict = CalcConsensusRegionLen(@cRegList,$alnObj);
    #Output the results
    print "#ID $Opts{i}\n";
    for(my $i = 0; $i < @cRegList; $i++){
        printf("#REG %d %d %0.3f\n",$cRegList[$i]->start,$cRegList[$i]->end,$cRegHListDict{_Consensus_}->[$i]);
    }
    foreach my $sID (sort keys %regions){
        my $lenStr = join(";",@{$cRegLenListDict{$sID}});
        my @hList = @{$cRegHListDict{$sID}};
        for(my $i = 0; $i < @hList; $i++){
            if(defined $hList[$i]){
                $hList[$i] = sprintf("%0.3f",$hList[$i]);
            } else {
                $hList[$i] = "NA";
            }
        }
        my $hStr = join(";",@hList);
        printf("%s\t%0.3f\t%s\t%s\t%s\n",$sID,$regions{$sID}->h,$cigars{$sID},$lenStr,$hStr);
    }
} main(shift(@ARGV));

############################################################
##Subroutines

#Calculates the shannon entropy for a vector of residue counts
#If the vector sum is 0, undefined is returned
sub CalcEntropy(@){
    my @sv = @_;
    my $w = eval(join("+",@sv));
    return undef unless($w);
    my $h = 0;
    foreach my $x (@sv){
        if($x < -0.001){
            $Logger->Log("Encountered negative count when Calculating Entropy","ERROR");
        }
        next unless($x > 0);
        my $p = $x/$w;
        $h -= $p * log($p)
    }
    return ($h / log(2));
}

#Given a protein sequence, calculates the minimum entropy of a window within that sequence
#Inputs - a protein sequence
#Ouput - an entropy value
sub FindMinEntropy($){
    my $seq = shift;
    $Logger->Log("Finding minimum entropy in sequence ....","DEBUG");
    my $min = 100;
    my @sv = (0) x @_AA;
    my @chars = map {(exists $_AAIDX{$_}) ? $_ : 'X'} split("",$seq);
    my $i = 0;
    $Logger->Log("filling window","DEBUG");
    #$sv[$_AAIDX{$chars[$i++]}]++ while($i < @chars and $i < $opts{w});
    while($i < @chars and $i < $Opts{w}){
        my @idx = @{$_AAIDX{$chars[$i]}};
        @sv[@idx] = map {$_ + 1/@idx} @sv[@idx];
    } continue{$i++}
    $Logger->Log("sliding window ($i)","DEBUG");
    #i will either be the length of the sequence or the length of the window at this point
    do{
        #$Logger->Log("$i: ".join(",",map {sprintf("%0.3f",$_)} @sv),"DEBUG");
        my $h = CalcEntropy(@sv);
        return undef unless(defined $h);
        $min = $h if($h < $min);
        if($i < @chars){
            my @posidx = @{$_AAIDX{$chars[$i]}};
            my @negidx = @{$_AAIDX{$chars[$i-$Opts{w}]}};
            @sv[@posidx] = map {$_ + 1/@posidx} @sv[@posidx];
            @sv[@negidx] = map {$_ - 1/@negidx} @sv[@negidx];
        }
    } while(++$i <= @chars);
    $Logger->Log("Minimum entropy of $min calculated","DEBUG");
    return $min;
}


#Given a sequence and the lowest entropy of a window in that sequence, run seg with that
#minimum entropy as K1; returns the lowest entropy 'lcr' found by seg with those parameters
#Inputs - A Bio::Seq Object
#       - the minimum entrpy of a window in the the sequence
#Output - An array of INTERVAL objects corresponding to each LCR
#           Note: segA intervals are 1 indexed
sub RunSeg($$){
    my ($seqObj, $minH) = @_;
    $Logger->Log("Running segA on @{[$seqObj->id]}...","DEBUG");
    my $tmpFile = File::Temp->new(TMPDIR => "/tmp", suffix => ".fna", TEMPLATE => "minHReg_XXXXXX",
        UNLINK => 1);
    $Logger->Log("tmpSequence in $tmpFile","DEBUG");
    my $oObj = Bio::SeqIO->newFh(-fh => $tmpFile, -format => 'fasta');
    print $oObj $seqObj;
    my $K1 = (int($minH * 1000) + 1) / 1000;
    my $K2 = $K1 + $Opts{k};
    my $segACmd = "$_SEGA @{[$tmpFile->filename]} $Opts{w} $K1 $K2 -l";
    $Logger->Log("Running segA with: $segACmd","DEBUG");
    my @regObjList;
    if(open(my $fh, "$segACmd |")){
        my $iObj = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
        while(my $lcrSeqObj = $iObj->next_seq){
            $Logger->Log("seg gave @{[$lcrSeqObj->id]} @{[$lcrSeqObj->desc]}","DEBUG");
            my ($id,$start,$end) = split(/\(|-|\)/,$lcrSeqObj->id);
            my (undef,$h) = split(/=| /,$lcrSeqObj->desc);
            $Logger->Log("parsed: $id,$start,$end,$h","DEBUG");
            my $regObj = INTERVAL->new(start => $start, end => $end);
            push(@regObjList,$regObj);
        }
    } else {
        $Logger->Log($seqObj->id." segA error: $!","WARNING");
    }
    $Logger->Log("seqA found @{[scalar(@regObjList)]} LCRs in @{[$seqObj->id]}","DEBUG");
    return @regObjList;
}

#Given a path to a fasta formatted file of protein sequences
#   runs mafft to get an MSA
#Inputs - A path to a fasta[.gz] file
#Output - A Bio::SimpleAln Object;
sub GetAln($){
    my $file = shift;
    $Logger->Log("Aligning sequences in $file ...","INFO");
    my $catCmd = ($file =~ /\.gz$/) ? "zcat" : "cat";
    my $mafftCmd = "$catCmd $file | $_MAFFT - |";
    $Logger->Log("Running mafft with cmd: $mafftCmd","DEBUG");
    open(my $fh, $mafftCmd) or $Logger->Log("Could not run mafft on $file: $!","ERROR");
    my $in = Bio::AlignIO->new(-fh => $fh, -format => 'fasta');
    my $alnObj = $in->next_aln();
    $Logger->Log("Alignment complete","DEBUG");
    return $alnObj;
}

#Given a path to a fasta fromatted file of homologous protein sequences
#   Find the minimum entropy regions in each protein sequence
#Inputs - A path the a fasta[.gz] file
#Output - An sID keyed Hash of PROTINFO Objects
sub GetRegions($){
    $Logger->Log("Identifying LCRs ...","INFO");
    my $fh = OpenFileHandle(shift,"fasta","ERROR");
    my $in = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
    my %sInfo;
    while(my $seqObj = $in->next_seq){
        my ($sID) = split(/:/,$seqObj->id);
        my $h = FindMinEntropy($seqObj->seq);
        my @ivList = RunSeg($seqObj,$h);
        my $infoObj = PROTINFO->new(h => $h, regions => \@ivList);
        $sInfo{$sID} = $infoObj;
    }
    $Logger->Log("LCRs Identified","DEBUG");
    return %sInfo;
}

#Given an alignment, generates a CIGAR string for each sequence
#The cigar strings with contain only D and M entries, which assumes the consensus
#   has a residue at every position
#Inputs - A Bio::SimpleAln Object
#Outputs - An sID keyed hash of cigar strings
sub GetCIGAR($) {
    my $alnObj = shift;
    $Logger->Log("Generating CIGAR strings from alignment ...","INFO");
    my %cigars;
    foreach my $seqObj ($alnObj->each_seq){
        my ($sID) = split(/:/,$seqObj->id);
        my @state = map {($_ eq $alnObj->gap_char)} split('',$seqObj->seq);
        my $curState = $state[0];
        my $count = 1;
        my $str;
        for(my $i = 1; $i < @state; $i++){
            if($curState == $state[$i]){
                $count++;
                next;
            }
            $str .= $count . (($curState) ? 'D' : 'M');
            $curState = $state[$i];
            $count = 1;
        }
        $str .= $count . (($curState) ? 'D' : 'M');
        $cigars{$sID} = $str;
    }
    $Logger->Log("CIGAR strings generated","DEBUG");
    return %cigars;
}

#Given an interval in protein coordinates, and a cigars string describing that proteins 
#   alignment to the consensus, convert to consensus coordinates
#Inputs - an INTERVAL Object
#       - a CIGAR string
#Output - an INTERVAL object
#           If the cigar string cannot convert the coords (a cigar for the the wrong prot)
#               undef is returned
sub ConvertCoords($$){
    my ($protIV, $cigar) = @_;
    unless(defined $protIV->start and defined $protIV->end and defined $cigar){
        $Logger->Log("Attempt to convert coords with undefined location or cigar","WARNING");
        return undef;
    }
    if($protIV->start < 1){
        $Logger->Log("Attempt to convert coords with interval which begins before the alignment","WARNING");
        return undef;
    }
    if($protIV->start > $protIV->end){
        $Logger->Log("Attempt to convert coords with an improper interval","WARNING");
        return undef;
    }
    $Logger->Log("Converting @{[$protIV->start,$protIV->end]} with $cigar","DEBUG");
    my $curPos = 1;
    my $conIV = INTERVAL->new(start => $protIV->start, end => $protIV->end);
    my @stateLen = split(/[MD]/,$cigar);
    my @stateTotals = (0,0);
    $stateTotals[0] = eval(join("+",(Unzip(@stateLen))[0 .. int(@stateLen/2)]));
    $stateTotals[1] = eval(join("+",(Unzip(@stateLen))[int(@stateLen/2) + 1 .. $#stateLen]));
    my $protLen = $stateTotals[1]; #assume starting on a D
    #Adjust start position if the alnignment starts with a match
    $cigar =~ m/([MD])/;
    if($1 eq "M"){
        my $len = shift(@stateLen);
        $curPos += $len;
        $protLen = $stateTotals[0];
    }
    if($protIV->end > $protLen){
        $Logger->Log("Attempt to convert coords for an interval which extends beyond the alignment","WARNING");
        return undef;
    }
    #On each loop update the interval by the deletion, then the curPos by the match
    while($curPos <= $protIV->end){
        unless(@stateLen){
            $Logger->Log("Attempt to convert coords of an lcr interval with a mismatched cigar","WARNING");
            return undef;
        }
        my $len = shift(@stateLen);
        $conIV->start($conIV->start + $len) if($curPos <= $protIV->start);
        $conIV->end($conIV->end + $len);
        $curPos += shift(@stateLen) if(@stateLen); #Aln might end with a del
    }
    return $conIV;
}

#Given a set of protein coordinate LCRs, and the cigars describing the msa
#   returns a set of intervals in consensus coordinates
#Inputs - An sID keyed hash of PROTINFO objects
#       - An sID keyed hash of CIGAR strings
#Output - A sorted array of non-overlapping INTERVAL objects
sub GetConsensusRegions(\%\%){
    my ($rRegions,$rCigars) = @_;
    $Logger->Log("Determining Consensus LCR Regions...","INFO");
    my @uSIDs = Unique(keys %{$rRegions}, keys %{$rCigars});
    my @missingIDs = grep {! exists $rRegions->{$_}} (@uSIDs);
    push(@missingIDs,grep {! exists $rCigars->{$_}} (@uSIDs));
    @missingIDs = Unique(@missingIDs);
    if(@missingIDs){
        $Logger->Log("The samples (@missingIDs) do not have both LCRs and CIGARS","ERROR");
    }
    my @cRegList;
    while(my ($sID,$infoObj) = each %{$rRegions}){
        foreach my $reg (@{$infoObj->regions}){
            my $cReg = ConvertCoords($reg,$rCigars->{$sID});
            unless(defined $cReg){
                $Logger->Log("Could not convert lcr coords for $sID(@{[$reg->start,$reg->end]}):Skipping","WARNING");
                next;
            }
            push(@cRegList,$cReg);
        }
    }
    #Collapse to unique, non-overlapping regions
    @cRegList = sort {$a->start <=> $b->start} @cRegList;
    for(my $i = 1; $i < @cRegList; $i++){
        if($cRegList[$i]->start <= $cRegList[$i-1]->end){
            if($cRegList[$i]->end > $cRegList[$i-1]->end){
                $cRegList[$i-1]->end($cRegList[$i]->end);
            }
            splice(@cRegList,$i--,1);
        }
    }
    $Logger->Log("Found @{[scalar(@cRegList)]} consensus LCR regions","INFO");
    return @cRegList;
}

#Given a set of intervals along the alignment, calculates the entropy in the consensus
#sequence for theose intervals
#Inputs - an array of INTERVAL objects
#       - a Bio::SimpleAln Object
#Output - an sID keyed hash of arrays of entropy values
#           Also has the special key of _Consensus_
#           Entropy values may be undefined, ie. all gaps at given range
sub CalcConsensusEntropy(\@$){
    my ($rIVList,$alnObj) = @_;
    $Logger->Log("Calculating entropy of consensus LCRs ...","INFO");
    my %seqDict = ("_Consensus_" => $alnObj->consensus_string());
    foreach my $seqObj ($alnObj->each_seq){
        my ($sID) = split(/:/,$seqObj->id);
        $seqDict{$sID} = $seqObj->seq;
    }
    my %hListDict;
    while(my ($sID,$seq) = each %seqDict){
        my @chars = map {(exists $_AAIDX{$_} or $_ eq $alnObj->gap_char) ? $_ : 'X'} split("",$seq);
        my @hList;
        foreach my $iv (@{$rIVList}){
            my @sv = (0) x @_AA;
            for(my $i = $iv->start - 1; $i < $iv->end; $i++){
                next if($chars[$i] eq $alnObj->gap_char);
                my @idx = @{$_AAIDX{$chars[$i]}};
                @sv[@idx] = map {$_ + 1/@idx} @sv[@idx];
            }
            push(@hList,CalcEntropy(@sv));
        }
        $hListDict{$sID} = \@hList;
    }
    $Logger->Log("Entropies calculated","DEBUG");
    return %hListDict;
}

#Given a set of intervals along an alignment, counts the number of non gap chars
#   in the interval for each protein
#Inputs - an array of INTERVAL objects
#       - a Bio::SimpleAln Object
#Output - an sID keyed hash of arrays of region lengths
sub CalcConsensusRegionLen(\@$){
    my ($rIVList,$alnObj) = @_;
    $Logger->Log("Determining sample specific consensus region lengths...","INFO");
    my %LenDict;
    foreach my $seqObj ($alnObj->each_seq){
        my ($sID) = split(/:/,$seqObj->id);
        my @chars = split(//,$seqObj->seq);
        my @lenList;
        foreach my $iv (@{$rIVList}){
            my @valid = grep {$_ ne $alnObj->gap_char} @chars[$iv->start - 1 .. $iv->end - 1];
            push(@lenList,scalar(@valid));
        }
        $LenDict{$sID} = \@lenList;
    }
    $Logger->Log("Lengths determined","DEBUG");
    return %LenDict;
}
