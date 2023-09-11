#!/usr/bin/perl
package RANGE;

sub new {
    my $package = shift;
    my $left = shift;
    my $right = shift || $left;
    my $self = {l => $left, r => $right};

    bless $self, $package;
    return $self;
}

sub copy {
    my $self = shift;
    return RANGE->new($self->{l},$self->{r});
}

sub left {
    my $self = shift;
    if(@_){
        $self->{l} = shift;
        $self->{l} = $self->{r} if($self->{l} > $self->{r});
    }
    return $self->{l};
}

sub right {
    my $self = shift;
    if(@_){
        $self->{r} = shift;
        $self->{r} = $self->{l} if($self->{r} < $self->{l});
    }
    return $self->{r};
}

sub length {
    my $self = shift;
    return $self->{r} - $self->{l} + 1;
}

sub intersect {
    my ($smol, $lorge) = sort {$a->{l} <=> $b->{l} or $a->{r} <=> $b->{r}} @_;
    my $left = $lorge->{l};
    my $right = $smol->{r};
    return undef if($left > $right);
    return RANGE->new($left,$right);
}

sub union {
    my ($self, $other) = @_;
    my ($left) = sort{$a <=> $b} ($self->{l},$other->{l});
    my ($right) = sort{$b <=> $a} ($self->{r},$other->{r});
    return RANGE->new($left,$right);
}

sub each_value {
    my $self = shift;
    return ($self->{l} .. $self->{r});
}

sub is_equal {
    my $self = shift;
    my $other = shift;
    return 0 if($self->{l} != $other->{l});
    return 0 if($self->{r} != $other->{r});
    return 1;
}

sub is_subrange {
    my $self = shift;
    my $other = shift;
    return 0 if($self->{l} < $other->{l});
    return 0 if($self->{r} > $other->{r});
    return 1;
}

sub contains {
    my $self = shift;
    my $val = shift;
    return ($val >= $self->{l} and $val <= $self->{r});
}

sub to_str {
    my $self = shift;
    return join("-",@{$self}{qw(l r)});
}

package main;
############################################################
# Description
# Given a tree and information on the lengths at the tips of the tree
# uses maximum parsimony to reconstruct the events that would explain the observed
# lengths
# Algorithm used is based on Fitch (1971) Society of Systematic Biologists
#   10.1093/sysbio/20.4.406

use warnings;
use strict;
use Bio::TreeIO;
use Class::Struct;
use File::Basename;
use POSIX qw(ceil);
use Zach::Logger;
use Zach::Util::Array qw(Unique);
use Zach::Util::File qw(OpenFileHandle);

############################################################
#Declarations

struct (CREG => {start => '$', end => '$', h => '$', lens => '%'});
struct (PROTINFO => {id => '$', regions => '@'});
struct (EVENT => {parent => '$', child => '$', pLen => '$', cLen => '$'});

sub LoadTree($);
sub ParseCReg($);
sub CalcMeanInterLeafDist($);
sub CalcMADLen($);
sub CalcMedianMADLen(@);
sub CalcIndelRate($@);
sub CheckForFixation($);
sub CheckForAutapomorphy($);
sub InitTree($$);
sub GetPrelimSet($);
sub FinalizeSet($);
sub CalcDeviationProb($$$);
sub DetermineEvents($$$);
sub LogNodeTags($;$);
sub RunParsimony($$$);
sub OutputEvents($$@);

#sub PrintTaggedTree($;$){
#sub PrintTaggedNode($$;$){


############################################################
#Handling Input

if(@ARGV < 2){
    die "Usage: ".basename($0). " Tree.newick regions.creg > Events.tsv\n";
}

#my $Logger = Zach::Logger->new(level => "WARNING");
#my $Logger = Zach::Logger->new(level => "INFO");
my $Logger = Zach::Logger->new(level => "DEBUG");

my %FileDict;
@FileDict{qw(tree creg)} = @ARGV[0 .. 2];
while (my ($type,$file) = each %FileDict){
    next if(-f $file);
    $Logger->Log("$type file ($file) does not exist","ERROR");
}

############################################################
# Main

sub main {
    my %fileDict = @_;
    my $treeObj = LoadTree($fileDict{tree});
    my @protInfoList = ParseCReg($fileDict{creg});
    my $indelRate = CalcIndelRate($treeObj,@protInfoList);
    print "ID\tRegion\tParent\tChild\tAncestral\tState\n";
    foreach my $protInfo (@protInfoList){
        $Logger->title($protInfo->id);
        my $counter = 0;
        foreach my $creg (@{$protInfo->regions}){
            my @events = RunParsimony($treeObj,$creg,$indelRate);
            OutputEvents($protInfo->id,$counter++,@events);
        }
    }
} main(%FileDict);

############################################################
# Subroutines

#Given a path to a newick formatted file, loads and returns the tree
#Inputs - A string represting a path to a newick formatted file
#Output - A Bio::Tree Object;
sub LoadTree($){
    $Logger->Log("Loading Tree...","DEBUG");
    my $in = Bio::TreeIO->new(-fh => OpenFileHandle(shift,"tree","ERROR"), -format => 'newick');
    my $tree = $in->next_tree;
    my $bValid = 1;
    my %Seen;
    my $nodeCount = 0;
    foreach my $node ($tree->get_nodes){
        unless(defined $node->id){
            $bValid = 0;
            last;
        }
        if($Seen{$node->id}){
            $bValid = 0;
            last;
        }
        $Seen{$node->id} = 1;
        $nodeCount++;
    }
    $Logger->Log("All nodes of the input tree must be uniquely labelled","ERROR") unless($bValid);
    $Logger->Log("Tree loaded with $nodeCount nodes","DEBUG");
    return $tree;
}

#Given a path to a creg file, loads information for each protein
#Inputs - A string representing a path to a creg formatted file (Output of get_consensus_lcrs.pl)
#Output - An array of PROTINFO objects
sub ParseCReg($){
    $Logger->Log("Parsing cregs...","DEBUG");
    my $fh = OpenFileHandle(shift,"creg","ERROR");
    my @infoList;
    my $curID = undef;
    my @curRegions = ();
    my $bInvalid = 0;
    while(my $line = <$fh>){
        chomp($line);
        next if($line =~ m/^$/);
        if($line =~ m/^#ID/){
            if(defined $curID){
                push(@infoList,PROTINFO->new(id => $curID, regions => [@curRegions]));
            }
            (undef, $curID) = split(/ /,$line);
            @curRegions = ();
            $bInvalid = 0;
            next;
        }
        next if($bInvalid);
        if($line =~ m/^#REG/){
            my (undef, $start,$end,$h) = split(/ /,$line);
            push(@curRegions,CREG->new(start => $start, end => $end, h => $h));
            next;
        }
        my ($sID,$minH,$cigar,$lenStr,$hStr) = split(/\t/,$line);
        my @lenList = split(/;/,$lenStr);
        if(@lenList != @curRegions){
            my $msg = sprintf("Length list for %s (%d) and region list for %s (%d) are non-conformable)",$sID,scalar @lenList,$curID,scalar @curRegions);
            $Logger->Log("$msg: Skipping $curID","WARNING");
            $bInvalid = 0;
            $curID = undef;
            next;
        }
        for(my $i = 0; $i < @lenList; $i++){
            $curRegions[$i]->lens->{$sID} = $lenList[$i];
        }
    }
    close($fh);
    if(defined $curID){
        push(@infoList,PROTINFO->new(id => $curID, regions => [@curRegions]));
    }
    $Logger->Log("Loaded information for @{[scalar(@infoList)]} proteins","INFO");
    return @infoList;
}

#Given a tree, calculate the average distance between leaves
sub CalcMeanInterLeafDist($){
    $Logger->Log("Calculating mean inter-leaf distance ...","DEBUG");
    my $treeObj = shift;
    my @leaves = $treeObj->get_leaf_nodes;
    my $total = 0;
    my $count = 0;
    for(my $i = 0; $i < $#leaves; $i++){
        for(my $j = $i+1; $j < @leaves; $j++){
            $total += $treeObj->distance(-nodes => [@leaves[($i,$j)]]);
            $count++;
        }
    }
    my $MILD = $total / $count;
    $Logger->Log("Mean Inter-leaf distance: $MILD","INFO");
    return $MILD;
}

#Given a creg, calculate the mean aboslute deviation in region length
sub CalcMADLen($){
    $Logger->Log("Calculating mean absolute deviation in length ...","DEBUG");
    my $creg = shift;
    my @x = values %{$creg->lens};
    my $mean = 0;
    $mean += $_ foreach (@x);
    $mean /= scalar(@x);
    my $MAD = 0;
    $MAD += abs($_ - $mean) foreach(@x);
    $MAD /= scalar(@x);
    $Logger->Log("Mean absolute deviation: $MAD","DEBUG");
    return $MAD;
}

#struct (CREG => {start => '$', end => '$', h => '$', lens => '%'});
#struct (PROTINFO => {id => '$', regions => '@'});
sub CalcMedianMADLen(@){
    $Logger->Log("Calculating Median MAD ...","DEBUG");
    my @MADList;
    foreach my $protInfo (@_){
        push(@MADList, map {CalcMADLen($_)} @{$protInfo->regions});
    }
    @MADList = sort {$a <=> $b} @MADList;
    my $n = scalar(@MADList)+1;
    my $median = ($MADList[int($n/2)-1] + $MADList[ceil($n/2)-1]) / 2;
    $median = 1/$n unless($median);
    $Logger->Log("Median MAD: $median","INFO");
    return $median;
}

#Given a tree and the region info, estimate the mutation rate in the protein length
#Assume an exponential distubution where λ = μt
#   set t to be the average distance between leaves of the tree
#   set the mean to be 1/2 of the Median MAD Len
#   λ is 1/mean across all regions in all proteins
#   solve for μ
#Inputs - A Bio::Tree::Tree Object
#       - an array PROTINFO Objects
#Output - an indel rate
sub CalcIndelRate($@){
    $Logger->Log("Calculating Indel Rate (2/t/mMAD)...","DEBUG");
    my ($treeObj, @protInfoList) = @_;
    my $rate = CalcMedianMADLen(@protInfoList) / CalcMeanInterLeafDist($treeObj);
    $Logger->Log("IndelRate: $rate","INFO");
    return $rate; 
}

#Given a creg object checks if all samples have the same length
#Inputs - A CREG object
#Output - If all samples have the same length, that length
#           else undefined
sub CheckForFixation($){
    my $creg = shift;
    my @lenList = Unique(values %{$creg->lens});
    if(@lenList == 1){
        return shift @lenList;
    }
    return undef;
}

#Given a creg object checks if all samples have the same length, except for one
#Inputs - A CREG Object
#Output - If there is an autapomorphy, the length of that sample
#           else undefined
sub CheckForAutapomorphy($){
    my $creg = shift;
    my %sIDByLen;
    while (my ($sID, $len) = each %{$creg->lens}) {
        $sIDByLen{$len} = [] unless(exists $sIDByLen{$len});
        push(@{$sIDByLen{$len}},$sID);
    }
    if(scalar(keys %sIDByLen) == 2){
        my @SIDListList = values %sIDByLen;
        foreach my $sIDList (@SIDListList){
            return shift(@{$sIDList}) if(scalar(@{$sIDList} == 1));
        }
    }
    return undef;
}

#Given a tree object, initializes the character states to be empty for internal nodes
#   and the known values at the tips
#Actually any node with a value in the provided input is labeled as such, so internal nodes
#   can be non-empty
#Inputs - A Bio::Tree::Tree Object
#       - a reference to an sID keyed hash of lengths
#Output - None, directly modifies the given tree
sub InitTree($$){
    $Logger->Log("Initializing tree ...","DEBUG");
    my ($tree, $rLenDict) = @_;
    foreach my $node ($tree->get_nodes){
        $node->remove_all_tags();
        $node->set_tag_value('min_change',0);
        if(exists $rLenDict->{$node->id}){
            $node->set_tag_value('range',RANGE->new($rLenDict->{$node->id}));
            $node->set_tag_value('defined',1);
        } else {
            $node->set_tag_value('range',undef);
        }
    }
    $Logger->Log("Tree Initialized","DEBUG");
}

#Recursivley determines the preliminary set of states for a node
#   The set for a node is the either the intersection of its children's sets if there is an
#   intersection, or the union of its children's sets
#   If a child does not have a defined set, then its set is determined first (recursion)
#Also counts the minimum number of changes required to account for the extant diversity
#Inputs - A Bio::Tree::Node object
#Outputs - a set of states for the child; Also modifies the tags of the node
sub GetPrelimSet($){
    my $node = shift;
    $Logger->Log("Determining Preliminary range for @{[$node->id]}...","DEBUG");
    my $range = $node->get_tag_values('range');
    return $range if(defined $range); #Node set is already defined
    my $minChange = 0;
    my @cRangeList;
    foreach my $child ($node->each_Descendent){
        my $cRange = GetPrelimSet($child);
        $minChange += ($child->get_tag_values('min_change'))[0];
        push(@cRangeList,$cRange) if(defined $cRange);
    }
    return undef unless(@cRangeList);
    my $union = shift @cRangeList;
    my $intersect = $union->copy;
    print STDERR $node->id, " \t";
    foreach my $cRange (@cRangeList){
        print STDERR $union->to_str(), " + ", $cRange->to_str(), " = ";
        $union = $union->union($cRange);
        print STDERR $union->to_str();
        next unless(defined $intersect);
        $intersect = $intersect->intersect($cRange);
    }
    print STDERR "\n";
    $range = $intersect;
    if(!defined $intersect){
        $range = $union;
        $node->set_tag_value('from_union',1);
        $minChange++;
    }
    $node->set_tag_value('min_change',$minChange);
    $node->set_tag_value('range',$range);
    $Logger->Log("Preliminary ranges determined for @{[$node->id]}","DEBUG");
    return $range;
}


#Given a non-root node, which has a preliminary set of states finalizes the set of
#   potential node states
#I - if the set is a superset of the parent -> II, else III
#II - Reduce to the intersection with the parent -> Done
#III - if formed from a union -> IV else -> V
#IV - Set to union with parent -> Done
#V - Add any states from the parent that are in the union of the children -> Done
#Tracking of prelim vs final range is done by the number of range tags, the first is the prelim, the
#second is the final 
#Inputs - A Bio::Tree::Tree Object
#Outputs - None, Modifies the node's tags
sub FinalizeSet($){
    my $node = shift;
    $Logger->Log("Determining final ranges for @{[$node->id]}...","DEBUG");
    my ($pRange) = $node->ancestor->get_tag_values('range');
    if($node->has_tag('defined')){
        my ($range) = $node->get_tag_values('range');
        $node->add_tag_value('range',$range->copy);
    } else {
        unless(defined $pRange){
            $Logger->Log("Too few leaves are defined to estimate ancestral states","ERROR");
        }
        my ($range) = $node->get_tag_values('range');
        if(defined $range and $pRange->is_subrange($range)){ # rule II
            $node->add_tag_value('range',$range->intersect($pRange));
        } elsif(!defined $range){ #Special case of rule IV if missing info on all descendents
            $node->add_tag_value('range',$pRange->copy);
        } elsif($node->has_tag('from_union')){ # rule IV
            $node->add_tag_value('range',$range->union($pRange));
        } else { #Rule 5
            my @cRangeList = map {($_->get_tag_values('range'))} $node->each_Descendent;
            @cRangeList = grep {defined $_} @cRangeList;
            my $cUnion = shift @cRangeList;
            foreach my $cRange (@cRangeList){
                $cUnion = $cUnion->union($cRange);
            }
            my $potentialAdds = $pRange->intersect($cUnion);
            $node->add_tag_value('range',$range->union($potentialAdds));
        }
    }
    foreach my $child ($node->each_Descendent){
        FinalizeSet($child);
    }
    $Logger->Log("Final ranges determined for @{[$node->id]}","DEBUG");
}

#Given parent and child lengths and a branch length, calculates the probability of observing that change
#in length assuming a poisson distribution
#Inputs - A difference in region length between parent and child
#       - A branch length
#       - The indel Rate
#Output - the probability of observing cLen given pLen after bLen time
sub CalcDeviationProb($$$){
    my ($deviation,$bLen,$indelRate) = @_;
    my $lambda = $bLen * $indelRate;
    my $p = $lambda**$deviation * exp(-$lambda);
    for(my $i = 1; $i < $deviation; $i++){
        $p /= $i;
    }
    return $p
}

#Given a proposed state for a node, and the node itself, calculates the best series of events creating
#that event, and the "probability" of the event chain
#Inputs - A proposed state for the node
#       - a Bio::Tree::Node object
#Output - A probability
#       - An array of EVENT objects
#           Note: The output is a single array with the probability as the first element
sub DetermineEvents($$$){
    my $state = shift;
    my $node = shift;
    my $indelRate = shift;
    $Logger->Log("Determining the best set of events if @{[$node->id]} had a length of $state ...","DEBUG");
    if($node->is_Leaf){
        return 1;
    }
    my ($prelimRange,$finalRange) = $node->get_tag_values('range');
    my $bPrelim = $prelimRange->contains($state);
    my $probability = 1;
    my @events;
    foreach my $child ($node->each_Descendent){
        my ($cPrelimRange,$cFinalRange) = $child->get_tag_values('range');
        my @cStateList;
        if($cPrelimRange->contains($state)){
            @cStateList = ($state);
        } else {
            @cStateList = $cPrelimRange->each_value;
            push(@cStateList,$state) if($cFinalRange->contains($state));
        }
        my %best = (prob => -1, events => undef);
        foreach my $cState (@cStateList){
            my ($prob, @cEvents) = DetermineEvents($cState,$child,$indelRate);
            $prob *= CalcDeviationProb(abs($state-$cState),$child->branch_length,$indelRate);
            unless($state eq $cState){
                 unshift(@cEvents,EVENT->new(parent => $node->id, child => $child->id, pLen => $state, cLen => $cState));
            }

            if($prob > $best{prob}){
                $best{prob} = $prob;
                $best{events} = \@cEvents;
            }
        }
        $probability *= $best{prob};
        unless($probability > -1){
            $Logger->Log("Failure to process children of non-leaf node (@{[$node->id]})","ERROR");
        }

        push(@events,@{$best{events}}) if(defined $best{events});
    }
    $Logger->Log("(@{[$node->id]}:$state) Found a set of @{[scalar @events]} events with -log(p)=@{[sprintf('%0.2f',-log($probability))]}","DEBUG") if($probability > 0);
    return ($probability, @events);
}

#Prints to STDERR the tags of a sub-tree in depth first order
sub LogNodeTags($;$){
    return unless($Logger->level eq "DEBUG");
    my $node = shift;
    my $level = shift || 0;
    my $str = ($level) ? " " x ($level - 1) . "∟" : "";
    $str .= $node->id;
    #    min_change, range, defined, from_union
    my @parens = ($node->has_tag('defined')) ? qw([ ]) : qw(( ));
    my $arrow  = $node->has_tag('from_union') ? '=>' : '->';
    my ($prelim,$final) = $node->get_tag_values('range');
    my $pStr = (defined $prelim) ? $prelim->to_str : "_";
    my $fStr = (defined $final) ? $final->to_str : "_";
    $str .= join(" ",($parens[0], $pStr, $arrow, $fStr,$parens[1]));
    print STDERR $str,"\n";
    foreach my $child ($node->each_Descendent){
        LogNodeTags($child,$level+1);
    }
}

#Given a tree and the lengths of extent regions, output the events which led to this state
#Inputs - A Bio::Tree::Tree Object
#       - A CREG Object
#       - an indel rate
#Output - An array of EVENT objects
sub RunParsimony($$$){
    my ($tree, $creg, $indelRate) = @_;
    $Logger->Log("Running Parsimony ...","DEBUG");
    if(my $fixedLen = CheckForFixation($creg)){
        return (EVENT->new(child => $tree->get_root_node->id, cLen => $fixedLen));
    }
    if(my $sID = CheckForAutapomorphy($creg)){
        my ($node) = $tree->find_node(-id => $sID);
        my $pLen;
        foreach my $otherSID (keys %{$creg->lens}){
            next if($otherSID eq $sID);
            $pLen = $creg->lens->{$otherSID};
            last;
        }
        my @events = EVENT->new(child => $tree->get_root_node->id, cLen => $pLen);
        push(@events,EVENT->new(parent => $node->ancestor->id, child => $sID, pLen => $pLen, cLen => $creg->lens->{$sID}));
        return @events;
    }
    InitTree($tree,$creg->lens);
    LogNodeTags($tree->get_root_node);
    GetPrelimSet($tree->get_root_node);
    LogNodeTags($tree->get_root_node);
    $tree->get_root_node->add_tag_value('range',($tree->get_root_node->get_tag_values('range'))[0]);
    foreach my $child ($tree->get_root_node->each_Descendent){
        FinalizeSet($child);
    }
    LogNodeTags($tree->get_root_node);
    my ($prelimRange,$finalRange) = $tree->get_root_node->get_tag_values('range');
    my %best = (prob => -1,state => undef, events => undef);
    foreach my $state ($finalRange->each_value){
        my ($prob,@events) = DetermineEvents($state,$tree->get_root_node,$indelRate);
        if($prob > $best{prob}){
            $best{prob} = $prob;
            $best{state} = $state;
            $best{events} = \@events;
        }
    }
    my @events = @{$best{events}};
    unshift(@events,EVENT->new(child =>$tree->get_root_node->id,cLen => $best{state}));
    $Logger->Log("Final set of @{[scalar @events]} events with -log(p)=@{[sprintf('%0.2f',-log($best{prob}))]}","INFO") if($best{prob} > 0);
    return @events;
}

#Given the parsimony events for a region in a protein, outputs them to stdout
#Inputs - a string represting an id for the protein
#       - a number indicating the index of the region within the protein
#       - an array of EVENT objects
#Outputs - None, writes to STDOUT
sub OutputEvents($$@){
    my ($tID,$idx,@events) = @_;
    foreach my $event (@events){
        unless(defined $event->parent){
            $event->parent("ANCESTRAL");
            $event->pLen($event->cLen);
        }
        print join("\t",($tID,$idx,$event->parent,$event->child,$event->pLen,$event->cLen)),"\n";
    }
}

##Given a tree object, prints it out in a way that show tag information at nodes
##Inputs - a Bio::Tree obj
##       - an optional filehandle; defaults to stderr
##Output - None, writes to filehandle
#sub PrintTaggedTree($;$){
#    my $tree = shift;
#    my $fh = shift || *STDERR;
#    my @depthList = (0);
#    my @nodeList = ($tree->get_root_node());
#    while(@nodeList){
#        my $depth = shift @depthList;
#        my $node = shift @nodeList;
#        my @children = PrintTaggedTree($node,$depth,$fh);
#        push(@depthList,($depth+1) x scalar(@children));
#        push(@nodeList,@children);
#    }
#}
#
##Given a node object, prints it out with all information
##Inputs - A Bio::Node object
##       - a current node depth
##       - an optional filehandle; defaults to stderr
##Output - An array of node objects, contains the children of this node
#sub PrintTaggedNode($$;$){
#    my $node = shift;
#    my $depth = shift;
#    my $fh = shift || *STDERR;
#    my $str = "-" x $depth;
#    return $node->each_Descendent;
#}


