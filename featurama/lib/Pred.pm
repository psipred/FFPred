#!/usr/bin/perl -w

package Pred;

use strict;
use base 'Tool';
use Data::Dumper;


sub new
{
    my ($class, $aa, $id, $md5, $path, $name) = @_;

    my $self = {
                 'aa'         => $aa,
                 'id'         => $id,
                 'md5'        => $md5,
                 'path'       => $path,
                 'name'       => $name,
                 'results'    => {},
                 'terminus'   => 20,
                 'termBuffer' => 5,
                 'minSegment' => 2,
                 'err'        => 0
               };

    $self->{'len'} = length($self->{'aa'});

    # bless $self, $class;
    return $self;
}

sub aa
{
    my ($self, $aa) = @_;

    $self->{'aa'} = $aa if (defined($aa));
    return $self->{'aa'};
}

sub addCterm
{
    my ($self, $RES, $name) = @_;
    my $Cterm_length = $self->terminus(); # Number of a.a.s in C-terminus.

    $RES->{"${name}_cterm"} += 1/$Cterm_length;
}

sub addMidSegment
{
    my ($self, $RES, $segments, $internal, $name) = @_;
    my $md5 = $self->md5();

    my $good_segment;
    foreach my $segment (sort {$a <=> $b} keys %$segments)
    {
        my $segment_length = $segments->{$segment};
        die "Bad segment composition for $md5, $name, internal $internal !\n" unless(defined($segment_length));
        $internal -= $segment_length;
        $good_segment = $segment and last unless ($internal > 0);
    }

    die "Bad (too high) internal position passed to sub addMidSegment for $md5, $name, internal $internal !\n" unless (defined($good_segment));

    $RES->{"${name}_S${good_segment}"} += 1/$segments->{$good_segment};
}

sub addNterm
{
    my ($self, $RES, $name) = @_;
    my $Nterm_length = $self->terminus(); # Number of a.a.s in N-terminus.

    $RES->{"${name}_nterm"} += 1/$Nterm_length;
}

sub addResidue
{
    my ($self, $RES, $idx, $name) = @_;
    # Note that any residue is part of either one or two segments, including
    # N-terminus, internal segments (which may not exist for short sequences),
    # and C-terminus - no residue is excluded and no residue is allowed to be
    # part of more than 2 segments.

    my $len = $self->len();
    my $terminus = $self->terminus();

    $self->addNterm($RES, $name) if ($idx <= $terminus);
    $self->addCterm($RES, $name) if ($idx > ($len - $terminus));

    # Do nothing more if the sequence is too short to have internal segments.
    my $segments = $self->{'segments'};
    return unless (defined($segments));

    # Do nothing more if the a.a. is too extremal to be included in one of the
    # internal segments. Note that for short sequences, part of the N- and
    # C-terminus may overlap with internal segments.
    my $min_inside = $self->numSegments() * $self->minSegment(); # Minimum number of a.a.s in internal segments.
    my $overlap = ($len < (2 * $terminus + $min_inside)); # Condition that is true only for short sequences where segments overlap with N- or C-terminus.
    my $Nterm = $overlap ? int(($len-$min_inside)/2) : $terminus; # Number of a.a.s at N-terminus not belonging to an internal segment.
    my $Cterm = $overlap ? ($len-$min_inside-$Nterm) : $terminus; # Number of a.a.s at C-terminus not belonging to an internal segment.

    return if (($idx <= $Nterm) || ($idx > ($len - $Cterm)));

    # If all is well, add the a.a. to the appropriate internal segment.
    my $internal = $idx - $Nterm;
    $self->addMidSegment($RES, $segments, $internal, $name);
}

sub addThresholdedRegion
{
    my ($self, $RES, $CFG, $name, $thresholds, $start, $length, $score) = @_;
    $score = 0 unless (defined($score));

    my $good_region = "";

    if (@$thresholds > 1)
    {
        for (my $i = 0; $i < @$thresholds; $i++)
        {
            if (($length >= $thresholds->[$i]) && (($i == @$thresholds-1) || ($length < $thresholds->[$i+1])))
            {
                $RES->{"num_${name}_$thresholds->[$i]+"}++;
                $good_region = 1;
            }
        }
    }
    else
    {
        if ($length >= $thresholds->[0])
        {
            $RES->{"num_${name}_regions"}++;
            $good_region = 1;
        }
    }

    if ($good_region)
    {
        my $n = scalar @$CFG;
        $CFG->[$n]{'type'} = $name;
        $CFG->[$n]{'from'} = $start;
        $CFG->[$n]{'to'} = ($start + $length - 1);
        $CFG->[$n]{'score'} = $score;
    }
}

sub createSegments
{
    my ($self, $num_segments) = @_;
    my $md5 = $self->md5();
    my $len = $self->len();
    my $terminus = $self->terminus();
    my $min_segment_length = $self->minSegment();

    my $segments = {}; # This hash will contain the lengths of all internal segments.

    # Do nothing if sequence is too short to have internal segments.
    my $min_termini = 2 * $terminus; # Minimum length of a sequence whose N- and C-terminus are not overlapping.
    my $min_inside = $num_segments * $min_segment_length; # Minimum number of a.a.s in internal segments.
    my $terminus_buffer = $self->termBuffer(); # Minimum distance of internal segments from ends of sequence.
    return undef if (($len < $min_termini) || ($len < ($min_inside + 2 * $terminus_buffer)));

    # Calculate lengths of the various segments.
    my $overlap = ($len < ($min_termini+$min_inside)); # Condition that is true only for short sequences where segments overlap with N- or C-terminus.
    my $short_segment = $overlap ? $min_segment_length : int(($len-$min_termini)/$num_segments); # Length of the shortest segment(s).
    my $remainder = $overlap ? 0 : (($len-$min_termini) % $num_segments); # Number of extra a.a.s after giving $short_segment a.a.s to each segment.

    for (my $counter = 1; $counter <= $num_segments; $counter++)
    {
        $segments->{$counter} = $short_segment;
    }
    # Here, adjust length of segments when number of internal residues is not a
    # multiple of the number of segments. Start increasing central segments as
    # $len increases.
    my $increase = 0;
    while (++$increase <= $remainder)
    {
        # This formula allows to consistently increase central segments before
        # end segments, and central segments on the N-terminus side before the
        # corresponding segments on the C-terminus side.
        my $segments_parity = $num_segments % 2;
        my $segment = (($num_segments+$segments_parity)/2 + ((-1)**($increase+$segments_parity)) * int($increase/2));
        die "Error when creating segments for $md5 !\n" unless (exists($segments->{$segment}));
        $segments->{$segment}++;
    }

    return $segments;
}

sub len
{
    my ($self, $len) = @_;

    $self->{'len'} = $len if (defined($len));
    return $self->{'len'};
}

sub md5
{
    my ($self, $md5) = @_;

    $self->{'md5'} = $md5 if (defined($md5));
    return $self->{'md5'};
}

sub minSegment
{
    my ($self, $minSegment) = @_;

    $self->{'minSegment'} = $minSegment if (defined($minSegment));
    return $self->{'minSegment'};
}

sub name
{
    my ($self, $name) = @_;

    $self->{'name'} = $name if (defined($name));
    return $self->{'name'};
}

sub numSegments
{
    my ($self, $numSegments) = @_;

    $self->{'numSegments'} = $numSegments if (defined($numSegments));
    return $self->{'numSegments'};
}

sub normalise
{
    my ($self) = @_;

    # Default case: do nothing and report success.
    return 1;
}

sub path
{
    my ($self, $path) = @_;

    $self->{'path'} = $path if (defined($path));
    return $self->{'path'};
}

sub print_config
{
    my ($self, $output_file) = @_;
    my $name = $self->name();
    my $config = $self->{$name}; # This weird organisation comes from legacy code.

    open(OUT, ">>", $output_file);

    if (ref($config) eq 'HASH')
    {
        foreach my $key (keys %$config)
        {
            my $val = defined($config->{$key}) ? $config->{$key} : 0;
            print OUT "$name\t$key\t$val\n";
        }

        if ($name eq 'SF')
        {
            foreach my $key (keys %{$self->{'AA'}})
            {
                my $val = defined($self->{'AA'}{$key}) ? $self->{'AA'}{$key} : 0;
                print OUT "AA\t$key\t$val\t" . (100 * $val) . "\n";
            }
        }
    }
    elsif (ref($config) eq 'ARRAY')
    {
        foreach my $idx (@$config)
        {
            my $score = exists($idx->{'score'}) ? $idx->{'score'} : 0;
            my $type  = exists($idx->{'type'}) ? $idx->{'type'} : '---';
            print OUT "$name\t" .
                      "$type\t" .
                      "$idx->{'from'}\t" .
                      "$idx->{'to'}\t" .
                      "$score\n";
        }
    }

    close(OUT);
}

sub print_keys
{
    my ($self, $output_file, $program) = @_;
    my $name = $self->name();
    my $results = $self->results();

    open(OUTPUT, ">>", $output_file);
    print OUTPUT "ID  \\  $program\t";

    foreach my $key (sort keys %$results)
    {
        if (ref($results->{$key}) eq 'HASH')
        {
            foreach my $sub (sort {$a <=> $b} keys %{$results->{$key}})
            {
                print OUTPUT "$name.$key.$sub\t";
            }
        }
        else
        {
            print OUTPUT "$name.$key\t";
        }
    }

    close(OUTPUT);
}

sub print_results
{
    my ($self, $output_file, $type) = @_;
    my $results = $self->results();

    open(OUTPUT, ">>", $output_file);
    print OUTPUT "${type}_" . $self->md5() . "\t";

    foreach my $key (sort keys %$results)
    {
        if (ref($results->{$key}) eq 'HASH')
        {
            foreach my $sub (sort {$a <=> $b} keys %{$results->{$key}})
            {
                print OUTPUT "$results->{$key}{$sub}\t";
            }
        }
        else
        {
            print OUTPUT "$results->{$key}\t";
        }
    }

    close(OUTPUT);
}

sub results
{
    my ($self, $results) = @_;

    $self->{'results'} = $results if (defined($results));
    return $self->{'results'};
}

sub termBuffer
{
    my ($self, $termBuffer) = @_;

    $self->{'termBuffer'} = $termBuffer if (defined($termBuffer));
    return $self->{'termBuffer'};
}

sub terminus
{
    my ($self, $terminus) = @_;

    $self->{'terminus'} = $terminus if (defined($terminus));
    return $self->{'terminus'};
}


1;
