#!/usr/bin/perl -w

package Pred;

use strict;
use base 'Tool';



sub new
{
    my ($class, $aa, $id, $md5, $path, $name) = @_;
    
    my $self = {
                 'aa'       => $aa,
                 'id'       => $id,
                 'md5'      => $md5,
                 'path'     => $path,
                 'name'     => $name, 
                 'results'  => {},
                 'err'      => 0
               };
    
    bless $self, $class;
    
    $self->{'len'}  = length($self->{'aa'});
    $self->{'seg3'} = $self->{'len'} > 102 ? ($self->{'len'} - 100)/3 + 1 : 1;
    $self->{'seg8'} = $self->{'len'} > 107 ? ($self->{'len'} - 100)/8 + 1 : 1;
    
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

    $RES->{"${name}_cterm"}++;
}

sub addMidSegment
{
    my ($self, $RES, $idx, $name) = @_;

    return if ($self->len() < 108);
 
    my $seg = 1 + int(($idx-50) / $self->seg8());
    $RES->{"${name}_S$seg"}++;
}

sub addNterm
{
    my ($self, $RES, $name) = @_;

    $RES->{"${name}_nterm"}++;
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

sub name
{
    my ($self, $name) = @_;
    
    $self->{'name'} = $name if (defined($name));
    return $self->{'name'};
}

sub normalise
{
    my ($self) = @_;
    
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
    my ($self, $outfile) = @_;
    my $name = $self->name();
    my $config = $self->{$name}; # This weird organisation comes from legacy code.
    
    open(OUT, ">>", $outfile);
    
    if (ref($config) eq 'HASH')
    {
        foreach my $key (keys %$config)
        {
            my $val = exists($config->{$key}{'val'}) ? $config->{$key}{'val'} : 0;
            print OUT "$name\t$key\t$val\n";
        }
        
        if ($name eq 'SF')
        {
            foreach my $key (keys %{$self->{'AA'}})
            {
                my $val = exists($self->{'AA'}{$key}{'val'}) ? $self->{'AA'}{$key}{'val'} : 0;
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
    my ($self, $outfile, $program) = @_;
    my $name = $self->name();
    my $results = $self->results();
    
    open(OUTFILE, ">>", $outfile);
    print OUTFILE "ID  \\  $program\t";
    
    foreach my $key (sort keys %$results)
    {
        if (ref($results->{$key}) eq 'HASH')
        {
            foreach my $sub (sort {$a <=> $b} keys %{$results->{$key}})
            {
                print OUTFILE "$name.$key.$sub\t";
            }
        }
        else
        {
            print OUTFILE "$name.$key\t";
        }
    }
    
    close(OUTFILE); 
}

sub print_results
{
    my ($self, $outfile, $type) = @_;
    my $results = $self->results();
    
    open(OUTFILE, ">>", $outfile);
    print OUTFILE "${type}_" . $self->md5() . "\t";
    
    foreach my $key (sort keys %$results)
    {
        if (ref($results->{$key}) eq 'HASH')
        {
            foreach my $sub (sort {$a <=> $b} keys %{$results->{$key}})
            {
                print OUTFILE "$results->{$key}{$sub}\t";
            }
        }
        else
        {
            print OUTFILE "$results->{$key}\t";
        }
    }
    
    close(OUTFILE);
}

sub results
{
    my ($self, $results) = @_;
    
    $self->{'results'} = $results if (defined($results));
    return $self->{'results'};
}

sub seg3
{
    my ($self, $seg3) = @_;
    
    $self->{'seg3'} = $seg3 if (defined($seg3));
    return $self->{'seg3'};
}

sub seg8
{
    my ($self, $seg8) = @_;
    
    $self->{'seg8'} = $seg8 if (defined($seg8));
    return $self->{'seg8'};
}



1;
