#!/usr/bin/perl -w

package MakeMat;

use strict;
use base 'Tool';



sub new
{
    my ($class, $exe_path, $prefix, $chk_name, $fasta_name) = @_;
    
    my $self = {
                 'exe'        => "$exe_path/makemat",
                 'prefix'     => $prefix,
                 'chk_name'   => $chk_name,
                 'fasta_name' => $fasta_name,
                 'err'        => 0
               };
    
    bless $self, $class;
    
    $self->{'cmd'} = "$self->{'exe'} -P $self->{'prefix'}";
    
    return $self;
}

sub chk_name
{
    my ($self, $chk_name) = @_;
    
    $self->{'chk_name'} = $chk_name if (defined($chk_name));
    return $self->{'chk_name'};
}

sub fasta_name
{
    my ($self, $fasta_name) = @_;
    
    $self->{'fasta_name'} = $fasta_name if (defined($fasta_name));
    return $self->{'fasta_name'};
}

sub init
{
    my ($self) = @_;
    
    open(PN, ">", $self->prefix() . ".pn");
    print PN $self->chk_name() . "\n";
    close PN;
    
    open(SN, ">", $self->prefix() . ".sn");
    print SN $self->fasta_name() . "\n";
    close(SN);
}

sub prefix
{
    my ($self, $prefix) = @_;
    
    $self->{'prefix'} = $prefix if (defined($prefix));
    return $self->{'prefix'};
}

sub run
{
    my ($self, $error) = @_;
    $error = 0 unless ($error);
    
    $self->err(1) and return unless ( (-f $self->prefix() . '.pn') && 
                                      (-f $self->prefix() . '.sn') && 
                                      (-s $self->prefix() . '.chk') );
    
    print STDERR "Running makemat: ";
    $self->SUPER::run($error);
}



1;
