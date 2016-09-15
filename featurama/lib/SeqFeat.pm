#!/usr/bin/perl -w

package SeqFeat;

use strict;
use base 'Pred';



sub new
{
    my ($class, $aa, $id, $md5, $cfg) = @_;
    my $name = 'SF';
    
    my $self = $class->SUPER::new($aa, $id, $md5, $cfg->{'PATH'}, $name);
    
    $self->{'exe'} = "$cfg->{'SEQFEAT'}/features";
    $self->{'cmd'} = "$self->{'exe'}" . 
                     " -i $self->{'path'}/$self->{'md5'}.fsa >" . 
                     " $self->{'path'}/$self->{'md5'}.seqfeat "; 
    
    $self->{$self->name()} = {}; # This weird organisation comes from legacy code.
    $self->{'AA'}          = {};
    
    return $self;
}

sub normalise
{
    my ($self) = @_;

    #log frequency data types

    $self->{'results'}{'len'}        = (log($self->{'results'}{'len'})-log(2))/(log(34350)-log(2));
    $self->{'results'}{'mwt'}        = (log($self->{'results'}{'mwt'})-log(149.2000671387)) / (log(3816170)-log(149.2000671387));
    $self->{'results'}{'vol'}        = (log($self->{'results'}{'vol'})-log(162.8999938948)) / (log(4592520)-log(162.89999389648));
    $self->{'results'}{'surf'}       = (log($self->{'results'}{'surf'})-log(185)) / (log(6262791)-log(185));
    $self->{'results'}{'hydro'}      = ($self->{'results'}{'hydro'}+16072.599609375) / (7980.0200195312+16072.599609375); 
    $self->{'results'}{'ave_hydro'}  = log($self->{'results'}{'ave_hydro'}+4.372373) / log(65.41+4.372373);
    $self->{'results'}{'charge'}     = ($self->{'results'}{'charge'}+356.111) / (504.622+356.111);
    $self->{'results'}{'mol_ext'}    = (log(1+$self->{'results'}{'mol_ext'})-log(436)) / (log(3930261)-log(436)) if ($self->{'results'}{'mol_ext'} > 0);
    $self->{'results'}{'iso_pt'}     = ($self->{'results'}{'iso_pt'}-2.14645) / (14-2.14645);
    $self->{'results'}{'aliphatic'} /= 234;
    $self->{'results'}{'num_atoms'}  = log($self->{'results'}{'num_atoms'}-17) / log(539031-17);
}

sub parse
{
    my ($self) = @_;
    
    my $RES = {};
    my $CFG = $self->{$self->name()};
    
    my @AA  = ('A','R','N','D','C','Q','E','G','H',
               'I','L','K','M','F','P','S','T','W',
               'Y','V');
    my @ATOM= ('atomC','atomH','atomO','atomN','atomS');
    
    if (-s "$self->{'path'}/$self->{'md5'}.seqfeat")
    {
        open(SF, "< $self->{'path'}/$self->{'md5'}.seqfeat");
    
        while(<SF>)
        {
            chomp $_;
            next if $_ =~ /^\s*$/;
            
            my ($sid, $met, $len, $mwt, $vol, $surf, $hydro, $meanH, $charge, $mol_ext, $iso, $aliphat, $npos, $nneg, $ppos, $pneg) = split(/\s+/,$_);
            
            $RES->{'len'}       = $len;
            $RES->{'mwt'}       = $mwt;
            $RES->{'vol'}       = $vol;
            $RES->{'surf'}      = $surf;
            $RES->{'hydro'}     = $hydro;
            $RES->{'ave_hydro'} = $meanH;
            $RES->{'charge'}    = $charge;
            $RES->{'mol_ext'}   = $mol_ext;
            $RES->{'iso_pt'}    = $iso;
            $RES->{'aliphatic'} = $aliphat;
            $RES->{'npos'}      = $ppos;
            $RES->{'nneg'}      = $pneg;
            
            $CFG->{'molecular weight'}{'val'}             = $mwt;
            $CFG->{'hydrophobicity'}{'val'}               = $meanH;
            $CFG->{'molar extinction coefficient'}{'val'} = $mol_ext;
            $CFG->{'isoelectric point'}{'val'}            = $iso;
            $CFG->{'charge'}{'val'}                       = $charge;
            $CFG->{'aliphatic index'}{'val'}              = $aliphat;
            $CFG->{'percent positive residues'}{'val'}    = $ppos*100;
            $CFG->{'percent negative residues'}{'val'}    = $pneg*100;
            
            my @tmp = split(/\s+/,$_);
            my $x = 0;
            
            for(my $i = 17; $i < 56; $i += 2, $x++)
            {
                $RES->{$AA[$x]} = $tmp[$i];
                $self->{'AA'}{$AA[$x]}{'val'} = $tmp[$i];
            }
            
            $x = 0;
            for(my $i = 63; $i < 72; $i += 2, $x++)
            {
                $RES->{$ATOM[$x]} = $tmp[$i];
            }
            
            $RES->{'num_atoms'} = pop @tmp;
        }
        
        close(SF);
    }
    
    $self->{'results'} = $RES;
}



1;
