#!/usr/bin/perl -w
use strict;
use lib '/webdata/binaries/current/FFPred3';
use DrawUtils::Transmembrane;



$| = 1;

main:
{
    my $id = $ARGV[0];
    my $OB = {};
    
    getDataFromFile($OB, $id);
    createImage($OB, $id);
    printImage($OB, $id);
}

sub getDataFromFile
{
    my ($OBJ, $id) = @_;
    
    $OBJ->{'SIGNAL'}   = 0;
    $OBJ->{'ANCHOR'}   = 0;
    $OBJ->{'NTERM'}    = 'in';
    $OBJ->{'TM_TOPOL'} = "";
    
    my $i = 0;
    
    open(IN, "<", "$id.featcfg");
    while (<IN>)
    {
        chomp $_;
        my @tmp = split(/\t/,$_);
        
        $OBJ->{'SIGNAL'} = 1 if ($tmp[1] =~ /SIGNAL/);
        $OBJ->{'ANCHOR'} = 1 if ($tmp[1] =~ /ANCHOR/);    
        $OBJ->{'NTERM'}  = 'out' if (($tmp[1] =~ /SIGNAL/) || ($tmp[1] =~ /ANCHOR/));
        
        if ($tmp[0] =~ /TM/)
        {
            $OBJ->{'TM'} = 1;      
            $OBJ->{'TM_TOPOL'} .= "$i\.$tmp[2]\,$tmp[3]\;";
            $i++;
        }   
    }
    close(IN);
    
    chop $OBJ->{'TM_TOPOL'} if ($OBJ->{'TM_TOPOL'} =~ /\;$/);
}

sub createImage
{
    my ($OBJ, $id) = @_;
    
    if (defined($OBJ->{'TM'}))
    {
	$OBJ->{'IMG'} = new Transmembrane(
                                           -outside_label      => 'Extracellular',
                                           -inside_label       => 'Cytoplasm',
                                           -membrane_label     => 'Membrane',
                                           -vertical_padding   => 50,
                                           -horizontal_padding => 100,
                                           -n_terminal_offset  => 50,
                                           -n_terminal_height  => 220,
                                           -c_terminal_offset  => 30,
                                           -c_terminal_height  => 220,
		     		           -n_terminal         => $OBJ->{'NTERM'},
                                           -signal             => $OBJ->{'SIGNAL'},
					   -anchor             => $OBJ->{'ANCHOR'},
                                           -helix_height       => 60,
                                           -helix_width        => 30,
                                           -short_loop_limit   => 5,
                                           -long_loop_limit    => 10,
                                           -n_terminal_height  => 50,
                                           -c_terminal_height  => 50,
                                           -topology_string    => $OBJ->{'TM_TOPOL'},
                                           #-labels             => $OBJ->{'LABELS'},
                                           -helix_label        => 'helix',
                                           -color_scheme       => 'yellow'
                                         );
    }
}

sub printImage
{
    my ($OBJ, $id)= @_;
    my $fhPNG = new FileHandle("${id}_tm.png", 'w');
    
    # Next block edited in this web server version in order to be compatible
    # with FFPred 3, see "/webdata/binaries/current/FFPred3/README".
    if (defined($OBJ->{'IMG'}))
    {
        print $fhPNG $OBJ->{'IMG'}->png;
    }
    else
    {
        print $fhPNG "NO TRANSMEMBRANE SEGMENTS\n";
    }
    
    $fhPNG->close;
}
