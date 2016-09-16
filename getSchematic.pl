#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use DrawUtils::Schematic;

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

    open(IN, "<", "$id.featcfg") or die "can't open featfcfg";

    my ($seqid, $seq) = split(/\t/,<IN>);
    $OBJ->{'LEN'} = length $seq;

    while (<IN>)
    {
        chomp $_;
        my @tmp = split(/\t/,$_);

        # Next line edited in this web server version in order to be compatible
        # with FFPred 3, see "/webdata/binaries/current/FFPred3/README".
        if ((@tmp == 5) && ($tmp[0] !~ /SP/))
        {
            my ($name, $type, $from, $to, $score) = @tmp;
            my $i = defined($OBJ->{'FEAT_ARRAY'}{$name}) ? scalar(@{$OBJ->{'FEAT_ARRAY'}{$name}}) : 0;

            $OBJ->{'FEAT_ARRAY'}{$name}[$i]{'from'} = $from;
            $OBJ->{'FEAT_ARRAY'}{$name}[$i]{'to'}   = $to;
            $OBJ->{'FEAT_ARRAY'}{$name}[$i]{'val'}  = $score;
            $OBJ->{'FEAT_ARRAY'}{$name}[$i]{'type'} = $type;
            $OBJ->{'FEAT_ARRAY'}{$name}[$i]{'name'} = $name;
        }
    }

    close(IN);
}

sub createImage
{
    my ($OBJ, $id) = @_;

    $OBJ->{'IMG'} = new Schematic(
                                   -vertical_padding   => 50,
                                   -horizontal_padding => 8,
                                   -offset             => 0,
                                   -max_height         => 500,
                                   -max_width          => 720,
                                   -feat_hash          => $OBJ->{'FEAT_ARRAY'},
                                   -len                => $OBJ->{'LEN'},
                                   -color_scheme       => 'default',
                                   -fonts              => $FindBin::Bin.'/DrawUtils/fonts/',
                                   -ttf_font           => 'Verdana',
                                   -ttf_font_size      =>  8
                                 );
}

sub printImage
{
    my ($OBJ, $id)= @_;
    my $fhPNG = new FileHandle("${id}_sch.png", 'w');

    print $fhPNG $OBJ->{'IMG'}->png;

    $fhPNG->close;
}
