package Transmembrane;

use strict;
use warnings;
use GD;
our @ISA = "GD"; 
our @EXPORT = qw(new add_tmhmm_feat png);
our $VERSION = '1.13';

sub new {

	my $class = shift;
 	my %options = @_;
	
	my $self = {

		## general paramaters
		'topology_array'  => \@{$options{-topology}},
		'topology_string' => exists $options{-topology_string} ? $options{-topology_string} : 0,
                'signal'          => exists $options{-signal} ? $options{-signal} : 0,
                'anchor'          => exists $options{-anchor} ? $options{-anchor} : 0,
		'n_term'          => exists $options{-n_terminal} ? $options{-n_terminal} : 'out', 
		'title'           => exists $options{-title} ? $options{-title} : '',
		'in'              => exists $options{-inside_label} ? $options{-inside_label} : "Cytoplasmic",
		'out'             => exists $options{-outside_label} ? $options{-outside_label} : "Extracellular",
		'membrane'        => exists $options{-membrane_label} ? $options{-membrane_label} : "Plasma Membrane",
	
		## dimensions
                'h_sig'              => exists $options{-h_sig} ? $options{-h_sig} : 15,
                'w_sig'              => exists $options{-w_sig} ? $options{-w_sig} : 15,
		'helix_height'       => exists $options{-helix_height} ? $options{-helix_height} : 130,
		'helix_width'        => exists $options{-helix_width} ? $options{-helix_width} : 50,
		'loop_width'         => exists $options{-loop_width} ? $options{-loop_width} : 20,
		'vertical_padding'   => exists $options{-vertical_padding} ? $options{-vertical_padding} : 140,
		'horizontal_padding' => exists $options{-horizontal_padding} ? $options{-horizontal_padding} : 150,
		'offset'             => exists $options{-membrane_offset} ? $options{-membrane_offset} : 6,
	
		## loop lengths and limits
		'short_loop'         => exists $options{-short_loop_height} ? $options{-short_loop_height} : 50,
		'medium_loop'        => exists $options{-medium_loop_height} ? $options{-medium_loop_height} : 60,
		'long_loop'          => exists $options{-long_loop_height} ? $options{-long_loop_height} : 70,
		'short_loop_limit'   => exists $options{-short_loop_limit} ? $options{-short_loop_limit} : 50,
		'long_loop_limit'    => exists $options{-long_loop_limit} ? $options{-long_loop_limit} : 70,
		'n_terminal_height'  => exists $options{-n_terminal_height} ? $options{-n_terminal_height} : 150,
		'c_terminal_height'  => exists $options{-c_terminal_height} ? $options{-c_terminal_height} : 80,
		'loop_heights'       => \%{$options{-loop_heights}},
		'n_term_offset'      => exists $options{-n_terminal_offset} ? $options{-n_terminal_offset} : 0,
		'c_term_offset'      => exists $options{-c_terminal_offset} ? $options{-c_terminal_offset} : 0,

		## colour scheme & display options
		'labels'             => exists $options{-show_labels} ? $options{-show_labels} : 'on',
		'bold_helices'       => exists $options{-bold_helices} ? $options{-bold_helices} : 1,
		'bold_labels'        => exists $options{-bold_labels} ? $options{-bold_labels} : 0,
		'scheme'             => exists $options{-colour_scheme} ? $options{-colour_scheme} : 'yellow',
		'draw_cytosol'       => exists $options{-draw_cytosol} ? $options{-draw_cytosol} : 0,
		'draw_bilayer'       => exists $options{-draw_bilayer} ? $options{-draw_bilayer} : 1,
		'draw_loops'         => exists $options{-draw_loops} ? $options{-draw_loops} : 1,
		'draw_terminai'      => exists $options{-draw_terminai} ? $options{-draw_terminai} : 1,
		'draw_helices'       => exists $options{-draw_helices} ? $options{-draw_helices} : 1,

		## labeling options
        	'loop_labels' => \%{$options{-labels}},
		'text_offset' => exists $options{-text_offset} ? $options{-text_offset} : 0,
		'helix_label' => exists $options{-helix_label} ? $options{-helix_label} : 'S',
		'n_term_label' => exists $options{'-n_term_label'} ? $options{-n_term_label} : 'N-Term',
		'c_term_label' => exists $options{'-c_term_label'} ? $options{-c_term_label} : 'C-Term',
		
		'dontsort' => exists $options{-dontsort} ? $options{-dontsort} : 0,
		'ttf_font' => exists $options{-ttf_font} ? $options{-ttf_font} : 0,
		'ttf_font_size' => exists $options{-ttf_font_size} ? $options{-ttf_font_size} : 8,

	};

	$self->{'loop_count'} = 1;

  	bless ($self,$class);
	
	return $self;

}

sub png {

 	my $self = shift;

	my @numeric = ('helix_height','helix_width','loop_width','vertical_padding','horizontal_padding','short_length','medium_loop_length','long_loop_length','short_loop_limit','long_loop_limit','n_terminal_height','membrane_offset','text_offset','n_term_offset','c_term_offset');

	foreach (@numeric){
		die "\nParameter $_ must be numeric.\n\n" if exists $self->{$_} && $self->{$_} =~ /-{?}\D+/;
	}

	foreach (keys %{$self->{'loop_labels'}}){
		die "\nLabel position $_ must be numeric.\n\n" if $_ =~ /\D+/;
	}

	foreach (keys %{$self->{'loop_heights'}}){
		die "\nLoop number $_ must be numeric.\n\n" if $_ =~ /\D+/;
	}

	foreach (values %{$self->{'loop_heights'}}){
		die "\nLoop height $_ must be numeric.\n\n" if $_ =~ /\D+/;
	}

	## n-terminal defaults to outside in it's not in,inside,out,outside
	$self->{'n_term'} = 'out' if (($self->{'n_term'} ne 'in')&&($self->{'n_term'} ne 'inside')&&($self->{'n_term'} ne 'out')||$self->{'n_term'} eq 'outside');
	$self->{'n_term'} = 'in' if $self->{'n_term'} eq 'inside';

	if ($self->{'topology_string'}){
	       	$self->{'topology_string'} =~ s/\d+\.//g;
		$self->{'topology_string'} =~ s/;/,/g;
		@{$self->{'topology_array'}} = split(/,/,$self->{'topology_string'});
	}

	## check to make sure we have pairs of helix boundaries and that data is numeric otherwise quit
	if (scalar @{$self->{'topology_array'}} % 2){
		die "\nUneven number of helix boundaries.\n\n";
	}

	foreach (@{$self->{'topology_array'}}){
		if ($_ =~ /\D/){
			die "\nTopology data is not numeric. $_\n\n";
			
		}
	}

	## check to make sure the TTF font exists, otherwise use gdSmallFont
	if ($self->{'ttf_font'}){
		unless (-e $self->{'ttf_font'}){
			print "\nCan't find font ".$self->{'ttf_font'}.".\n";
			$self->{'ttf_font'} = 0;
		}
	}

	my @sorted_topology = sort {$a <=> $b} @{$self->{'topology_array'}};
	
	## Don't automatically sort the topology array
	@sorted_topology = @{$self->{'topology_array'}} if $self->{'dontsort'};
	
	$self->{'helix_count'} = scalar @{$self->{'topology_array'}} / 2;
	
	unless ($self->{'helix_count'}){
		die "\nNo topology data found.\n\n";
	}

	## put helix start/stop points in $self->{'helix_span'} and loop lengths in $self->{'loop_length'}
	foreach (0..($self->{'helix_count'} - 1)){
		my $count = $_ * 2;
		$self->{'helix_span'}{$_ + 1}{'start'} = $sorted_topology[$count];
		$self->{'helix_span'}{$_ + 1}{'stop'} = $sorted_topology[$count + 1];
		$self->{'loop_length'}{$_ + 1} = scalar ($sorted_topology[$count + 2] - $sorted_topology[$count + 1]) unless ($_ + 1 == $self->{'helix_count'});
	}
	
	$self->{'width'} = ($self->{'horizontal_padding'} * 2) + ($self->{'helix_width'} * $self->{'helix_count'}) + ($self->{'loop_width'} * ($self->{'helix_count'} - 1));
	$self->{'height'} = $self->{'helix_height'} + ($self->{'vertical_padding'} * 2);
	
	## create a new image
	$self->{'im'} = new GD::Image($self->{'width'},$self->{'height'});
	
	$self->{'black'} = $self->{'im'}->colorAllocate(0,0,0);
	$self->{'white'} = $self->{'im'}->colorAllocate(255,255,255);
	$self->{'im'}->fill(0,0,$self->{'white'});
	
	## write title
	if ($self->{'ttf_font'}){
		$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,4,12,$self->{'title'},{linespacing=>0.6,charmap  => 'Unicode',}) if $self->{'title'};
	}else{
		$self->{'im'}->string(gdSmallFont,4,3,$self->{'title'},$self->{'black'}) if $self->{'title'};
	}

	$self->draw_cytosol  if $self->{'draw_cytosol'};
	$self->draw_bilayer  if $self->{'draw_bilayer'};
	$self->draw_loops    if $self->{'draw_loops'};
	$self->draw_terminai if $self->{'draw_terminai'};
        $self->draw_signal   if $self->{'signal'};
        $self->draw_anchor   if $self->{'anchor'};
	$self->draw_helices  if $self->{'draw_helices'};

	## use GD to convert to png
	return $self->{'im'}->GD::Image::png;
	
}


sub add_tmhmm_feat {

	my $self = shift;
	my $feat = shift;
	
	#print Dumper $feat;
	
	## add a helix from a tmhmm feature	
	if ($feat->{'_primary_tag'} eq 'transmembrane'){
		push @{$self->{'topology_array'}},$feat->{'_location'}{'_start'};
		push @{$self->{'topology_array'}},$feat->{'_location'}{'_end'};
	}

	## i've made a few changes to TmHmm.pm to include the inside/outside loops.
	## this bit looks for the topology of the 1st residue so we can now position the n-terminal 
	if ($feat->{'_location'}{'_start'} == 1){
		if ($feat->{'_primary_tag'} =~ /(\w+)_loop/){
			$self->{'n_term'} = $1; 
		}
	}
	
	return $self;

}

sub draw_cytosol {

  	my $self = shift;
	
	my $cytosol_offset = 5;
	my $light_grey = $self->{'im'}->colorAllocate(164,164,164);

	## draw cytosol
	$self->{'im'}->filledRectangle(($self->{'horizontal_padding'} / 3) ,($self->{'vertical_padding'} + $self->{'helix_height'} - $cytosol_offset),($self->{'width'} - ($self->{'horizontal_padding'} / 3)),($self->{'vertical_padding'} + ($self->{'helix_height'} * 2 ) + $cytosol_offset),$light_grey);

	return $self;
}



sub draw_bilayer {

  	my $self = shift;
	
	my $dark_grey  = $self->{'im'}->colorAllocate(40,40,40);
	my $dark_grey1 = $self->{'im'}->colorAllocate(50,50,50);
	my $dark_grey2 = $self->{'im'}->colorAllocate(60,60,60);
	my $dark_grey3 = $self->{'im'}->colorAllocate(70,70,70);

	## label either side of membrane
	if ($self->{'ttf_font'}){

		$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,($self->{'horizontal_padding'} / 3) + 2,($self->{'vertical_padding'} + $self->{'offset'} - 10),$self->{'out'},{linespacing=>0.6,charmap  => 'Unicode',}) if $self->{'labels'};
		$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,($self->{'horizontal_padding'} / 3) + 2,($self->{'vertical_padding'} + $self->{'helix_height'} - $self->{'offset'} + 22),$self->{'in'},{linespacing=>0.6,charmap  => 'Unicode',}) if $self->{'labels'};
	}else{
		$self->{'im'}->string(gdSmallFont,($self->{'horizontal_padding'} / 3) + 2,($self->{'vertical_padding'} + $self->{'offset'} - 24),$self->{'out'},$self->{'black'}) if $self->{'labels'};
		$self->{'im'}->string(gdSmallFont,($self->{'horizontal_padding'} / 3) + 2,($self->{'vertical_padding'} + $self->{'helix_height'} - $self->{'offset'} + 5),$self->{'in'},$self->{'black'}) if $self->{'labels'};
	}



	## draw membrane with graded fill
	$self->{'im'}->filledRectangle(($self->{'horizontal_padding'} / 3),($self->{'vertical_padding'} + $self->{'offset'}),($self->{'width'} - $self->{'horizontal_padding'} / 3),($self->{'vertical_padding'} + $self->{'helix_height'} - $self->{'offset'}),$dark_grey);
	$self->{'im'}->filledRectangle(($self->{'horizontal_padding'} / 3) + 1,($self->{'vertical_padding'} + $self->{'offset'}) + 1,($self->{'width'} - $self->{'horizontal_padding'} / 3) - 1,($self->{'vertical_padding'} + $self->{'helix_height'} - $self->{'offset'}) - 1,$dark_grey1);
	$self->{'im'}->filledRectangle(($self->{'horizontal_padding'} / 3) + 2,($self->{'vertical_padding'} + $self->{'offset'}) + 2,($self->{'width'} - $self->{'horizontal_padding'} / 3) - 2,($self->{'vertical_padding'} + $self->{'helix_height'} - $self->{'offset'}) - 2,$dark_grey2);
	$self->{'im'}->filledRectangle(($self->{'horizontal_padding'} / 3) + 3,($self->{'vertical_padding'} + $self->{'offset'}) + 3,($self->{'width'} - $self->{'horizontal_padding'} / 3) - 3,($self->{'vertical_padding'} + $self->{'helix_height'} - $self->{'offset'}) - 3,$dark_grey3);
			
	## label membrane
	if ($self->{'ttf_font'}){
		$self->{'im'}->stringFT($self->{'white'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,($self->{'horizontal_padding'} + ($self->{'helix_count'} * $self->{'helix_width'}) + (($self->{'helix_count'} - 1) * $self->{'loop_width'}) + 4) + 1,($self->{'vertical_padding'} + $self->{'helix_height'} - $self->{'offset'} - 3),$self->{'membrane'},{linespacing=>0.6,charmap  => 'Unicode',}) if $self->{'labels'};
	}else{
		$self->{'im'}->string(gdSmallFont,($self->{'horizontal_padding'} + ($self->{'helix_count'} * $self->{'helix_width'}) + (($self->{'helix_count'} - 1) * $self->{'loop_width'}) + 4) + 1,($self->{'vertical_padding'} + $self->{'helix_height'} - $self->{'offset'} - 14),$self->{'membrane'},$self->{'white'}) if $self->{'labels'};
	}

	
	return $self;

}


sub draw_anchor 
{

    my $self = shift;

    #signal anchor lime 
    my $colour  = $self->{'im'}->colorAllocate(5,240,0);
    my $colour1 = $self->{'im'}->colorAllocate(5,180,0);

    if ($self->{'n_term'} eq 'out')
    {
        $self->{'im'}->filledArc(
                                  $self->{'horizontal_padding'}  - $self->{'n_term_offset'}+6,
                                  $self->{'cy'} - $self->{'n_terminal_height'} / 2 ,
                                  $self->{'w_sig'}+2,
                                  $self->{'h_sig'}-4,
                                  0,
                                  360,
                                  $colour
                                 );

        $self->{'im'}->arc(
                                  $self->{'horizontal_padding'}  - $self->{'n_term_offset'}+6,
                                  $self->{'cy'} - $self->{'n_terminal_height'} / 2 ,
                                  $self->{'w_sig'}+2,
                                  $self->{'h_sig'}-4,
                                  0,
                                  360,
                                  $colour1
                                 );

       

        if ($self->{'ttf_font'})
        {

	    $self->{'im'}->stringFT(
                                     $self->{'black'},
                                     $self->{'ttf_font'},
                                     $self->{'ttf_font_size'},
                                     0,
                                     $self->{'horizontal_padding'}-$self->{'n_term_offset'}+$self->{'w_sig'},
                                     $self->{'cy'} - $self->{'n_terminal_height'} +4,
                                     'Anchor',
                                    {linespacing=>0.6, charmap  => 'Unicode'}
                                   );

	}else{
		$self->{'im'}->string(
                                      gdSmallFont,
                                      $self->{'horizontal_padding'}-$self->{'n_term_offset'}+$self->{'w_sig'},
                                      $self->{'cy'} - $self->{'n_terminal_height'}+4 ,
                                      'Anchor',
                                      $self->{'black'}
                                     );

	     }


    }
    else
    {
       
	$self->{'im'}->filledArc(
                                  $self->{'horizontal_padding'}  - $self->{'n_term_offset'}+6,
                                  $self->{'cy'}+$self->{'helix_height'}/2 - 5,
                                  $self->{'w_sig'}+2,
                                  $self->{'h_sig'}-4,
                                  0,
                                  360,
                                  $colour
                                );

       	$self->{'im'}->arc(
                                  $self->{'horizontal_padding'}  - $self->{'n_term_offset'}+6,
                                  $self->{'cy'}+$self->{'helix_height'}/2 - 5,
                                  $self->{'w_sig'}+2,
                                  $self->{'h_sig'}-4,
                                  0,
                                  360,
                                  $colour1
                                );

        
        if ($self->{'ttf_font'})
        {

	    $self->{'im'}->stringFT(
                                     $self->{'black'},
                                     $self->{'ttf_font'},
                                     $self->{'ttf_font_size'},
                                     0,
                                     $self->{'horizontal_padding'}-$self->{'n_term_offset'}+$self->{'w_sig'},
                                     $self->{'cy'} + $self->{'helix_height'} / 2,
                                     'Anchor',
                                    {linespacing=>0.6, charmap  => 'Unicode'}
                                   );

 	}else{
		$self->{'im'}->string(
                                      gdSmallFont,
                                      $self->{'horizontal_padding'}-$self->{'n_term_offset'}+$self->{'w_sig'},
                                      $self->{'cy'} + $self->{'helix_height'} / 2 ,
                                      'Anchor',
                                      $self->{'black'}
                                     );

	     }



    }

}

sub draw_signal 
{

    my $self = shift;

    #signal anchor red
    my $colour  = $self->{'im'}->colorAllocate(240,0,0);
    my $colour1 = $self->{'im'}->colorAllocate(180,0,0);
   

    if ($self->{'n_term'} eq 'out')
    {
        $self->{'im'}->filledArc($self->{'horizontal_padding'}  - $self->{'n_term_offset'}+5,
                                 $self->{'cy'} - $self->{'n_terminal_height'} / 2 ,
                                 $self->{'w_sig'},
                                 $self->{'h_sig'},
                                 0,
                                 360,
                                 $colour
                                );

        $self->{'im'}->arc($self->{'horizontal_padding'}  - $self->{'n_term_offset'}+5,
                                 $self->{'cy'} - $self->{'n_terminal_height'} / 2 ,
                                 $self->{'w_sig'},
                                 $self->{'h_sig'},
                                 0,
                                 360,
                                 $colour1
                                );

       	if ($self->{'ttf_font'})
        {

	    $self->{'im'}->stringFT(
                                    $self->{'black'},
                                    $self->{'ttf_font'},
                                    $self->{'ttf_font_size'},
                                    0,
                                    $self->{'horizontal_padding'},
                                    $self->{'cy'} - $self->{'n_terminal_height'} / 2,
                                    'Signal',
                                   {linespacing=>0.6, charmap  => 'Unicode'}
                                   );

	}else{
		$self->{'im'}->string(
                                      gdSmallFont,
                                       $self->{'horizontal_padding'}- $self->{'n_term_offfset'}+7+$self->{'w_sig'},
                                       $self->{'cy'} - $self->{'n_terminal_height'} / 2,
                                       'Signal',
                                       $self->{'black'}
                                     );

	    }


    }
    else
    {
	$self->{'im'}->filledArc(
                                  $self->{'horizontal_padding'}  - $self->{'n_term_offset'}+5,
                                  $self->{'cy'}+$self->{'helix_height'},
                                  $self->{'w_sig'},
                                  $self->{'h_sig'},
                                  0,
                                  360,
                                  $colour
                                );

	$self->{'im'}->arc      (
                                  $self->{'horizontal_padding'}  - $self->{'n_term_offset'}+5,
                                  $self->{'cy'}+$self->{'helix_height'},
                                  $self->{'w_sig'},
                                  $self->{'h_sig'},
                                  0,
                                  360,
                                  $colour1
                                 );

        if ($self->{'ttf_font'})
        {

	    $self->{'im'}->stringFT(
                                     $self->{'black'},
                                     $self->{'ttf_font'},
                                     $self->{'ttf_font_size'},
                                     0,
                                     $self->{'horizontal_padding'}-$self->{'n_term_offset'}+$self->{'w_sig'},
                                     $self->{'cy'} + $self->{'helix_height'} / 2,
                                     'Signal',
                                    {linespacing=>0.6, charmap  => 'Unicode'}
                                   );

 	}else{
		$self->{'im'}->string(
                                      gdSmallFont,
                                      $self->{'horizontal_padding'}-$self->{'n_term_offset'}+$self->{'w_sig'},
                                      $self->{'cy'} + $self->{'helix_height'} / 2 ,
                                      'Signal',
                                      $self->{'black'}
                                     );

	     }


    }
}

sub draw_helices {

	my $self = shift;

	my $x = $self->{'horizontal_padding'};
	my $y = $self->{'vertical_padding'};

	my ($colour,$colour1,$colour2,$colour3,$colour4,$colour5,$colour6);
		
	if($self->{'scheme'} eq 'blue'){

		$colour  = $self->{'im'}->colorAllocate(90,160,255);
		$colour1 = $self->{'im'}->colorAllocate(80,150,255);
		$colour2 = $self->{'im'}->colorAllocate(70,140,255);
		$colour3 = $self->{'im'}->colorAllocate(60,130,255);
		$colour4 = $self->{'im'}->colorAllocate(50,120,255);
		$colour5 = $self->{'im'}->colorAllocate(40,110,255);
		$colour6 = $self->{'im'}->colorAllocate(30,100,255);	
	
	}elsif($self->{'scheme'} eq 'pink'){

		$colour  = $self->{'im'}->colorAllocate(255,1,255);
		$colour1 = $self->{'im'}->colorAllocate(240,1,255);
		$colour2 = $self->{'im'}->colorAllocate(230,1,255);
		$colour3 = $self->{'im'}->colorAllocate(220,1,255);
		$colour4 = $self->{'im'}->colorAllocate(200,1,255);
		$colour5 = $self->{'im'}->colorAllocate(180,1,255);
		$colour6 = $self->{'im'}->colorAllocate(160,1,255);

	}elsif($self->{'scheme'} eq 'green'){

		$colour  = $self->{'im'}->colorAllocate(5,240,0);
		$colour1 = $self->{'im'}->colorAllocate(5,230,0);
		$colour2 = $self->{'im'}->colorAllocate(5,220,0);
		$colour3 = $self->{'im'}->colorAllocate(5,205,0);
		$colour4 = $self->{'im'}->colorAllocate(5,195,0);
		$colour5 = $self->{'im'}->colorAllocate(5,185,0);
		$colour6 = $self->{'im'}->colorAllocate(5,155,0);	
		
	}elsif($self->{'scheme'} eq 'red'){

		$colour  = $self->{'im'}->colorAllocate(240,0,0);
		$colour1 = $self->{'im'}->colorAllocate(230,0,0);
		$colour2 = $self->{'im'}->colorAllocate(220,0,0);
		$colour3 = $self->{'im'}->colorAllocate(205,0,0);
		$colour4 = $self->{'im'}->colorAllocate(190,0,0);
		$colour5 = $self->{'im'}->colorAllocate(170,0,0);
		$colour6 = $self->{'im'}->colorAllocate(150,0,0);
			
	}elsif($self->{'scheme'} eq 'white'){

		$colour  = $self->{'white'};
		$colour1 = $self->{'white'};
		$colour2 = $self->{'white'};
		$colour3 = $self->{'white'};
		$colour4 = $self->{'white'};
		$colour5 = $self->{'black'};
		$colour6 = $self->{'black'};
			
	}else{

		## default is yellow

		$colour  = $self->{'im'}->colorAllocate(255,235,55);
		$colour1 = $self->{'im'}->colorAllocate(255,230,50);
		$colour2 = $self->{'im'}->colorAllocate(255,220,40);
		$colour3 = $self->{'im'}->colorAllocate(255,210,30);
		$colour4 = $self->{'im'}->colorAllocate(255,200,20);
		$colour5 = $self->{'im'}->colorAllocate(255,190,10);
		$colour6 = $self->{'im'}->colorAllocate(255,180,0);
	}


	for (1..$self->{'helix_count'}){

		## draw helix, with graduated fill
		$self->{'im'}->filledRectangle($x,$y,($x + $self->{'helix_width'})-1,($y + $self->{'helix_height'})-1,$colour6);	
		$self->{'im'}->filledRectangle($x+1,$y+1,($x + $self->{'helix_width'})-1,($y + $self->{'helix_height'})-1,$colour5);
		$self->{'im'}->filledRectangle($x+2,$y+2,($x + $self->{'helix_width'})-2,($y + $self->{'helix_height'})-2,$colour4);
		$self->{'im'}->filledRectangle($x+3,$y+3,($x + $self->{'helix_width'})-3,($y + $self->{'helix_height'})-3,$colour3);
		$self->{'im'}->filledRectangle($x+4,$y+4,($x + $self->{'helix_width'})-4,($y + $self->{'helix_height'})-4,$colour2);
		$self->{'im'}->filledRectangle($x+5,$y+5,($x + $self->{'helix_width'})-5,($y + $self->{'helix_height'})-5,$colour1);
		$self->{'im'}->filledRectangle($x+6,$y+6,($x + $self->{'helix_width'})-6,($y + $self->{'helix_height'})-6,$colour);
	
		## draw a white box around it
		if ($self->{'bold_helices'}){
			$self->{'im'}->rectangle($x,$y,($x + $self->{'helix_width'}),($y + $self->{'helix_height'}),$self->{'black'});
			$self->{'im'}->rectangle($x - 1,$y - 1,($x + $self->{'helix_width'} + 1),($y + $self->{'helix_height'} + 1),$self->{'white'});
		}else{
			$self->{'im'}->rectangle($x,$y,($x + $self->{'helix_width'}),($y + $self->{'helix_height'}),$self->{'white'});
		}
		
		## this is the text on each helix
		my $text = substr($self->{'helix_label'},0,1).$_;

		## draw a white box in the centre and label the helix
		my $x_offset = 5;
		$x_offset = 8 if $_ >= 10;
	
		my $white_box = 12;
		$white_box = 17 if $_ >= 10;

		$self->{'im'}->filledRectangle(($x + ($self->{'helix_width'} / 2) - $x_offset) - 5,($y + ($self->{'helix_height'} / 2) - 7) - 3,$white_box + ($x + ($self->{'helix_width'} / 2) - $x_offset) + 3,12 + ($y + ($self->{'helix_height'} / 2) - 4),$colour1) if $self->{'labels'};
		$self->{'im'}->filledRectangle(($x + ($self->{'helix_width'} / 2) - $x_offset) - 4,($y + ($self->{'helix_height'} / 2) - 7) - 2,$white_box + ($x + ($self->{'helix_width'} / 2) - $x_offset) + 2,12 + ($y + ($self->{'helix_height'} / 2) - 5),$colour2) if $self->{'labels'};
		$self->{'im'}->filledRectangle(($x + ($self->{'helix_width'} / 2) - $x_offset) - 3,($y + ($self->{'helix_height'} / 2) - 7) - 1,$white_box + ($x + ($self->{'helix_width'} / 2) - $x_offset) + 1,12 + ($y + ($self->{'helix_height'} / 2) - 6),$colour3) if $self->{'labels'};
		
		
		$self->{'im'}->filledRectangle(($x + ($self->{'helix_width'} / 2) - $x_offset) - 2,($y + ($self->{'helix_height'} / 2) - 7),$white_box + ($x + ($self->{'helix_width'} / 2) - $x_offset),12 + ($y + ($self->{'helix_height'} / 2) - 7),$self->{'white'}) if $self->{'labels'};

		if ($self->{'ttf_font'}){
			$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,($x + ($self->{'helix_width'} / 2) - $x_offset - 1),($y + ($self->{'helix_height'} / 2) + 4),$text,{linespacing=>0.6,charmap  => 'Unicode',}) if $self->{'labels'};
		}else{
			$self->{'im'}->string(gdSmallFont,($x + ($self->{'helix_width'} / 2) - $x_offset),($y + ($self->{'helix_height'} / 2) - 7),$text,$self->{'black'}) if $self->{'labels'};
		}

		## label start and end positions of helices
		
		$self->{'x'} = $x;
		$self->{'y'} = $y;
		
		if ($self->{'labels'}){
		
			if (($self->{'n_term'} eq 'out')&&($_ % 2)){
				
				$self->label_helix_o_i();
				
			}elsif($self->{'n_term'} eq 'out'){
				
				$self->label_helix_i_o();
				
			}elsif(($self->{'n_term'} eq 'in')&&($_ % 2)){
				
				$self->label_helix_i_o();
				
			}else{
				
				$self->label_helix_o_i();
				
			}
		}
		
		$x = $self->{'horizontal_padding'} + ($_ * ($self->{'helix_width'} + $self->{'loop_width'}));
	
	}

	if ($self->{'labels'}){

		$x = $self->{'horizontal_padding'};
		$y = $self->{'vertical_padding'};

		for (1..$self->{'helix_count'}){
	
			my $y_mod = 0;
			
			foreach my $l (sort {$b <=> $a} keys %{$self->{'loop_labels'}}){

				if (($l >= $self->{'helix_span'}{$_}{'start'})&&($l <= $self->{'helix_span'}{$_}{'stop'})){

					my $label_length = 0;
					if ($self->{'ttf_font'}){
						## Might need to fiddle with this
						my $size_dif = $self->{'ttf_font_size'} - 8;
						$label_length = 6 + (6 * (length $self->{'loop_labels'}{$l}) + (8 * $size_dif));
					}else{					
						$label_length = 9 + (6 * length $self->{'loop_labels'}{$l});
					}
				
					if ($_ % 2){

						if ($self->{'bold_labels'}){
							$self->{'im'}->filledRectangle(($x + ($self->{'helix_width'} / 2)) + 5,($y + ($self->{'helix_height'} / 2) - 30) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)) + 2,12 + ($y + ($self->{'helix_height'} / 2) - 24) + $y_mod,$self->{'black'});
							my $p = new GD::Polygon;
        						$p->addPt(($x + ($self->{'helix_width'} / 2)) + 4,($y + ($self->{'helix_height'} / 2) - 30) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)));
        						$p->addPt(($x + ($self->{'helix_width'} / 2)) - 5,($y + ($self->{'helix_height'} / 2) - 21) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)));
       							$p->addPt(($x + ($self->{'helix_width'} / 2)) + 4,($y + ($self->{'helix_height'} / 2) - 12) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)));
        						$self->{'im'}->filledPolygon($p,$self->{'black'});
						}

						## add darker box
						$self->{'im'}->filledRectangle(($x + ($self->{'helix_width'} / 2)) + 6,($y + ($self->{'helix_height'} / 2) - 29) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)) + 1,12 + ($y + ($self->{'helix_height'} / 2) - 25) + $y_mod,$colour6);
						## add white box
						$self->{'im'}->filledRectangle(($x + ($self->{'helix_width'} / 2)) + 7,($y + ($self->{'helix_height'} / 2) - 28) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)),12 + ($y + ($self->{'helix_height'} / 2) - 26) + $y_mod,$self->{'white'});


						## draw darker arrowhead
				        	my $poly = new GD::Polygon;
        					$poly->addPt(($x + ($self->{'helix_width'} / 2)) + 5,($y + ($self->{'helix_height'} / 2) - 29) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)));
        					$poly->addPt(($x + ($self->{'helix_width'} / 2)) - 3,($y + ($self->{'helix_height'} / 2) - 21) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)));
       						$poly->addPt(($x + ($self->{'helix_width'} / 2)) + 5,($y + ($self->{'helix_height'} / 2) - 13) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)));
        					$self->{'im'}->filledPolygon($poly,$colour6);

						## draw white arrowhead
				        	my $poly2 = new GD::Polygon;
        					$poly2->addPt(($x + ($self->{'helix_width'} / 2)) + 6,($y + ($self->{'helix_height'} / 2) - 28) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)));
        					$poly2->addPt(($x + ($self->{'helix_width'} / 2)) - 1,($y + ($self->{'helix_height'} / 2) - 21) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)));
       						$poly2->addPt(($x + ($self->{'helix_width'} / 2)) + 6,($y + ($self->{'helix_height'} / 2) - 14) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)));
        					$self->{'im'}->filledPolygon($poly2,$self->{'white'});

						## add label
						if ($self->{'ttf_font'}){
							$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,($x + ($self->{'helix_width'} / 2)) + 10,($y + ($self->{'helix_height'} / 2) - 16) + $y_mod,$self->{'loop_labels'}{$l},{linespacing=>0.6,charmap  => 'Unicode',});
						}else{
							$self->{'im'}->string(gdSmallFont,($x + ($self->{'helix_width'} / 2)) + 9,($y + ($self->{'helix_height'} / 2) - 27) + $y_mod,$self->{'loop_labels'}{$l},$self->{'black'});
						}
			
						$y_mod = $y_mod - 19;	
							
					}else{

						if ($self->{'bold_labels'}){
							$self->{'im'}->filledRectangle(($x + ($self->{'helix_width'} / 2)) + 5,($y + ($self->{'helix_height'} / 2) + 9) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)) + 2,13 + ($y + ($self->{'helix_height'} / 2) + 15) + $y_mod,$self->{'black'});
				        		my $b = new GD::Polygon;
        						$b->addPt(($x + ($self->{'helix_width'} / 2)) + 4,($y + ($self->{'helix_height'} / 2) + 10) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)));
        						$b->addPt(($x + ($self->{'helix_width'} / 2)) - 5,($y + ($self->{'helix_height'} / 2) + 19) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)));
       							$b->addPt(($x + ($self->{'helix_width'} / 2)) + 4,($y + ($self->{'helix_height'} / 2) + 28) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)));
        						$self->{'im'}->filledPolygon($b,$self->{'black'});
						}

						## add darker box
						$self->{'im'}->filledRectangle(($x + ($self->{'helix_width'} / 2)) + 6,($y + ($self->{'helix_height'} / 2) + 11) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)) + 1,12 + ($y + ($self->{'helix_height'} / 2) + 15) + $y_mod,$colour6);

						## add white box
						$self->{'im'}->filledRectangle(($x + ($self->{'helix_width'} / 2)) + 7,($y + ($self->{'helix_height'} / 2) + 12) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)),12 + ($y + ($self->{'helix_height'} / 2) + 14) + $y_mod,$self->{'white'});

						## draw darker arrowhead
				        	my $poly = new GD::Polygon;
        					$poly->addPt(($x + ($self->{'helix_width'} / 2)) + 5,($y + ($self->{'helix_height'} / 2) + 11) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)));
        					$poly->addPt(($x + ($self->{'helix_width'} / 2)) - 3,($y + ($self->{'helix_height'} / 2) + 19) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)));
       						$poly->addPt(($x + ($self->{'helix_width'} / 2)) + 5,($y + ($self->{'helix_height'} / 2) + 27) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)));
        					$self->{'im'}->filledPolygon($poly,$colour6);

						## draw white arrowhead
				        	my $poly2 = new GD::Polygon;
        					$poly2->addPt(($x + ($self->{'helix_width'} / 2)) + 6,($y + ($self->{'helix_height'} / 2) + 12) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)));
        					$poly2->addPt(($x + ($self->{'helix_width'} / 2)) - 1,($y + ($self->{'helix_height'} / 2) + 19) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)));
       						$poly2->addPt(($x + ($self->{'helix_width'} / 2)) + 6,($y + ($self->{'helix_height'} / 2) + 26) + $y_mod,$label_length + ($x + ($self->{'helix_width'} / 2)));
        					$self->{'im'}->filledPolygon($poly2,$self->{'white'});

						## add label
						if ($self->{'ttf_font'}){
							$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,($x + ($self->{'helix_width'} / 2)) + 10,($y + ($self->{'helix_height'} / 2) + 24) + $y_mod,$self->{'loop_labels'}{$l},{linespacing=>0.6,charmap  => 'Unicode',});
						}else{
							$self->{'im'}->string(gdSmallFont,($x + ($self->{'helix_width'} / 2)) + 9,($y + ($self->{'helix_height'} / 2) + 12) + $y_mod,$self->{'loop_labels'}{$l},$self->{'black'});
						}

						$y_mod = $y_mod + 19;	
					}
				}
			}
		
			$x = $self->{'horizontal_padding'} + ($_ * ($self->{'helix_width'} + $self->{'loop_width'}));
		}
	}

	return $self;
}

sub draw_terminai {

	my $self = shift;
	
	my $loop_number = ($self->{'helix_count'} - 1);

	## width of terminal
	$self->{'w'}  = $self->{'helix_width'} + $self->{'loop_width'};
	$self->{'cx'} = $self->{'horizontal_padding'} - ($self->{'loop_width'} / 2);
	$self->{'cy'} = $self->{'vertical_padding'};

	## draw N-terminal
	if ($self->{'n_term'} eq 'out'){

		$self->{'im'}->arc(($self->{'cx'} - $self->{'n_term_offset'}),
                                    $self->{'cy'},($self->{'w'} + (2 * $self->{'n_term_offset'})),
                                    $self->{'n_terminal_height'},270,360,$self->{'black'});	

		if ($self->{'ttf_font'}){
			$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,
                                                ($self->{'horizontal_padding'} - ($self->{'w'} / 2) - 33 - $self->{'n_term_offset'}),
                                                ($self->{'cy'} - ($self->{'n_terminal_height'} / 2) + 5),
                                                 $self->{'n_term_label'},
                                                 {linespacing=>0.6,charmap  => 'Unicode',}) if $self->{'labels'};
		}else{
			$self->{'im'}->string(gdSmallFont,
                                             ($self->{'horizontal_padding'} - ($self->{'w'} / 2) - 40 - $self->{'n_term_offset'}) +20,
                                             ($self->{'cy'} - ($self->{'n_terminal_height'} / 2) - 6),
                                              $self->{'n_term_label'},$self->{'black'}) if $self->{'labels'};
		    }
		
		
		
		## label n-terminal
		my $y_mod = 0;
		foreach (sort {$b <=> $a} keys %{$self->{'loop_labels'}}){
			if ($_ <= $self->{'helix_span'}{1}{'start'}){
					
				if ($self->{'ttf_font'}){
					$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,($self->{'horizontal_padding'} - ($self->{'w'} / 2) - 33 - $self->{'n_term_offset'}),($self->{'cy'} - ($self->{'n_terminal_height'} / 2) - 6) - 4 + $y_mod,$self->{'loop_labels'}{$_},{linespacing=>0.6,charmap  => 'Unicode',}) if $self->{'labels'};
				}else{
					$self->{'im'}->string(gdSmallFont,($self->{'horizontal_padding'} - ($self->{'w'} / 2) - 40 - $self->{'n_term_offset'}),($self->{'cy'} - ($self->{'n_terminal_height'} / 2) - 6) - 15 + $y_mod,$self->{'loop_labels'}{$_},$self->{'black'}) if $self->{'labels'};
				}

				$y_mod = $y_mod - 15;
			}
		}

	}else{
		$self->{'cy'} = $self->{'cy'} + $self->{'helix_height'};
		$self->{'im'}->arc(($self->{'cx'} - $self->{'n_term_offset'}),$self->{'cy'},($self->{'w'} + (2 * $self->{'n_term_offset'})),$self->{'n_terminal_height'},0,90,$self->{'black'});
		
		if ($self->{'ttf_font'}){
			$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,
                                               ($self->{'horizontal_padding'} - ($self->{'w'} / 2) - 33 - $self->{'n_term_offset'}),
                                               ($self->{'cy'} + ($self->{'n_terminal_height'} / 2) + 5),
                                                $self->{'n_term_label'},{linespacing=>0.6,charmap  => 'Unicode',}) if $self->{'labels'};
		}else{
			$self->{'im'}->string(gdSmallFont,
                                             ($self->{'horizontal_padding'} - ($self->{'w'} / 2) - 40 - $self->{'n_term_offset'}) + 20,
                                             ($self->{'cy'} + ($self->{'n_terminal_height'} / 2) - 6),
                                              $self->{'n_term_label'},$self->{'black'}) if $self->{'labels'};
		}

		## label n-terminal
		my $y_mod = 0;
		foreach (sort {$a <=> $b} keys %{$self->{'loop_labels'}}){
			if ($_ <= $self->{'helix_span'}{1}{'start'}){
				
			if ($self->{'ttf_font'}){
				$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,($self->{'horizontal_padding'} - ($self->{'w'} / 2) - 33 - $self->{'n_term_offset'}),($self->{'cy'} + ($self->{'n_terminal_height'} / 2) - 6) + 26 + $y_mod,$self->{'loop_labels'}{$_},{linespacing=>0.6,charmap  => 'Unicode',}) if $self->{'labels'};
			}else{
				$self->{'im'}->string(gdSmallFont,($self->{'horizontal_padding'} - ($self->{'w'} / 2) - 40 - $self->{'n_term_offset'}),($self->{'cy'} + ($self->{'n_terminal_height'} / 2) - 6) + 15 + $y_mod,$self->{'loop_labels'}{$_},$self->{'black'}) if $self->{'labels'};
			}	

				$y_mod = $y_mod + 15;
			}
		}

	}

	$self->{'cx'} = ($self->{'helix_count'} * $self->{'helix_width'}) + (($self->{'helix_count'} - 1) * $self->{'loop_width'}) + $self->{'horizontal_padding'} + ($self->{'loop_width'} / 2);

	## draw C-terminal
	if (($self->{'n_term'} eq 'out')&&($loop_number % 2)){
		
		$self->draw_ext_c_term;

	}elsif($self->{'n_term'} eq 'out'){
		
		$self->draw_int_c_term;

	}elsif(($self->{'n_term'} eq 'in')&&($loop_number % 2)){
		
		$self->draw_int_c_term;

	}elsif($self->{'n_term'} eq 'in'){

		$self->draw_ext_c_term;

	}

	return $self;
}

sub draw_loops {

	my $self = shift;
	
	$self->{'x'} = $self->{'horizontal_padding'} + ($self->{'helix_width'} / 2);
	$self->{'h'} = $self->{'medium_loop'};
	$self->{'w'} = $self->{'helix_width'} + $self->{'loop_width'};

	for (1..($self->{'helix_count'} - 1)){

		## Alter loop height according to its actual length
		if ($self->{'loop_length'}{$_} < $self->{'short_loop_limit'}){
			$self->{'h'} = $self->{'short_loop'};
		}elsif($self->{'loop_length'}{$_} > $self->{'long_loop_limit'}){
			$self->{'h'} = $self->{'long_loop'};
		}

		$self->{'l_start'} = $self->{'helix_span'}{$_}{'stop'};
		$self->{'l_stop'} = $self->{'helix_span'}{$_ + 1}{'start'};

		if (($self->{'n_term'} eq 'out')&&($_ % 2)){

			$self->draw_int_loop;

		}elsif($self->{'n_term'} eq 'out'){
			$self->draw_ext_loop;

		}elsif(($self->{'n_term'} eq 'in')&&($_ % 2)){
			$self->draw_ext_loop;

		}elsif($self->{'n_term'} eq 'in'){
			
			$self->draw_int_loop;

		}
		
		$self->{'x'} = $self->{'x'} + $self->{'helix_width'} + $self->{'loop_width'};
		$self->{'x'} = $self->{'horizontal_padding'} if $_ == ($self->{'helix_count'} - 1);
		
	}
	
	return $self;
}

sub label_helix_o_i {

	my $self = shift;

	my $offset = 6;
	$offset = 3 if $self->{'helix_span'}{$_}{'start'} < 100;
	if ($self->{'ttf_font'}){
		$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,($self->{'x'} + ($self->{'helix_width'} / 2) - $offset) - 1,$self->{'y'} + 13,$self->{'helix_span'}{$_}{'start'},{linespacing=>0.6,charmap  => 'Unicode',});
	}else{
		$self->{'im'}->string(gdSmallFont,($self->{'x'} + ($self->{'helix_width'} / 2) - $offset) - 1,$self->{'y'} + 1,$self->{'helix_span'}{$_}{'start'},$self->{'black'});
	}	

	$offset = 3 if $self->{'helix_span'}{$_}{'stop'} < 100;
	$offset = 6 if $self->{'helix_span'}{$_}{'stop'} >= 100;

	if ($self->{'ttf_font'}){
		$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,($self->{'x'} + ($self->{'helix_width'} / 2) - $offset) - 1,($self->{'y'} + $self->{'helix_height'} - 3),$self->{'helix_span'}{$_}{'stop'},{linespacing=>0.6,charmap  => 'Unicode',});
	}else{
		$self->{'im'}->string(gdSmallFont,($self->{'x'} + ($self->{'helix_width'} / 2) - $offset) - 1,($self->{'y'} + $self->{'helix_height'} - 14),$self->{'helix_span'}{$_}{'stop'},$self->{'black'});
	}	
	
	return $self;
}

sub label_helix_i_o {

	my $self = shift;

	my $offset = 6;
	$offset = 3 if $self->{'helix_span'}{$_}{'stop'} < 100;
	if ($self->{'ttf_font'}){
		$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,($self->{'x'} + ($self->{'helix_width'} / 2) - $offset) - 1,$self->{'y'} + 13,$self->{'helix_span'}{$_}{'stop'},{linespacing=>0.6,charmap  => 'Unicode',});
	}else{
		$self->{'im'}->string(gdSmallFont,($self->{'x'} + ($self->{'helix_width'} / 2) - $offset) - 1,$self->{'y'} + 1,$self->{'helix_span'}{$_}{'stop'},$self->{'black'});
	}	

	$offset = 3 if $self->{'helix_span'}{$_}{'start'} < 100;

	if ($self->{'ttf_font'}){
		$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,($self->{'x'} + ($self->{'helix_width'} / 2) - $offset) - 1,($self->{'y'} + $self->{'helix_height'} - 3),$self->{'helix_span'}{$_}{'start'},{linespacing=>0.6,charmap  => 'Unicode',});
	}else{
		$self->{'im'}->string(gdSmallFont,($self->{'x'} + ($self->{'helix_width'} / 2) - $offset) - 1,($self->{'y'} + $self->{'helix_height'} - 14),$self->{'helix_span'}{$_}{'start'},$self->{'black'});
	}	
	
	return $self;
}

sub draw_int_c_term {

	my $self = shift;
	
	## draw internal c-terminal
	$self->{'im'}->arc(($self->{'cx'} + $self->{'c_term_offset'}),($self->{'vertical_padding'} + $self->{'helix_height'}),($self->{'w'} + (2 * $self->{'c_term_offset'})),$self->{'c_terminal_height'},90,180,$self->{'black'});	

	if ($self->{'ttf_font'}){
		$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,($self->{'cx'} + 4 + $self->{'c_term_offset'}),(($self->{'vertical_padding'} + $self->{'helix_height'}) + ($self->{'c_terminal_height'} / 2) + 5),$self->{'c_term_label'},{linespacing=>0.6,charmap  => 'Unicode',}) if $self->{'labels'};
	}else{
		$self->{'im'}->string(gdSmallFont,($self->{'cx'} + 3 + $self->{'c_term_offset'}),(($self->{'vertical_padding'} + $self->{'helix_height'}) + ($self->{'c_terminal_height'} / 2) - 6),$self->{'c_term_label'},$self->{'black'}) if $self->{'labels'};
	}	

	## label terminal
	if ($self->{'labels'}){
		my $y_mod = 0;
		foreach (sort {$a <=> $b} keys %{$self->{'loop_labels'}}){
			if ($_ >= $self->{'helix_span'}{$self->{'helix_count'}}{'stop'}){
							
				if ($self->{'ttf_font'}){
					$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,($self->{'cx'} + 4 + $self->{'c_term_offset'}),(($self->{'vertical_padding'} + $self->{'helix_height'}) + ($self->{'c_terminal_height'} / 2) - 6) + 26 + $y_mod,$self->{'loop_labels'}{$_},{linespacing=>0.6,charmap  => 'Unicode',}) if $self->{'labels'};
				}else{
					$self->{'im'}->string(gdSmallFont,($self->{'cx'} + 3 + $self->{'c_term_offset'}),(($self->{'vertical_padding'} + $self->{'helix_height'}) + ($self->{'c_terminal_height'} / 2) - 6) + 15 + $y_mod,$self->{'loop_labels'}{$_},$self->{'black'});
				}
				$y_mod = $y_mod + 15;
			}
		}
	}
	return $self;
}

sub draw_ext_c_term {

	my $self = shift;
	
	## draw external c-terminal
	$self->{'im'}->arc(($self->{'cx'} + $self->{'c_term_offset'}),$self->{'vertical_padding'},($self->{'w'} + (2 * $self->{'c_term_offset'})),$self->{'c_terminal_height'},180,270,$self->{'black'});
	
	if ($self->{'ttf_font'}){
		$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,,($self->{'cx'} + 3 + $self->{'c_term_offset'}),($self->{'vertical_padding'} - ($self->{'c_terminal_height'} / 2) + 5),$self->{'c_term_label'},{linespacing=>0.6,charmap  => 'Unicode',}) if $self->{'labels'};
	}else{
		$self->{'im'}->string(gdSmallFont,($self->{'cx'} + 3 + $self->{'c_term_offset'}),($self->{'vertical_padding'} - ($self->{'c_terminal_height'} / 2) - 6),$self->{'c_term_label'},$self->{'black'}) if $self->{'labels'};
	}


	## label terminal
	if ($self->{'labels'}){
		my $y_mod = 0;
		foreach (sort {$b <=> $a} keys %{$self->{'loop_labels'}}){
			if ($_ >= $self->{'helix_span'}{$self->{'helix_count'}}{'stop'}){

				if ($self->{'ttf_font'}){
					$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,($self->{'cx'} + 3 + $self->{'c_term_offset'}),($self->{'vertical_padding'} - ($self->{'c_terminal_height'} / 2) - 6) - 4 - $y_mod,$self->{'loop_labels'}{$_},{linespacing=>0.6,charmap  => 'Unicode',}) if $self->{'labels'};
				}else{
					$self->{'im'}->string(gdSmallFont,($self->{'cx'} + 3 + $self->{'c_term_offset'}),($self->{'vertical_padding'} - ($self->{'c_terminal_height'} / 2) - 6) - 15 - $y_mod,$self->{'loop_labels'}{$_},$self->{'black'});
				}
				$y_mod = $y_mod + 15;
			}
		}
	}
	return $self;
}

sub draw_int_loop {

	my $self = shift;
	
	## draw internal loop
	$self->{'cx'} = $self->{'x'} + ($self->{'helix_width'} / 2) + ($self->{'loop_width'} / 2);

	## this sets the height to the value given by the loop_heights hash
	foreach (sort {$a <=> $b} keys %{$self->{'loop_heights'}}){
		$self->{'h'} = $self->{'loop_heights'}{$_} if $_ == $self->{'loop_count'};
	}

	$self->{'im'}->arc($self->{'cx'},($self->{'vertical_padding'} + $self->{'helix_height'}),$self->{'w'},$self->{'h'},0,180,$self->{'black'});

	## label loop
	if ($self->{'labels'}){
		my $y_mod = 0;
		foreach (sort {$a <=> $b} keys %{$self->{'loop_labels'}}){
			if (($_ >= $self->{'l_start'})&&($_ <= $self->{'l_stop'})){
				$self->{'im'}->line($self->{'cx'},($self->{'vertical_padding'} + $self->{'helix_height'} + ($self->{'h'} / 2)),$self->{'cx'},($self->{'vertical_padding'} + $self->{'helix_height'} + ($self->{'h'} / 2) + 5),$self->{'black'}) unless $y_mod;

				if ($self->{'ttf_font'}){
					$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,($self->{'cx'} + $self->{'text_offset'} - 3),($self->{'vertical_padding'} + $self->{'helix_height'} + ($self->{'h'} / 2) + 17) + $y_mod,$self->{'loop_labels'}{$_},{linespacing=>0.6,charmap  => 'Unicode',}) if $self->{'labels'};
				}else{
					$self->{'im'}->string(gdSmallFont,($self->{'cx'} + $self->{'text_offset'}),($self->{'vertical_padding'} + $self->{'helix_height'} + ($self->{'h'} / 2) + 6) + $y_mod,$self->{'loop_labels'}{$_},$self->{'black'}) ;
				}
				$y_mod = $y_mod + 15;
			}
		}
	}
	$self->{'loop_count'}++;
	return $self;
}

sub draw_ext_loop {

	my $self = shift;
	
	## draw external loop
	$self->{'cx'} = $self->{'x'} + ($self->{'helix_width'} / 2) + ($self->{'loop_width'} / 2);
	
	## this sets the height to the value given by the loop_heights hash
	foreach (sort {$a <=> $b} keys %{$self->{'loop_heights'}}){
		$self->{'h'} = $self->{'loop_heights'}{$_} if $_ == $self->{'loop_count'};
	}
	
	$self->{'im'}->arc($self->{'cx'},$self->{'vertical_padding'},$self->{'w'},$self->{'h'},180,360,$self->{'black'});

	## label loop
	if ($self->{'labels'}){
		my $y_mod = 0;
		foreach (sort {$b <=> $a} keys %{$self->{'loop_labels'}}){
			if (($_ >= $self->{'l_start'})&&($_ <= $self->{'l_stop'})){
				$self->{'im'}->line($self->{'cx'},($self->{'vertical_padding'} - ($self->{'h'} / 2)),$self->{'cx'},($self->{'vertical_padding'} - ($self->{'h'} / 2) - 5),$self->{'black'}) unless $y_mod;

				if ($self->{'ttf_font'}){
					$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,($self->{'cx'} + $self->{'text_offset'} - 3),($self->{'vertical_padding'} - ($self->{'h'} / 2) - 8) + $y_mod,$self->{'loop_labels'}{$_},{linespacing=>0.6,charmap  => 'Unicode',}) if $self->{'labels'};
				}else{
					$self->{'im'}->string(gdSmallFont,($self->{'cx'} + $self->{'text_offset'}),($self->{'vertical_padding'} - ($self->{'h'} / 2) - 19) + $y_mod,$self->{'loop_labels'}{$_},$self->{'black'});
				}

				$y_mod = $y_mod - 15;
			}	
		}
	}
	$self->{'loop_count'}++;
	return $self;
}

1;

=head1 NAME

LOCAL::lib::DrawMembrane.pm - draw a cartoon of an Alpha-helical transmembrane protein.

=head1 SYNOPSIS

use DrawTransmbrane;

my @topology = (20,45,59,70,86,109,145,168,194,220);

## Simple use - -topology is the only option that is required

my $im = DrawTransmembrane->new(-title => 'This is a cartoon displaying transmembrane helices.',
                                -topology => \@topology);

## More advanced use

my %labels = (5 => '5 - Sulphation Site',
              21 => '1st Helix',
              47 => '40 - Mutation',
              60 => 'Voltage Sensor',
              72 => '72 - Mutation 2',
              73 => '73 - Mutation 3',
              138 => '138 - Glycosylation Site',
              170 => '170 - Phosphorylation Site',
              200 => 'Last Helix');

my $im = DrawTransmembrane->new(-n_terminal=> 'out',
                                -topology => \@topology,
                                -bold_helices=> 1,
                                -labels=> \%labels,
                                -text_offset=> -15,
                                -outside_label=>'Lumen',
                                -inside_label=>'Cytoplasm',
                                -membrane_label=>'Membrane',
                                -vertical_padding=> 155);


## Parse Tmhmm data

use Bio::Tools::Tmhmm;

my $im = DrawTransmembrane->new(-title=>'Let\'s parse some Tmhmm output...',
                                -bold_helices=> 1);

open(FILE, 'tmhmm.out');  

my $parser = new Bio::Tools::Tmhmm(-fh => \*FILE );

while(my $tmhmm_feat = $parser->next_result ) {

	## Load features into DrawTransmembrane object
	$im->add_tmhmm_feat($tmhmm_feat);

}

close FILE;

## Now write the image to a .png file

open(OUTPUT, ">output.png");

binmode OUTPUT;

print OUTPUT $im->png; 

close OUTPUT;

=head1 DESCRIPTION

A module to draw a cartoon of an alpha-helical transmembrane protein. It uses GD and allows the image to be written to a .png file.

The options are a set of tag/value pairs as follows:

  Option              Value                                         Default
  ------              -----                                         -------

  -topology           Array containing transmembrane helix          none
	              boundaries. This is the only option that 
		      is required

  -topology_string    Alternative to -topology, provide a string    none
                      containing the topology data in the form
		      A.11,31;B.41,59;C.86,107;D.145,166
	    
  -n_terminal         Location of the N-terminal of the sequence,   out
  	              either 'in' or 'out'

  -title              Title to add to the image                     none
  
  -inside_label       Label for the inside of the membrane          Cytoplasmic

  -outside_label      Label for the outside of the membrane         Extracellular

  -membrane_label     Label for the membrane                        Plasma Membrane

  -colour_scheme      Colour scheme to use. Current choices are     blue
                      blue, yellow, red, green, pink or white. 
		    
  -labels             Label loops and helices using data from a     none
                      hash, e.g.
		      
		      %labels = (138 => 'Glycosylation Site',
		                 190 => 'Binding Site');
      
		      The hash key must be numeric, ranges are 
		      not allowed.
		      
  -bold_helices       Draws black boxes round helices               1
  
  -bold_labels        Draws black boxes round labels                0

  -text_offset        Shift the text labeling the loops. Use a      0 
                      negative value to shift it left, a positive
		      value to shift it right

  -helix_height       Transmembrane helix height                    130

  -helix_width        Transmembrane helix width                     50

  -loop_width         Loop width                                    20

  -vertical_padding   Vertical padding                              140

  -horizontal_padding Horizontal Padding                            150

  -membrane_offset    Offest between helix end and membrane         6

  -short_loop_height  Height of short loops                         90
  
  -medium_loop_height Height of medium loops                        120
  
  -long_loop_height   Height of long loops                          150
  
  -short_loop_limit   Length in residues below which a loop is      15
                      classed as short
  
  -long_loop_limit    Length in residues above which a loop is      30
                      classed as long
  
  -loop_heights       Explicitly set heights of each loop, e.g.
                    
                      %loop_heights = (1 => 45,
                                       2 => 220,
                                       3 => 50,
                                       4 => 220,
                                       9 => 70);
				     
                      The key corresponds to the loop number. Both
                      key and value must be numeric. If you use
                      -loop_height and there is a defined height
                      for the current loop then other height values
                      will be overridden		      				     
  
  -n_terminal_height  Height of N-terminal                          150
  
  -c_terminal_height  Height of C-terminal                          80

  -n_terminal_offset  Shift the N-terminal left by this amount      0
		      
  -c_terminal_offset  Shift the C-terminal right by this amount     0

  -helix_label        Change the 'S' label on each helix. Only 1    S
                      character is allowed

  -show_labels        Display text labels                           on
		      
  -draw_cytosol       Show the cytosol                              false

  -draw_bilayer       Show the membrane                             true
  
  -draw_loops         Show the loops                                true

  -draw_terminai      Show the terminai                             true

  -draw_helices       Show the helices                              true

  -dontsort           Don't automatically sort the topology array   0

  -ttf_font           Path to TTF font, e.g.                        none 
                      /usr/share/fonts/msttcorefonts/arial.ttf
		      
  -ttf_font_size      Default size for TTF font. Use 7-9 with       8
                      Arial for best results  
		
		
Height, width, padding and other numerical values can gernerally be left alone. They are useful if your labels consists of a lot of text as this may lead to them overlapping. In this case try increasing the loop_width or helix_width options. -text_offset is also very useful for avoiding overlapping.

=head1 AUTHOR

Tim Nugent E<lt>timnugent@gmail.comE<gt>

=cut
