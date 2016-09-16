package Schematic;

use strict;
use warnings;
use GD;
our @ISA = "GD"; 
our @EXPORT = qw(new png);
our $VERSION = '1.13';

sub new {

	my $class = shift;
 	my %options = @_;
	
	my $self = {

		## general paramaters
		'feat_hash'          => $options{-feat_hash},
		'feat_string'        => exists $options{-feat_string}        ? $options{-feat_string}    : 0,
		'title'              => exists $options{-title}              ? $options{-title}          : '',
                'len'                => exists $options{-len}                ? $options{-len}            : undef,
	
		## dimensions
		'max_height'         => exists $options{-max_height}         ? $options{-max_height}         : 700,
		'max_width'          => exists $options{-max_width}          ? $options{-max_width}          : 900,
		'vertical_padding'   => exists $options{-vertical_padding}   ? $options{-vertical_padding}   : 140,
		'horizontal_padding' => exists $options{-horizontal_padding} ? $options{-horizontal_padding} : 150,
		'offset'             => exists $options{-offset}             ? $options{-offset}             : 6,
                'ss_height'          => exists $options{-ss_height}          ? $options{-ss_height}          : 15,
	
		## colour scheme & display options
		'labels'             => exists $options{-show_labels}        ? $options{-show_labels}        : 0,
		'bold_labels'        => exists $options{-bold_labels}        ? $options{-bold_labels}        : 0,
		'scheme'             => exists $options{-colour_scheme}      ? $options{-colour_scheme}      : 'default',


		## labeling options
		'text_offset'        => exists $options{-text_offset}        ? $options{-text_offset}        : 0,
		'n_term_label'       => exists $options{-n_term_label}       ? $options{-n_term_label}       : 'N-Terminal',
		'c_term_label'       => exists $options{-c_term_label}       ? $options{-c_term_label}       : 'C-Terminal',
		
                ## font settings 
                'fonts'              => exists $options{-fonts}              ? $options{-fonts}              : 0,
		'ttf_font'           => defined($options{-ttf_font})         ? $options{-ttf_font}           : 0,
		'ttf_font_size'      => exists $options{-ttf_font_size}      ? $options{-ttf_font_size}      : 8,

	};

  	bless ($self,$class);
	
	return $self;

}

sub png 
{

 	my $self = shift;

	my @numeric = ('len','vertical_padding','horizontal_padding','text_offset','n_term_offset','c_term_offset');

	foreach (@numeric)
        {
		die "\nParameter $_ must be numeric.\n\n" if exists $self->{$_} && $self->{$_} =~ /-{?}\D+/;
	}

        ## convert string of features into array
	if ($self->{'feat_string'})
        {
	       	$self->{'feat_string'}   =~ s/^\d+\.//g;
		$self->{'feat_string'}   =~ s/;/,/g;
	      @{$self->{'feat_array'}}   =  split(/,/,$self->{'feat_string'});
	}

	## check to see if feat array data is numeric otherwise quit
	foreach (keys %{$self->{'feat_hash'}})
        {
	    foreach my $feat ( @{$self->{'feat_hash'}->{$_}} )
            {
		if ($feat->{'from'} =~ /\D/ || $feat->{'to'} =~ /\D/)
                {
			die "\nFeature data is not numeric. $_ $feat->{'f'} $feat->{'t'}\n\n";
			
		}
	    }
	}

	## check to make sure the TTF font exists, otherwise use gdSmallFont
	if ($self->{'ttf_font'})
	{
	#	print "$self->{fonts}/$self->{ttf_font}.ttf\n";
		if(-e "$self->{fonts}/$self->{ttf_font}.ttf")
		{
			#print "Font found\n";
		}
		else
		{
			print "\nCan't find font ".$self->{ttf_font}.".\n";
			$self->{'ttf_font'} = 0;
		}
		#unless (! -e "$self->{fonts}/$self->{ttf_font}.ttf")
        #{
        #	print "\nCan't find font ".$self->{ttf_font}.".\n";
		#	$self->{'ttf_font'} = 0;
		#}
                 
	}
	
	unless (keys %{$self->{'feat_hash'}} > 0)
        {
		die "\nNo feature data found.\n\n";
	}


	$self->{'phos'} = exists($self->{'feat_hash'}->{PH}) ? 1 : 0;
        $self->{'glyc'} = exists($self->{'feat_hash'}->{OG}) ? 1 : 0;
	$self->{'lowc'} = exists($self->{'feat_hash'}->{LC}) ? 1 : 0;
        $self->{'diso'} = exists($self->{'feat_hash'}->{DI}) ? 1 : 0;
        $self->{'pest'} = exists($self->{'feat_hash'}->{PE}) ? 1 : 0;
        $self->{'coils'}= exists($self->{'feat_hash'}->{CC}) ? 1 : 0;

        my $height = ($self->{'glyc'}*50+
                      $self->{'lowc'}*20+
                      $self->{'phos'}*50+
                      $self->{'pest'}*20+
                      $self->{'coils'}*20+
                      $self->{'diso'}*20)+100;

         

        $self->{'max_height'} = $height  <  $self->{'max_height'} ? $height : $self->{'max_height'};

	## create a new image
	$self->{'im'} = new GD::Image($self->{'max_width'},$self->{'max_height'});
	
        $self->{'cy'}    = $self->{'glyc'}*50+$self->{'phos'}*50+20;
        
	$self->{'black'} =  $self->{'im'}->colorAllocate(0,0,0);
        $self->{'green'} =  $self->{'im'}->colorAllocate(0,100,0);
        $self->{'lime'}  =  $self->{'im'}->colorAllocate(0,180,0);
        $self->{'brown'} =  $self->{'im'}->colorAllocate(165,42,42);
	$self->{'white'} =  $self->{'im'}->colorAllocate(255,255,255);
        $self->{'gray'}  =  $self->{'im'}->colorAllocate(200,200,200); 
        

	$self->{'im'}->fill(0,0,$self->{'white'});
        $self->{'secstr'}= exists($self->{'feat_hash'}->{SS}) ? 1 : 0;

	## write title
	if ($self->{'ttf_font'})
        {
		$self->{'im'}->stringFT(
                                         $self->{'black'},
                                         $self->{'ttf_font'},
                                         $self->{'ttf_font_size'},
                                         0,
                                         4,
                                         12,
                                         $self->{'title'},
                                         {linespacing=>0.6,charmap  => 'Unicode',}
                                        ) if $self->{'title'};
	}else{
		$self->{'im'}->string(
                                      gdSmallFont,
                                      4,
                                      3,
                                      $self->{'title'},
                                      $self->{'black'}
                                     ) if $self->{'title'};
	}

	## draw baseline for features
        $self->{ny}=0;
        $self->draw_baseline();
        $self->draw_secstr();
        $self->draw_disorder() if $self->{'diso'};
        $self->draw_signal()   if $self->{'anchor'} || $self->{'signal'};
        $self->draw_lowc()     if $self->{'lowc'};
        $self->draw_pest()     if $self->{'pest'};
        $self->draw_coils()    if $self->{'coils'};
        $self->draw_ptm()      if $self->{'phos'} || $self->{'glyc'}; 
        $self->draw_key();

	## use GD to convert to png
	return $self->{'im'}->GD::Image::png;
	
}

sub draw_key
{
       my $self = shift;

       my $y    = $self->{'max_height'} - 15;
       my $x    = 50;

       my ($hcolor,$ecolor);

	$hcolor->[0] = $self->{'im'}->colorAllocate(24,116,205);
	$hcolor->[1] = $self->{'im'}->colorAllocate(34,146,255);
	$hcolor->[2] = $self->{'im'}->colorAllocate(30,136,245);
	$hcolor->[3] = $self->{'im'}->colorAllocate(30,126,235);



	$ecolor->[0] = $self->{'im'}->colorAllocate(205,38,38);
	$ecolor->[1] = $self->{'im'}->colorAllocate(245,8,8);
	$ecolor->[2] = $self->{'im'}->colorAllocate(235,8,8);
	$ecolor->[3] = $self->{'im'}->colorAllocate(235,18,18);



       #color key for secondary structure

        my $width=50; my $height=15;

       if($self->{'secstr'})
       { 
           my $color=$hcolor;

           $self->{'im'}->filledRectangle(
                                            $x-20,
                                            $y+7,
                                            $x+$width+20,
                                            $y+9,
					    $self->{black}
		 			    );
	

           $self->{'im'}->filledRectangle(
                                            $x,
                                            $y,
                                            $x+$width,
                                            $y+$height,
                                            $color->[0]
		 			    );
	
            
             
            $self->{'im'}->filledArc      (
                                            $x+$width,
                                            $y+$height/2,
                                            7,
                                            $height,
                                            270,
                                            90,
                                            $color->[0]
		 			    );
  
            for(my $i=1; $i < 4; $i++)
            {
            $self->{'im'}->filledRectangle(
                                           $x,
                                           $y+$i-1,
                                           $x+$width+$i-1,
                                           $y+$i,
                                           $color->[$i]    
                                          );
            }

            
            $self->{'im'}->filledArc      (
                                            $x,
                                            $y+$height/2,
                                            5,
                                            $height,
                                            270,
                                            90,
					    $self->{'white'}
		 			    );
            
            $self->{im}->stringFT(
                                  $self->{'black'},
                                  "$self->{fonts}/$self->{'ttf_font'}.ttf",
                                  $self->{ttf_font_size},
                                  0,
                                  $x+80,
                                  $y+13,
                                  "Predicted helix",
                                  {linespacing => 0.8, charmap => 'Unicode'}
                             ); 

            $color=$ecolor;
            $x+= 400;

             my $poly = new GD::Polygon;

                $poly->addPt($x,$y+3);
                $poly->addPt($x+$width-5,$y+3);
                $poly->addPt($x+$width-5,$y);
                $poly->addPt($x+$width,$y+$height/2);
                $poly->addPt($x+$width-5,$y+$height);
                $poly->addPt($x+$width-5,$y+$height-3);
                $poly->addPt($x,$y+$height-3);

                
                $self->{'im'}->filledRectangle(
                                               $x-20,
                                               $y+7,
                                               $x+$width+20,
                                               $y+9,
				  	       $self->{black}
		 			    );


                $self->{'im'}->filledPolygon(
                                            $poly,
                                            $color->[3]
                                            );

                $self->{'im'}->openPolygon(
                                            $poly,
                                            $color->[0]
                                           );

           
                $self->{im}->stringFT(
                                  $self->{'black'},
                                  "$self->{fonts}/$self->{'ttf_font'}.ttf",
                                  $self->{ttf_font_size},
                                  0,
                                  $x+80,
                                  $y+13,
                                  "Predicted sheet",
                                  {linespacing => 0.8, charmap => 'Unicode'}
                             ); 
   
       }
}

sub draw_coils
{
	my $self = shift;

	my $y = $self->{'cy'}+$self->{'ny'};
        my $height = 5;

	my $color = $self->{'im'}->colorAllocate(128,205,255);


 
        foreach my $ss (@{$self->{'feat_hash'}->{'CC'}})
        {
            my $width = int (($ss->{to} - $ss->{'from'}) *  (($self->{'max_width'}-$self->{'horizontal_padding'}*2)/$self->{'len'}));
            my $x     = $self->{'horizontal_padding'} + 
                        int ($ss->{'from'} *  (($self->{'max_width'}-$self->{'horizontal_padding'}*2)/$self->{'len'}));
           
            $self->{'im'}->filledRectangle(
                                            $x,
                                            $y,
                                            $x+$width,
                                            $y+$height,
                                            $color
		 			    );


            $self->{'im'}->filledArc      (
                                            $x+$width,
                                            $y+2,
                                            3,
                                            $height,
                                            270,
                                            90,
                                            $color
		 			    );
	
        }

        $self->{im}->stringFT(
                               $self->{'black'},
                               "$self->{fonts}/$self->{'ttf_font'}.ttf",
                               $self->{ttf_font_size},
                               0,
                               $self->{max_width}-$self->{horizontal_padding}-70,
                               $y+$height+10,
                               "Coiled Coils",
                               {linespacing => 0.8, charmap => 'Unicode'}
                             ); 

	return $self;

}

sub draw_pest
{
	my $self = shift;

	my $y = $self->{'cy'}+$self->{'ny'};
        my $height = 5;

	my ($bcolor,$ycolor);

        $bcolor = $self->{'im'}->colorAllocate(228,162,244);
	$ycolor = $self->{'im'}->colorAllocate(228,112,214);

 
        foreach my $ss (@{$self->{'feat_hash'}->{'PE'}})
        {
            my $width = int (($ss->{to} - $ss->{'from'}) *  (($self->{'max_width'}-$self->{'horizontal_padding'}*2)/$self->{'len'}));
            my $x     = $self->{'horizontal_padding'} + 
                        int ($ss->{'from'} *  (($self->{'max_width'}-$self->{'horizontal_padding'}*2)/$self->{'len'}));
           
            my $color = defined($ss->{'val'}) && $ss->{'val'} > 0 ? $ycolor : $bcolor;

            $self->{'im'}->filledRectangle(
                                            $x,
                                            $y,
                                            $x+$width,
                                            $y+$height,
                                            $color
		 			    );


            $self->{'im'}->filledArc      (
                                            $x+$width,
                                            $y+2,
                                            3,
                                            $height,
                                            270,
                                            90,
                                            $color
		 			    );
	
             
        }

        $self->{im}->stringFT(
                               $self->{'black'},
                               "$self->{fonts}/$self->{'ttf_font'}.ttf",
                               $self->{ttf_font_size},
                               0,
                               $self->{max_width}-$self->{horizontal_padding}-35,
                               $y+$height+12,
                               "PEST",
                               {linespacing => 0.8, charmap => 'Unicode'}
                             ); 

        $self->{'ny'} += 20;

	return $self;

}

sub draw_lowc
{
	my $self = shift;

	my $y = $self->{'cy'}+$self->{'ny'};
        my $height = 5;

	my ($bcolor,$ycolor);

	$ycolor = $self->{'im'}->colorAllocate(255,165,0);
 
        foreach my $ss (@{$self->{'feat_hash'}->{'LC'}})
        {
            my $width = int (($ss->{to} - $ss->{'from'}) *  (($self->{'max_width'}-$self->{'horizontal_padding'}*2)/$self->{'len'}));
            my $x     = $self->{'horizontal_padding'} + 
                        int ($ss->{'from'} *  (($self->{'max_width'}-$self->{'horizontal_padding'}*2)/$self->{'len'}));
           
            $self->{'im'}->filledRectangle(
                                            $x,
                                            $y,
                                            $x+$width,
                                            $y+$height,
                                            $ycolor
		 			    );
	
             
         }

        $self->{im}->stringFT(
                               $self->{'black'},
                               "$self->{fonts}/$self->{'ttf_font'}.ttf",
                               $self->{ttf_font_size},
                               0,
                               $self->{max_width}-$self->{horizontal_padding}-90,
                               $y+$height+10,
                               "Low complexity",
                               {linespacing => 0.8, charmap => 'Unicode'}
                             ); 
        $self->{'ny'} +=20;

	return $self;

}

sub draw_ptm
{
  my $self = $_[0];

  ##-- draw letter P for phos, G for glyc N for nglyc --##
  ##-- draw border for P G and N                      --##

  

  my ($x,$y) = (0,0,0,0);

  my $midx=0;
  my $basey = $self->{'cy'}-20;
  
  my $scale=($self->{'max_width'}-($self->{'horizontal_padding'}*2))/$self->{'len'};

  ##draw second baseline ##

  $self->{'im'}->filledRectangle(
                                 $self->{horizontal_padding},
                                 $self->{cy}-20,
                                 $self->{max_width}-$self->{horizontal_padding},
                                 $self->{cy}-19,
                                 $self->{'black'}   
                                );

  my $maxy = 35;

  if($self->{'phos'})
  {

  
  $self->{'im'}->line(
                        $self->{horizontal_padding},
                        $basey-($maxy*0.5),
                        $self->{max_width}-$self->{horizontal_padding},
                        $basey-($maxy*0.5),
                        $self->{'gray'}   
                      );

  $self->{'im'}->line(
                        $self->{horizontal_padding},
                        $basey-$maxy,
                        $self->{max_width}-$self->{horizontal_padding},
                        $basey-$maxy,
                        $self->{'gray'}   
                      );

  $self->{'im'}->line(
                        $self->{horizontal_padding},
                        $basey,
                        $self->{horizontal_padding},
                        $basey-$maxy,
                        $self->{'gray'}   
                      );

  $self->{'im'}->line(
                        $self->{max_width}-$self->{horizontal_padding},
                        $basey,
                        $self->{max_width}-$self->{horizontal_padding},
                        $basey-$maxy,
                        $self->{'gray'}   
                      );

  foreach my $p(@ {$self->{feat_hash}->{PH}})
  {
     $x = $p->{'from'}*$scale+$self->{horizontal_padding};
     $y = $p->{'val'}*$maxy;

     $self->{'im'}->filledRectangle(
                                    $x,
                                    $basey-$y,
                                    $x+1,
                                    $basey,
                                    $self->{'brown'}
                                   );
            
  }

	#print $self->{ttf_font}."\n";
  $self->{'im'}->stringFT(
                        $self->{'black'},
                        "$self->{fonts}/$self->{ttf_font}.ttf",
                        $self->{ttf_font_size},
                        0,
                        $self->{max_width}-$self->{horizontal_padding}-90,
                        $basey+15,
                        "Phosphorylation",
                        {linespacing => 0.8, charmap=>"Unicode" }
                       ) or die "can't use ttf $self->{tt_font} $!\n";  
 


  $self->{'im'}->stringFT(
                        $self->{'black'},
                        "$self->{fonts}/$self->{ttf_font}.ttf",
                        $self->{ttf_font_size}-2,
                        0,
                        $self->{max_width}-$self->{horizontal_padding}+4,
                        $basey+2,
                        "0",
                        {linespacing => 0.8, charmap=>"Unicode" }
                       ) or die "can't use ttf $self->{tt_font} $!\n";  
 


  $self->{'im'}->stringFT(
                        $self->{'black'},
                        "$self->{fonts}/$self->{ttf_font}.ttf",
                        $self->{ttf_font_size}-2,
                        0,
                        $self->{max_width}-$self->{horizontal_padding}+4,
                        $basey-$maxy,
                        "1",
                        {linespacing => 0.8, charmap=>"Unicode" }
                       ) or die "can't use ttf $self->{tt_font} $!\n";  
 
  $basey -= 55;

  }

  if($self->{'glyc'})
  {

  $self->{'im'}->line(
                        $self->{horizontal_padding},
                        $basey-$maxy,
                        $self->{max_width}-$self->{horizontal_padding},
                        $basey-$maxy,
                        $self->{'gray'}   
                      );

  
   
  $self->{'im'}->line(
                        $self->{horizontal_padding},
                        $basey-($maxy*0.5),
                        $self->{max_width}-$self->{horizontal_padding},
                        $basey-($maxy*0.5),
                        $self->{'gray'}   
                      );

   $self->{'im'}->line(
                        $self->{horizontal_padding},
                        $basey,
                        $self->{horizontal_padding},
                        $basey-$maxy,
                        $self->{'gray'}   
                      );

   $self->{'im'}->line(
                        $self->{max_width}-$self->{horizontal_padding},
                        $basey,
                        $self->{max_width}-$self->{horizontal_padding},
                        $basey-$maxy,
                        $self->{'gray'}   
                      );

  $self->{'im'}->filledRectangle(
                                 $self->{horizontal_padding},
                                 $basey,
                                 $self->{max_width}-$self->{horizontal_padding},
                                 $basey+1,
                                 $self->{black}
                                );     

 
     foreach my $o(@ {$self->{feat_hash}->{OG}})
     {
       $x= $o->{'from'}*$scale+$self->{horizontal_padding};
       $y= $o->{'val'}*$maxy;
   
       $self->{'im'}->filledRectangle(
                                      $x,
                                      $basey-$y,
                                      $x+1,
                                      $basey,
                                      $self->{'green'} 
                                     );

  }

   $self->{'im'}->stringFT(
                        $self->{'black'},
                        "$self->{fonts}/$self->{'ttf_font'}.ttf",
                        $self->{ttf_font_size}-2,
                        0, 
                        $self->{max_width}-$self->{horizontal_padding}+4,
                        $basey+2,
                        "0" 
                       );  

  $self->{'im'}->stringFT(
                        $self->{'black'},
                        "$self->{fonts}/$self->{'ttf_font'}.ttf",
                        $self->{ttf_font_size}-2,
                        0,
                        $self->{max_width}-$self->{horizontal_padding}+4,
                        $basey-$maxy,
                        "1" 
                       );  


  $self->{'im'}->stringFT(
                        $self->{'black'},
                        "$self->{fonts}/$self->{ttf_font}.ttf",
                        $self->{ttf_font_size},
                        0,
                        $self->{max_width}-$self->{horizontal_padding}-80,
                        $basey+15,
                        "Glycosylation",
                        {linespacing=>0.8, charmap => 'Unicode'}
                       );

 

  foreach my $n(@ {$self->{feat_hash}->{NG}})
  {
    $x= $n->{'from'}*$scale+$self->{horizontal_padding};
    $y= $n->{'val'}*$maxy;

   $self->{'im'}->filledRectangle(
                                  $x,
                                  $basey-$y,
                                  $x+1,
                                  $basey,
                                  $self->{'lime'}
                                 );
  }  

  
  }
}

sub draw_baseline
{

     my $self = $_[0];

        $self->{'im'}->filledRectangle(
                                       $self->{'horizontal_padding'},
                                       $self->{'cy'}+8,
                                       $self->{'max_width'}-$self->{'horizontal_padding'},
                                       $self->{'cy'}+11,
                                       $self->{'black'}
                                      );
         $self->{'im'}->stringFT(
                                     $self->{'black'},
                                     "$self->{fonts}/$self->{'ttf_font'}.ttf",
                                     $self->{'ttf_font_size'}-1,
                                     0,
                                     $self->{'horizontal_padding'},
                                     $self->{'max_height'}-25,
                                     '1',
                                    {linespacing=>0.6, charmap  => 'Unicode'}
                                   );

         $self->{'im'}->stringFT(
                                     $self->{'black'},
                                     "$self->{fonts}/$self->{'ttf_font'}.ttf",
                                     $self->{'ttf_font_size'}-1,
                                     0,
                                     $self->{'max_width'}-$self->{'horizontal_padding'}-15,
                                     $self->{'max_height'}-25,
                                     $self->{'len'},
                                    {linespacing=>0.6, charmap  => 'Unicode'}
                                   );

	
  

     my $scale=($self->{'max_width'}-($self->{'horizontal_padding'}*2))/$self->{'len'};

     $self->{'dash_height'} = $self->{'max_height'}-20;
  
     my $step = int($self->{'len'} / 900) * 50 + 50;

     for (my $i=$step; $i < $self->{'len'}-($step/2) ; $i+=$step)
     {
         	$self->{'im'}->stringFT(
                                      $self->{'black'},
                                      "$self->{fonts}/$self->{'ttf_font'}.ttf",
                                      $self->{ttf_font_size}-1,
                                      0,
                                      $self->{'horizontal_padding'}+($i*$scale)-8,
                                      $self->{'max_height'}-25,
                                      $i,
                                      {linespacing=>0.8, charmap => 'Unicode'}
                                     );
       
                #--  draw grey dashed lines  --#

                $self->{'im'}->setStyle(
                                        $self->{'white'},$self->{'white'},
                                        $self->{'gray'} ,$self->{'gray'}
                                       );

                $self->{'im'}->line(
                                    $self->{'horizontal_padding'}+($i*$scale),
                                    10,
                                    $self->{'horizontal_padding'}+($i*$scale),
                                    $self->{'dash_height'}-20,
                                    gdStyled                     
                                   );	     
     }

     $self->{ny} = 10;
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



sub draw_secstr
{

	my $self = shift;

	my $y = $self->{'cy'}-$self->{ss_height}/2 + 10;
        my $height = $self->{ss_height};

	my ($hcolor,$ecolor);

	$hcolor->[0] = $self->{'im'}->colorAllocate(24,116,205);
	$hcolor->[1] = $self->{'im'}->colorAllocate(34,146,255);
	$hcolor->[2] = $self->{'im'}->colorAllocate(30,136,245);
	$hcolor->[3] = $self->{'im'}->colorAllocate(30,126,235);



	$ecolor->[0] = $self->{'im'}->colorAllocate(205,38,38);
	$ecolor->[1] = $self->{'im'}->colorAllocate(245,8,8);
	$ecolor->[2] = $self->{'im'}->colorAllocate(235,8,8);
	$ecolor->[3] = $self->{'im'}->colorAllocate(235,18,18);




 
        foreach my $ss (@{$self->{'feat_hash'}->{'SS'}})
        {
            my $width = int (($ss->{to} - $ss->{'from'}) *  (($self->{'max_width'}-$self->{'horizontal_padding'}*2)/$self->{'len'}));
            my $x     = $self->{'horizontal_padding'} + 
                        int ($ss->{from} *  (($self->{'max_width'}-$self->{'horizontal_padding'}*2)/$self->{'len'}));

            my $color = $ss->{type} =~ /H/ ? $hcolor : $ecolor;
       
            if( $ss->{type} =~ /H/ )
            {
    
            $self->{'im'}->filledRectangle(
                                            $x,
                                            $y,
                                            $x+$width,
                                            $y+$height,
                                            $color->[0]
		 			    );
	
            
             
            $self->{'im'}->filledArc      (
                                            $x+$width,
                                            $y+$height/2,
                                            7,
                                            $height,
                                            270,
                                            90,
                                            $color->[0]
		 			    );
  
            for(my $i=1; $i < 4; $i++)
            {
            $self->{'im'}->filledRectangle(
                                           $x,
                                           $y+$i-1,
                                           $x+$width+$i-1,
                                           $y+$i,
                                           $color->[$i]    
                                          );
            }

	    $self->{'im'}->filledArc      (
                                            $x,
                                            $y+$height/2,
                                            5,
                                            $height,
                                            270,
                                            90,
					    $self->{'white'}
		 			    );
            
          }
           else
               {
                 # draw arrows for sheets
                my $poly = new GD::Polygon;
                   $poly->addPt($x,$y+3);
                   $poly->addPt($x+$width-5,$y+3);
                   $poly->addPt($x+$width-5,$y);
                   $poly->addPt($x+$width,$y+$height/2);
                   $poly->addPt($x+$width-5,$y+$height);
                   $poly->addPt($x+$width-5,$y+$height-3);
                   $poly->addPt($x,$y+$height-3);

                 $self->{'im'}->filledPolygon(
                                            $poly,
                                            $color->[3]
                                            );

                $self->{'im'}->openPolygon(
                                            $poly,
                                            $color->[0]
                                           );


               }



        }


        $self->{'im'}->stringFT(
                                 $self->{'black'},
                                 "fonts/$self->{ttf_font}.ttf",
                                 $self->{ttf_font_size},
                                 0,
                                 $self->{'max_width'}-$self->{'horizontal_padding'}-60,
                                 $self->{cy}+30,
                                 "Structure", 
                                 {linespacing=>0.8, charmap => 'Unicode'}

                               ); 

        $self->{'ny'} += 25;
	return $self;
}

sub draw_disorder
{
	my $self = shift;

	my $y = $self->{'cy'}+$self->{'ny'};
        my $height = 5;

	my ($bcolor,$ycolor);

	$bcolor = $self->{'im'}->colorAllocate(138,58,58);
	$ycolor = $self->{'im'}->colorAllocate(255,236,139);

 
        foreach my $ss (@{$self->{'feat_hash'}->{'DI'}})
        {
            my $width = int (($ss->{to} - $ss->{'from'}) *  (($self->{'max_width'}-$self->{'horizontal_padding'}*2)/$self->{'len'}));
            my $x     = $self->{'horizontal_padding'} + 
                        int ($ss->{'from'} *  (($self->{'max_width'}-$self->{'horizontal_padding'}*2)/$self->{'len'}));
           
            $self->{'im'}->filledRectangle(
                                            $x,
                                            $y,
                                            $x+$width,
                                            $y+$height,
                                            $ycolor
		 			    );
	
             
            $self->{'im'}->filledArc      (
                                            $x+$width,
                                            $y+2,
                                            3,
                                            $height,
                                            270,
                                            90,
                                            $ycolor
		 			    );

            
            
        }

        $self->{im}->stringFT(
                               $self->{'black'},
                               "fonts/$self->{'ttf_font'}.ttf",
                               $self->{ttf_font_size},
                               0,
                               $self->{max_width}-$self->{horizontal_padding}-55,
                               $y+$height+10,
                               "Disorder",
                               {linespacing => 0.8, charmap => 'Unicode'}
                             ); 


        $self->{'ny'} += 20;
	return $self;

}

1;

