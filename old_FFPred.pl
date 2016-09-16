#!/usr/bin/perl -w

# Runs FFPred3 jobs. NOTE: this is an edited web server version!
# 
# The version number must be specified by $version below.
# The subversion (date of SVM training) must be specified by $subversion below.
# 
# The user may need to modify $FFPred_dir as described below.
# 
# Also, if the temporary files produced by FFPred are not needed, they can be
# deleted automatically by removing the comment sign (#) from the appropriate
# line - look for the "Use the following line to automatically delete all
# temporary files" message below.
# 
# Currently takes just one input text file at a time ('FASTA_input_file' - see
# $usage below), containing one or more protein sequences in FASTA format.

use strict;
use Getopt::Long;
use Cwd;
use File::Copy;



my $version = '3';
my $subversion = '5Feb2015';

# IMPORTANT : set the following variable ($FFPred_dir) to the absolute path to
# the main FFPred directory.
my $FFPred_dir = '/webdata/binaries/current/FFPred3'; # Edited in the web server version.

# Start here.

my $usage = "\n$0 -i [FASTA_input_file] -o [output_folder]".          # Edited in the web server version.
            "\n$0 -a -i [FASTA_input_folder] -o [output_folder]\n\n"; # Edited in the web server version.
# ./FFPred.pl -i in/TEST.fsa -o out # Sample usage (check that input files/folders exist!).
# ./FFPred.pl -a -i in -o out # Sample usage for array job (check that input files/folders exist!).

my ($arrayjob, $fasta, $dirOut);

my $args = GetOptions(
                       "a"   => \$arrayjob,
                       "i=s" => \$fasta,
                       "o=s" => \$dirOut
                     );

die "$usage" unless (defined($fasta));
#$dirOut = "$FFPred_dir/executables/out" unless (defined($dirOut)); # Skipped in the web server version.

my $currentWD = cwd();
foreach my $path ($fasta, $dirOut)
{
    $path = "$currentWD/$path" unless ((substr($path, 0, 1) eq '/') || (substr($path, 0, 1) eq '~'));
}

if (defined($arrayjob))
{
    die "No job-array variable found in the environment !\n" unless (defined($ENV{'SGE_TASK_ID'}));
    my @jobarray = sort glob("$fasta/*"); # BEWARE! if somebody touches FASTA_input_folder while jobs are spawned ...
    $fasta = $jobarray[$ENV{'SGE_TASK_ID'}-1];
}

#`mkdir -p $dirOut` unless (-d $dirOut); # Skipped in the web server version.
#die "$usage" unless ((-T $fasta) && (-w $dirOut)); # Skipped in the web server version.

# Initialise useful variables.
my $submit_datetime = Timestamp();

#my $jobs_directory = ReadConfigPATH($FFPred_dir)            # Skipped in the web server version.
#    or die "\n----- error in sub ReadConfigPATH -----\n\n"; # Skipped in the web server version.
#system("mkdir -p $jobs_directory") == 0                                                                                      # Skipped in the web server version.
#    or die "ABORT - could not create directory '$jobs_directory' for temporary files as indicated on the CONFIG file: $!\n"; # Skipped in the web server version.
#chdir $jobs_directory; # This may be useful, as BLAST default error.log file will be written here. # Skipped in the web server version.

my $GOdomains = {
                  'cellular_component' => 'CC',
                  'molecular_function' => 'MF',
                  'biological_process' => 'BP'
                };

# This list must be identical to the one that was used for training the current
# version/subversion of FFPred. The list used for training can be found in the
# 'Train_FFPred_SVMs.pm' module within the training material (usually on the CS
# cluster, within the "FFPred" project space, subfolder "material/lib/") - see
# the definitions at the beginning of the module.
# Note that this list is not included in the downloadable standalone script.
my $feature_groups_lists = {
                             'COILS'            =>  [1..12] ,
                             'DISOPRED'         =>  [13..30] ,
                             'LOWC'             =>  [31..42] ,
                             'MEMSAT'           =>  [43..57] ,
                             'NETNGLYC'         =>  [58..68] ,
                             'NETOGLYC'         =>  [69..90] ,
                             'NETPHOS'          =>  [91..140] ,
                             'PEST'             =>  [141..152] ,
                             'PSIPRED_helices'  =>  [163..171,192..201] ,
                             'PSIPRED_sheets'   =>  [155..162,182..191] ,
                             'PSIPRED_rcoils'   =>  [153..154,172..181] ,
                             'PSORT'            =>  [202..213] ,
                             'SEQFEAT'          =>  [214..250] ,
                             'SIGP'             =>  [251..258]
                           };

# Web server version - the following bit that defines $inputs is modified so
# that the %$inputs hash actually contains just one entry, obtained from the
# input path $fasta : the value is the folder and the key is everything else
# (i.e. the filename root). Note that in this version $fasta is NOT the full
# path to the FASTA file.
my ($input_folder, $input_rootname) = ($1, $2) if ($fasta =~ m{(.*?)([^/]+)$});
$input_folder =~ s{/*$}{};
die "Bad input to FFPred.pl !\n" unless ($input_folder && $input_rootname);
my $inputs = {$input_rootname => $input_folder};
## Read input sequences into a hashref.                 # Skipped in the web server version.
#my $inputs = {};                                       # Skipped in the web server version.
#ReadFasta($fasta, $inputs)                             # Skipped in the web server version.
#    or die "\n----- error in sub ReadFasta -----\n\n"; # Skipped in the web server version.

# Parse SVM data.
my $GOterms = {};
foreach my $criteria_folder (glob "$FFPred_dir/SVMs/*") # Edited in the web server version.
{
    next unless (-d $criteria_folder);
    my $criteria = $1 if ($criteria_folder =~ m{([^/]+)$});
    
    # Read logistic regression values for this GO term.
    my $input_file = "$criteria_folder/${criteria}_ABfile";
    open(ABFILE, "<", $input_file) or die "Cannot open file $input_file !\n";
    while (defined(my $line = <ABFILE>))
    {
        next if ($line =~ /^#/);
        @{$GOterms->{$criteria}{$1}}{'Avalue', 'Bvalue'} = ($2, $3) if ($line =~ /\s*(\S+)\s+(\S+)\s+(\S+)/);
    }
    close(ABFILE) or print STDERR "Cannot close file $input_file !\n";
    
    # Read which feature groups are used by this GO term.
    $input_file = "$criteria_folder/${criteria}_feature_file";
    open(FEATURE, "<", $input_file) or die "Cannot open file $input_file !\n";
    while (defined(my $line = <FEATURE>))
    {
        next if ($line =~ /^#/);
        if ($line =~ /\s*(\S+)\s+(.+?)\s*$/)
        {
            next unless (exists($GOterms->{$criteria}{$1}));
            foreach my $feature_group (split(' ', $2))
            {
                $GOterms->{$criteria}{$1}{'features'}{$feature_group} = 1;
            }
        }
    }
    close(FEATURE) or print STDERR "Cannot close file $input_file !\n";
    
    # Read performance evaluation values for this GO term.
    $input_file = "$criteria_folder/${criteria}_long_summary";
    open(SUMM, "<", $input_file) or die "Cannot open file $input_file !\n";
    while (defined(my $line = <SUMM>))
    {
        next if ($line =~ /^#/);
        if ($line =~ /\s*(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+?)\s*$/)
        {
            next unless (exists($GOterms->{$criteria}{$1}));
            @{$GOterms->{$criteria}{$1}}{'MCC', 'sensitivity', 'specificity', 'precision', 'domain', 'description'} = ($2, $3, $4, $5, $GOdomains->{$6}, $7);
        }
    }
    close(SUMM) or print STDERR "Cannot close file $input_file !\n";
    
    foreach my $GOterm (keys %{$GOterms->{$criteria}})
    {
        $GOterms->{$criteria}{$GOterm}{'name'} = "GO:$1" if ($GOterm =~ /GO(\d+)/);
        
        # Criteria for deciding whether the GO term's SVM is reliable.
        my $safe_condition = ( ($GOterms->{$criteria}{$GOterm}{'MCC'}         >= 0.3) && 
                               ($GOterms->{$criteria}{$GOterm}{'sensitivity'} >= 0.3) && 
                               ($GOterms->{$criteria}{$GOterm}{'specificity'} >= 0.7) && 
                               ($GOterms->{$criteria}{$GOterm}{'precision'}   >= 0.3) );
        $GOterms->{$criteria}{$GOterm}{'safe'} = $safe_condition ? 1 : "";
    }
}

# Perform FFPred jobs.
foreach my $id (keys %$inputs)
{
    my ($error, $job, $pred) = (0, {}, {});
    
    $job->{'id'} = $id;
#    $job->{'seq'} = $inputs->{$id};        # Skipped in the web server version.
#    $job->{'len'} = length($job->{'seq'}); # Skipped in the web server version.
    $job->{'submitted'} = $submit_datetime;
    
#    # Find a root name for output files using this sequence's FASTA header. # Skipped in the web server version.
#    my $sanitise_name = ($job->{'id'} =~ /(\S{1,15})/) ? $1 : 'default';
#    $sanitise_name =~ s/[^\w\-]+/__/g;
#    $sanitise_name =~ s/_+$//;
#    my $output_name = my $output_rootname = $sanitise_name;
#    my $output_folder = "$dirOut/FFPred_${output_rootname}";
#    my $suffix_dir = 0;
#    while (!mkdir($output_folder))
#    {
#        die "Too many similar output folder names or folder creation error !\n" if (++$suffix_dir > 10000000);
#        $output_name = "${output_rootname}_${suffix_dir}";
#        $output_folder = "$dirOut/FFPred_${output_name}";
#    }
#    
#    $job->{'out'} = "$output_folder/$output_name";
    $job->{'out'} = "$inputs->{$id}/$id"; # Added in the web server version.
    
#    GetMD5($job);       # Skipped in the web server version.
    $job->{'md5'} = $id; # Added in the web server version.
    
#    print STDERR "\n".Timestamp()." - Started job $job->{'id'}.\n" . 
#                 "Input file $fasta\n" . 
#                 "Input sequence's md5 code $job->{'md5'}\n" . 
#                 "Output folder $output_folder\n";
    
#    if (($job->{'len'} < 15) || ($job->{'len'} > 2000))
#    {
#        print STDERR "\n\nBEWARE !\n" . 
#                     "$job->{'id'} is being processed, but its sequence length is outside FFPred's range.\n" . 
#                     "FFPred was optimised using sequences between 15aa and 2000aa - $job->{'len'} aa were submitted.\n\n";
#    }
#    elsif ($job->{'seq'} =~ /^[CTAG]+$/)
#    {
#        print STDERR "\n\nBEWARE !\n" . 
#                     "$job->{'id'} is being processed, but it looks like DNA.\n" . 
#                     "FFPred is intended for protein sequences.\n\n";
#    }
    
#    MkDir($jobs_directory, $job); # Skipped in the web server version.
    $job->{'dir'} = $inputs->{$id}; # Added in the web server version.
    
#    print STDERR "\n".Timestamp()." - Temporary folder $job->{'dir'} - Started featurama.\n\n"; # Skipped in the web server version.
#    $error += RunFeaturama($FFPred_dir, $job);                                                  # Skipped in the web server version.
#    print STDERR "\n".Timestamp()." - Temporary folder $job->{'dir'} - Finished featurama!\n";  # Skipped in the web server version.
#    
#    die "Error ($error) while running featurama - ABORTING !\n" if ($error);                    # Skipped in the web server version.
    
    Reformat($job);
#    Copy_features($job); # Skipped in the web server version.
    
    $job->{'features'} = Parse_features($job)
        or die "\n----- error in sub Parse_features -----\n\n";
    
    print STDERR "\n".Timestamp()." - Started using SVM library.";
    $error += RunSVM($FFPred_dir, $feature_groups_lists, $GOterms, $job, $pred);
    print STDERR "\n".Timestamp()." - Finished using SVM library!\n";
    
    die "Error ($error) while running SVMs - ABORTING !\n" if ($error);
    
    PrintSVMresults($GOterms, $job, $pred)
        or die "\n----- error in sub PrintSVMresults -----\n\n";
    print STDERR "\n".Timestamp()." - Finished printing SVM results!\n";
    
    # Use the following line to automatically delete all temporary files.
#    Clean_temp($job); # Skipped in the web server version.
    
    print STDERR "\n".Timestamp()." - Finished job $job->{'id'}.\n\n\n";
}

exit(0);



################################################## Subroutines after this line.



sub Timestamp
{
    my @MONTHS = ('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December');
    my ($sec, $min, $hour, $day, $month, $year) = (localtime)[0..5];
    
    return "$day $MONTHS[$month] " . ($year+1900) . ", $hour:$min:$sec";
}

sub ReadConfigPATH
{
    my ($rootdir) = @_;
    my $cfgPATH = "";
    
    open(CONFIG, "<", "$rootdir/CONFIG") or print STDERR "Cannot read CONFIG - $!\n" and return; # Edited in the web server version.
    
    while (defined(my $line = <CONFIG>))
    {
        next if ($line =~ /^#/);
        if ($line =~ /\bPATH\s+(\S+)/)
        {
            $cfgPATH = $1;
            $cfgPATH = "$rootdir/$cfgPATH" unless (substr($cfgPATH, 0, 1) eq '/');
        }
    }
    
    close(CONFIG) or print STDERR "Cannot close CONFIG - $!\n";
    
    print STDERR "'PATH' for temporary jobs not found in CONFIG file!\n" unless ($cfgPATH);
    return $cfgPATH;
}

sub ReadFasta
{
    my ($fasta_filepath, $inputs_hash) = @_;
    
    local $/;
    open(FASTA, "<", $fasta_filepath) or print STDERR "Cannot open file $fasta_filepath - $!\n" and return;
    my $content = <FASTA>;
    close(FASTA) or print STDERR "Cannot close file $fasta_filepath - $!\n";
    
    print STDERR "Invalid input FASTA file !\n" and return unless ($content =~ />/);
    $content =~ s/\A.*?^>/>/ms; # Avoid any initial lines without a valid FASTA header.
    my $key;
    
    foreach my $line (split("\n", $content))
    {
        if ($line =~ /^>/)
        {
            $line =~ s/^>\s*//;
            $line =~ s/\s*$//;
            
            if ($line ne "")
            {
                $key = $line;
            }
            else
            {
                print STDERR "Invalid FASTA header in input file !\n" and return;
            }
        }
        elsif ($line !~ /^\s*$/)
        {
            $line =~ s/\s//g;
            $line = uc $line;
            
            if ($line =~ /^[A-Z]+$/)
            {
                $inputs_hash->{$key} .= $line;
            }
            else
            {
                print STDERR "Invalid line in input file !\n" and return;
            }
        }
    }
    
    return 1;
}

sub GetMD5
{
    my ($job_hash) = @_;
    
    use Digest::MD5 qw(md5_hex);
    
    my %hextobinhash =   ('x0'=>'0000', 'x1'=>'0001', 'x2'=>'0010',
                          'x3'=>'0011', 'x4'=>'0100', 'x5'=>'0101',
                          'x6'=>'0110', 'x7'=>'0111', 'x8'=>'1000',
                          'x9'=>'1001', 'xa'=>'1010', 'xb'=>'1011',
                          'xc'=>'1100', 'xd'=>'1101', 'xe'=>'1110',
                          'xf'=>'1111');
    
    my %bintohex32hash = ('x00000'=>'0', 'x00001'=>'1', 'x00010'=>'2',
                          'x00011'=>'3', 'x00100'=>'4', 'x00101'=>'5',
                          'x00110'=>'6', 'x00111'=>'7', 'x01000'=>'8',
                          'x01001'=>'9', 'x01010'=>'a', 'x01011'=>'b',
                          'x01100'=>'c', 'x01101'=>'d', 'x01110'=>'e',
                          'x01111'=>'f', 'x10000'=>'g', 'x10001'=>'h',
                          'x10010'=>'i', 'x10011'=>'j', 'x10100'=>'k',
                          'x10101'=>'l', 'x10110'=>'m', 'x10111'=>'n',
                          'x11000'=>'o', 'x11001'=>'p', 'x11010'=>'q',
                          'x11011'=>'r', 'x11100'=>'s', 'x11101'=>'t',
                          'x11110'=>'u', 'x11111'=>'v');
    
    my $enc = md5_hex($job_hash->{'seq'});
    my $md5cut = substr($enc, -20, 20);
    
    my ($char, $binary, $hex_md5, $bin);
    
    for (my $x = 0; $x < 20; $x++)
    {
        $char = "x" . substr($md5cut, $x, 1);
        $binary .= $hextobinhash{$char};
    }
    
    $hex_md5 = "";
    
    for (my $y = 0; $y < 80; $y += 5)
    {
        $bin = "x" . substr($binary, $y, 5);
        $hex_md5 .= $bintohex32hash{$bin};
    }
    
    $job_hash->{'md5'} = $hex_md5;
}

sub MkDir
{
    my ($jobs_dir, $job_hash)= @_;
    
    $job_hash->{'dir'} = my $job_root = "$jobs_dir/" . substr($job_hash->{'md5'}, 0, 5);
    
    # This avoids race conditions when many jobs are run at the same time.
    my $suffix_tempdir = 0;
    while (!mkdir($job_hash->{'dir'}))
    {
        die "Too many similar temporary folder names or folder creation error !\n" if (++$suffix_tempdir > 10000000);
        $job_hash->{'dir'} = "${job_root}_${suffix_tempdir}";
    }
    
    open(OUT, ">", "$job_hash->{'dir'}/$job_hash->{'md5'}.fsa");
    print OUT ">$job_hash->{'md5'}\n$job_hash->{'seq'}\n";
    close(OUT);
}

sub RunFeaturama
{
    my ($rootdir, $job_hash) = @_;
    
    my $featurama = "$rootdir/featurama/features.pl " . # Edited in the web server version.
                    "-d $job_hash->{'dir'} " . 
                    "-i $job_hash->{'dir'}/$job_hash->{'md5'}.fsa " . 
                    "-o $job_hash->{'dir'}/$job_hash->{'md5'}.results " . 
                    "-fconfig $job_hash->{'dir'}/$job_hash->{'md5'}.featcfg";
    $featurama .= " -orphan" if ($job_hash->{'len'} <= 20);
    
    return system($featurama);
}

sub Reformat
{
    my ($job_hash) = @_;
    my $results = "$job_hash->{'dir'}/$job_hash->{'md5'}.results";
    my $features = "$job_hash->{'dir'}/$job_hash->{'md5'}.features";
    
    open(RESULTS, "<", $results) or die "Cannot open file $results !\n";
    
    my $header = <RESULTS>; #           First line of '.results' file.
    my $raw_output = <RESULTS>; #      Second line of '.results' file.
    my $scaled_output = <RESULTS>; #    Third line of '.results' file.
    
    close(RESULTS) or print STDERR "Cannot close file $results - $!\n";
    
    my @scaled_values = split(/\s+/, $scaled_output);
    
    open(FEATURES, ">", $features) or die "Cannot open file $features !\n";
    
    my $sequenceID = $scaled_values[0];
    $sequenceID =~ s/\w+_//;
    print FEATURES "$sequenceID\t";
    my $feature_number = 0;
    
    for (my $index = 0; $index < @scaled_values; $index++)
    {
        my $value = $scaled_values[$index];
        next if ($value =~ /$sequenceID/);
        
        $feature_number++;
        print FEATURES "${feature_number}:$value\t";
    }
    
    print FEATURES "\n";
    close(FEATURES) or print STDERR "Cannot close file $features - $!\n";
}

sub Copy_features
{
    my ($job_hash) = @_;
    
    copy("$job_hash->{'dir'}/$job_hash->{'md5'}.features", "$job_hash->{'out'}.features");
    copy("$job_hash->{'dir'}/$job_hash->{'md5'}.results", "$job_hash->{'out'}.results");
}

sub Parse_features
{
    my ($job_hash) = @_;
    my $filepath_features = "$job_hash->{'out'}.features";
    
    open(FEATURES, "<", $filepath_features) or print STDERR "Cannot open file $filepath_features - $!\n" and return;
    my $features = <FEATURES>;
    close(FEATURES) or print STDERR "Cannot close file $filepath_features - $!\n";
    
    $features =~ s/^.+?\s+//;
    chomp $features;
    my @list_in = split /\t/, $features;
    my @list_out = ();
    foreach my $item (@list_in)
    {
        $item =~ /\d+:(.*)/;
        push(@list_out, $item) unless ($1 == 0);
    }
    
    if (@list_out > 0)
    {
        return ' '.join(' ', @list_out);
    }
    else
    {
        print STDERR "\n\nNo features predicted for $job_hash->{'id'} - Aborting!\n" and return;
    }
}

sub RunSVM
{
    my ($rootdir, $feature_groups, $GOterms_hash, $job_hash, $pred_hash) = @_;
    my $svm_classify = "/webdata/binaries/current/svm_light-6.02/svm_classify"; # Edited in the web server version.
    my $error_SVMlight = 0;
    
    foreach my $criteria (keys %$GOterms_hash)
    {
        $pred_hash->{$criteria} = {};
        
        foreach my $GOterm (sort keys %{$GOterms_hash->{$criteria}})
        {
            my $GOhash = $GOterms_hash->{$criteria}{$GOterm};
            
            my $filepath_SVMinput = "$job_hash->{'dir'}/$job_hash->{'md5'}.${criteria}_${GOterm}_test";
            my $filepath_SVMmodel = "$rootdir/SVMs/$criteria/${criteria}_models/${criteria}_${GOterm}.model"; # Edited in the web server version.
            my $filepath_SVMpred = "$job_hash->{'dir'}/$job_hash->{'md5'}.${criteria}_${GOterm}_prediction";
            
            # Create a test file including only features used by this GO term.
            my @bad_features_list = ();
            foreach my $feature_group (keys %$feature_groups)
            {
                push @bad_features_list, @{$feature_groups->{$feature_group}} unless (exists($GOhash->{'features'}{$feature_group}));
            }
            
            my $bad_features = "";
            $bad_features = '(?:' . join('|', @bad_features_list) . ')' if (scalar @bad_features_list);
            
            my $used_features = $job_hash->{'features'};
            $used_features =~ s/\s$bad_features:\S+//g;
            
            open(INPUT, ">", $filepath_SVMinput) or die "Cannot open file $filepath_SVMinput - $!\n";
            print INPUT "0$used_features\n";
            close(INPUT) or print STDERR "Cannot close file $filepath_SVMinput - $!\n";
            
            # Run svm_classify to obtain the prediction.
            my $classify_command = "$svm_classify $filepath_SVMinput $filepath_SVMmodel $filepath_SVMpred";
            #print STDERR "$classify_command\n";
            $error_SVMlight += system($classify_command);
            
            # Use the prediction to obtain a posterior probability.
            open(PRED, "<", $filepath_SVMpred) or die "Cannot open file $filepath_SVMpred - $!\n";
            my $prediction = <PRED>;
            close(PRED) or print STDERR "Cannot close file $filepath_SVMpred - $!\n";
            unlink $filepath_SVMinput, $filepath_SVMpred; # Use this to selectively delete only temporary files used by SVM-Light.
            
            $prediction =~ s/\s//g;
            my $posterior = sprintf("%.3f", 1 / (1 + exp($GOhash->{'Avalue'}*$prediction + $GOhash->{'Bvalue'})));
            
            # Place the GO term in the correct set if it is predicted.
            # Also, separately save probabilities for all GO terms, including
            # those that do not satisfy the prediction condition, to be printed
            # out in a separate file.
            my $prediction_condition = ($posterior >= 0.5);
            
            if ($GOhash->{'safe'})
            {
                $pred_hash->{$criteria}{'safe'}{$GOterm} = $posterior if ($prediction_condition);
                $pred_hash->{$criteria}{'all_safe'}{$GOterm} = $posterior;
            }
            else
            {
                $pred_hash->{$criteria}{'unsafe'}{$GOterm} = $posterior if ($prediction_condition);
                $pred_hash->{$criteria}{'all_unsafe'}{$GOterm} = $posterior;
            }
        }
    }
    
    return $error_SVMlight;
}

sub PrintSVMresults
{
    my ($GOterms_hash, $job_hash, $pred_hash) = @_;
    
    foreach my $criteria (keys %$pred_hash)
    {
        my $filepath_formatted = "$job_hash->{'out'}.${criteria}_formatted";
        my $filepath_raw = "$job_hash->{'out'}.${criteria}_raw";
        my $filepath_all = "$job_hash->{'out'}.${criteria}_all";
        my @GOterms_safe = sort {$pred_hash->{$criteria}{'safe'}{$b} <=> $pred_hash->{$criteria}{'safe'}{$a}} keys %{$pred_hash->{$criteria}{'safe'}};
        my @GOterms_unsafe = sort {$pred_hash->{$criteria}{'unsafe'}{$b} <=> $pred_hash->{$criteria}{'unsafe'}{$a}} keys %{$pred_hash->{$criteria}{'unsafe'}};
        
        open(FORMATTED, ">", $filepath_formatted) or print STDERR "Cannot open file $filepath_formatted - $!\n" and return;
        open(RAW, ">", $filepath_raw) or print STDERR "Cannot open file $filepath_raw - $!\n" and return;
        
        print FORMATTED "--------------------------------------------------------------------------------\n" . 
                        "          \\\\ FFPred // version $version, subversion: $subversion training\n" . 
                        "                    Results for \"$job_hash->{'id'}\" - $criteria criteria\n" . 
                        "                         Job md5: $job_hash->{'md5'}\n" . 
                        "                    Submitted on: $job_hash->{'submitted'}\n" . 
                        "--------------------------------------------------------------------------------\n\n" . 
                        "                               Column content:\n\n" . 
                        "          Score       - Posterior probability for the annotation\n" . 
                        "          GO term     - GO term code\n" . 
                        "          RL          - Reliability Level for that GO term (High or Low)\n" . 
                        "          Domain      - Ontology domain for that GO term (CC, MF, BP)\n" . 
                        "          Description - Full GO term name\n\n" . 
                        "------------------------------- GO TERM RESULTS --------------------------------\n";
        
        print RAW "#   FFPred version $version, subversion: $subversion training\n" . 
                  "#   Results for \"$job_hash->{'id'}\" - $criteria criteria\n"; 
        
        if ((scalar(@GOterms_safe) > 0) || (scalar(@GOterms_unsafe) > 0))
        {
            print FORMATTED "Score\tGO term\t\tRL\tDomain\tDescription\n" . 
                            "--------------------------------------------------------------------------------\n";
            print RAW "#Score\tGO term\t\tRL\tDomain\tDescription\n";
        }
        else
        {
            print FORMATTED "\n-----------   NO GO TERMS CONFIDENTLY PREDICTED FOR THIS SEQUENCE !   ----------\n\n";
        }
        
        if (scalar(@GOterms_safe) > 0)
        {
            foreach my $GOterm (@GOterms_safe)
            {
                my $GOterm_result = "$pred_hash->{$criteria}{'safe'}{$GOterm}\t$GOterms_hash->{$criteria}{$GOterm}{'name'}\tH\t$GOterms_hash->{$criteria}{$GOterm}{'domain'}\t$GOterms_hash->{$criteria}{$GOterm}{'description'}\n";
                print RAW $GOterm_result;
                print FORMATTED $GOterm_result;
            }
        }
        else
        {
            print FORMATTED "\n-------   NO 'SAFE' GO TERMS CONFIDENTLY PREDICTED FOR THIS SEQUENCE !   -------\n\n";
        }
        
        if (scalar(@GOterms_unsafe) > 0)
        {
            print FORMATTED "\n----  BEWARE: the following terms are always predicted as mere speculation  ----\n\n";
            foreach my $GOterm (@GOterms_unsafe)
            {
                my $GOterm_result = "$pred_hash->{$criteria}{'unsafe'}{$GOterm}\t$GOterms_hash->{$criteria}{$GOterm}{'name'}\tL\t$GOterms_hash->{$criteria}{$GOterm}{'domain'}\t$GOterms_hash->{$criteria}{$GOterm}{'description'}\n";
                print RAW $GOterm_result;
                print FORMATTED $GOterm_result;
            }
        }
        else
        {
            print FORMATTED "\n------   NO 'UNSAFE' GO TERMS CONFIDENTLY PREDICTED FOR THIS SEQUENCE !   ------\n\n";
        }
        
        print FORMATTED "--------------------------------------------------------------------------------\n";
        
        close(RAW) or print STDERR "Cannot close file $filepath_raw - $!\n";
        close(FORMATTED) or print STDERR "Cannot close file $filepath_formatted - $!\n";
        
        # Separately print a file with probabilities for all GO terms.
        
        open(ALL, ">", $filepath_all) or print STDERR "Cannot open file $filepath_all - $!\n" and return;
        
        print ALL "#   FFPred version $version, subversion: $subversion training\n" . 
                  "#   Results for \"$job_hash->{'id'}\" - $criteria criteria\n" . 
                  "#   Submitted on: $job_hash->{'submitted'}\n" . 
                  "#\n" . 
                  "#   Column content:\n" . 
                  "#   Score       - posterior probability for the annotation\n" . 
                  "#   GO term     - GO term name\n" . 
                  "#   RL          - Reliability Level for that GO term (High or Low)\n" . 
                  "#   Domain      - Ontology domain for that GO term (CC, MF, BP)\n" . 
                  "#   Description - Full GO term name\n" . 
                  "#\n" . 
                  "#Score\tGO term\t\tRL\tDomain\tDescription\n";
        
        foreach my $GOterm (sort {$pred_hash->{$criteria}{'all_safe'}{$b} <=> $pred_hash->{$criteria}{'all_safe'}{$a}} keys %{$pred_hash->{$criteria}{'all_safe'}})
        {
            print ALL "$pred_hash->{$criteria}{'all_safe'}{$GOterm}\t$GOterms_hash->{$criteria}{$GOterm}{'name'}\tH\t$GOterms_hash->{$criteria}{$GOterm}{'domain'}\t$GOterms_hash->{$criteria}{$GOterm}{'description'}\n";
        }
        
        foreach my $GOterm (sort {$pred_hash->{$criteria}{'all_unsafe'}{$b} <=> $pred_hash->{$criteria}{'all_unsafe'}{$a}} keys %{$pred_hash->{$criteria}{'all_unsafe'}})
        {
            print ALL "$pred_hash->{$criteria}{'all_unsafe'}{$GOterm}\t$GOterms_hash->{$criteria}{$GOterm}{'name'}\tL\t$GOterms_hash->{$criteria}{$GOterm}{'domain'}\t$GOterms_hash->{$criteria}{$GOterm}{'description'}\n";
        }
        
        close(ALL) or print STDERR "Cannot close file $filepath_all - $!\n";
    }
    
    return 1;
}

sub Clean_temp
{
    my ($job_hash) = @_;
    
    unlink(glob("$job_hash->{'dir'}/$job_hash->{'md5'}" . "[._]*"));
    rmdir $job_hash->{'dir'} or print STDERR "Could not delete folder $job_hash->{'dir'} - $!\n";
}
