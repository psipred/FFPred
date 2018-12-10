#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/lib";
use PsiPred;
use DisoPred;
use LowComplex;
use Memsat;
use Coils;
use NetNGlyc;
use NetOGlyc;
use NetPhos;
use PEST;
use Psort;
use SeqFeat;
use SigP;
use Cwd;
use Data::Dumper;

my $root_directory = $FindBin::Bin; # This needs to be set to the full path to the main FFPred2 folder.
my $relative_path_CONFIG = '../CONFIG';

my ($folder, $fasta, $help, $orphan, $verbose, $fconfig, $out, $web_control);

my $args = GetOptions(
                       "d=s"       => \$folder,
                       "i=s"       => \$fasta,
                       "o=s"       => \$out,
                       "fconfig=s" => \$fconfig,
                       "orphan"    => \$orphan,
                       "w"         => \$web_control,
                     );
print STDERR "No args entered... $! \n" and exit(1) unless ($args);

open(FASTA, $fasta) or die "Can't open fasta $fasta: $!\n";
my $md5 = $1 if (<FASTA> =~ />(\S+)/); #  First line in FASTA file produced by the calling script.
my $aa = <FASTA>; #                      Second line in FASTA file produced by the calling script.
chomp $aa;
close(FASTA);

if($web_control)
{
  ## $fasta =~ /(.{8}-.{4}-.{4}-.{4}-.{12})\.fsa/;
  ## $md5 = $1;
  if($fasta_file =~ /(.{8}-.{4}-.{4}-.{4}-.{12})\.fsa/){$md5 = $1;}
  elsif($fasta_file =~ /^(.{8}-.{4}-.{4}-.{4}-.{12})\.sing/){$md5 = $1;}
}
my $id = $md5; # Weird choice inherited from legacy code.

my $cfg = readConfig($root_directory, $relative_path_CONFIG);

$cfg->{'PATH'} = $folder;
#print("PATH:".$cfg->{'PATH'}."\n");
my $PRED = init($aa, $id, $md5, $cfg);
if(!$web_control)
{
  run($PRED, $orphan);
}
parse($PRED, $out);
print_results($PRED, $aa, $id, $out, 'raw');
normalise($PRED);
print_results($PRED, $aa, $id, $out, 'norm', $fconfig);

# This loop is needed due to the legacy object-oriented design used.
foreach my $prog (keys %$PRED)
{
    # print Dumper $PRED->{$prog}->err();
    if ($PRED->{$prog}->err())
    {
        print STDERR "$id failed at $prog\n";
        exit(1);
    }
}

exit(0);



############################################################# Subroutines here.



sub readConfig
{
    my ($rootdir, $rel_path_CONFIG) = @_;
    my $cfg = {};

    print("$rootdir/$rel_path_CONFIG\n");
    open(CONFIG, "<", "$rootdir/$rel_path_CONFIG") or die "Cannot read CONFIG: $!\n";

    while (defined(my $line = <CONFIG>))
    {
        next if ($line =~ /^#/);
        if ($line =~ /(\S+)\s+(\S+)/)
        {
            my ($key, $value) = ($1, $2);
            $value = "$value" unless (substr($value, 0, 1) eq '/');
            $cfg->{$key} = $value;
        }
    }

    close(CONFIG);

    return $cfg;
}

sub init
{
    my ($aa, $id, $md5, $cfg) = @_;
    my $PROGS = {};

    $PROGS->{SIGP}     = new SigP($aa, $id, $md5, $cfg);
    if($web_control)
    {
      $PROGS->{DISOPRED} = new DisoPred($aa, $id, $md5, $cfg, 1);
      $PROGS->{MEMSAT}   = new Memsat($aa, $id, $md5, $cfg, 1);
      $PROGS->{PSIPRED}  = new PsiPred($aa, $id, $md5, $cfg, 1);
      $PROGS->{SEQFEAT}  = new SeqFeat($aa, $id, $md5, $cfg, 1);
    }
    else
    {
      $PROGS->{DISOPRED} = new DisoPred($aa, $id, $md5, $cfg, 0);
      $PROGS->{MEMSAT}   = new Memsat($aa, $id, $md5, $cfg, 0);
      $PROGS->{PSIPRED}  = new PsiPred($aa, $id, $md5, $cfg, 0);
      $PROGS->{SEQFEAT}  = new SeqFeat($aa, $id, $md5, $cfg, 0);
    }
    $PROGS->{NETNGLYC} = new NetNGlyc($aa, $id, $md5, $cfg);
    $PROGS->{NETOGLYC} = new NetOGlyc($aa, $id, $md5, $cfg);
    $PROGS->{NETPHOS}  = new NetPhos($aa, $id, $md5, $cfg);
    $PROGS->{PEST}     = new PEST($aa, $id, $md5, $cfg);
    $PROGS->{PSORT}    = new Psort($aa, $id, $md5, $cfg);
    $PROGS->{COILS}    = new Coils($aa, $id, $md5, $cfg);
    $PROGS->{LOWC}     = new LowComplex($aa, $id, $md5, $cfg);

    return $PROGS;
}

sub run
{
    my ($PRED, $orphan) = @_;
    my @blastPreds = ('PSIPRED', 'DISOPRED', 'MEMSAT');

    print STDERR "Running " . join(', ', @blastPreds) . "\n\n";

    unless ($orphan)
    {
        foreach my $blastPred (@blastPreds)
        {
            $PRED->{$blastPred}->runPsiblast();
            if ($PRED->{$blastPred}->err())
            {
                print STDERR "\n$blastPred\'s psiblast results: FATAL ERROR !\n";
                exit(1);
            }
        }

        print STDERR "The psiblast runs for Psipred, Disopred and Memsat finished successfully.\n";
    }

    foreach my $blastPred (@blastPreds)
    {
        $PRED->{$blastPred}->run();
        if ($PRED->{$blastPred}->err())
        {
            print STDERR "$blastPred\'s run failed !\n";
            exit(1);
        }
    }

    print STDERR "\n" . join(' / ', @blastPreds) . " : success\n\n";

    foreach my $prog (sort keys %$PRED)
    {
        if ($prog !~ /PSIPRED/ && $prog !~ /DISOPRED/  && $prog !~ /MEMSAT/)
        {
            print STDERR "Running $prog\n";
            $PRED->{$prog}->run();

            if ($PRED->{$prog}->err())
            {
                print STDERR "$prog failed !\n";
                exit(1);
            }
            print STDERR "$prog : success\n";
        }
    }
}

sub parse
{
    my ($PRED, $out) = @_;

    foreach my $prog (sort keys %$PRED)
    {
        print STDERR "Parsing $prog\n";
	      $PRED->{$prog}->parse();
        $PRED->{$prog}->print_keys($out, $prog);
    }

    open(OUT, ">>", $out);
    print OUT "\n";
    close(OUT);
}

sub print_results
{
    my ($PRED, $aa, $id, $out, $type, $fconf) = @_;

    if (defined($fconf))
    {
        open(CONFIG, ">", $fconf);
        print CONFIG "$id\t$aa\n";
        close(CONFIG);
    }

    foreach my $prog (sort keys %$PRED)
    {
        print STDERR "Printing $prog\n";
	      $PRED->{$prog}->print_results($out, $type);
        $PRED->{$prog}->print_config($fconf) if (defined($fconf));
    }

    open(OUT, ">>", $out);
    print OUT "\n";
    close(OUT);
}

sub normalise
{
    my ($PRED) = @_;

    foreach my $prog (keys %$PRED)
    {
	$PRED->{$prog}->normalise();
    }
}



# Old (legacy) versions of run and check subs, which were kept separate in
# order to allow use of children processes for the PSI-BLAST runs using 'fork'.
# They are here to provide a guide for possible future implementation of this.

sub run_OLD
{
    my ($PRED,$orphan) = @_;

    unless ($orphan)
    {
        print STDERR "Running PSIPRED, DISOPRED, MEMSAT\n\n";
        $PRED->{PSIPRED}->runPsiblast();
        $PRED->{DISOPRED}->runPsiblast();
        $PRED->{MEMSAT}->runPsiblast();
    }

    foreach my $prog (keys %$PRED)
    {
        if($prog !~ /PSIPRED/ && $prog !~ /DISOPRED/  && $prog !~ /MEMSAT/)
        {
            print STDERR "Running $prog\n";
            $PRED->{$prog}->run();

            if($PRED->{$prog}->err())
            {
                print STDERR "$prog failed\n";
                exit(1);
            }
            print STDERR "$prog : success\n";
        }
    }
}

sub check_OLD
{
    my ($PRED, $orphan, $fork) = @_;

    unless ($orphan)
    {
        # NOTE : Complete this. (TO DO!)
        if ($fork)
        {
            my $time = 0;

            # initial pause
            sleep(5);

            # NOTE : Improve this condition & avoid infinite loop: what if blast has finished, but with an error ? (TO DO!)
            while ( !($PRED->{PSIPRED}->checkPsiblast() && $PRED->{DISOPRED}->checkPsiblast() && $PRED->{MEMSAT}->checkPsiblast()) )
            {
                sleep(5);
                print STDERR ".....$time.....\n";
                $time +=5;
            }

            # Final sleep in case nfs slow update.
            sleep(5);
        }
        else
        {
            if ( $PRED->{PSIPRED}->err() || $PRED->{DISOPRED}->err() || $PRED->{MEMSAT}->err() )
            {
                print STDERR "\npsiblast results: FATAL ERROR !\n";
                exit(1);
            }
        }

        print STDERR "\nThe psiblast runs for Psipred, Disopred and Memsat finished successfully.\n";
    }

    $PRED->{PSIPRED}->run();
    $PRED->{DISOPRED}->run();
    #$PRED->{MEMSAT}->runMemsat3(); # Obsolete; runs memsat3.
    $PRED->{MEMSAT}->run();

    if ( $PRED->{PSIPRED}->err() || $PRED->{DISOPRED}->err() || $PRED->{MEMSAT}->err() )
    {
        print STDERR "PsiPred/DisoPred/MemsatSVM failed\n";
        exit(1);
    }

    print STDERR "PsiPred/DisoPred/MemsatSVM : success\n\n";

    return 1;
}
