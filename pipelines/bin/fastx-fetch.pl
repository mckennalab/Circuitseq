#!/usr/bin/env perl
use warnings;
use strict;

use Pod::Usage; ## uses pod documentation in usage code
use Getopt::Long qw(:config auto_help pass_through);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::Compress::Gzip; ## for creating demultiplexed files

our $VERSION = "0.1";

=head1 NAME

fastx-fetch.pl -- Retrieve sequences from a fastq/fasta file

=head1 SYNOPSIS

./fastx-fetch.pl <reads.fq> [options]

=head2 Options

=over 2

=item B<-help>

Only display this help message

=item B<-idfile> I<idFileName>

Retrieve sequence IDs from I<idFileName>

=item B<-quiet>

Don't report additional information to standard error

=item B<-reverse> or B<-v>

Invert selection logic: exclude specified sequences

=item B<-unique>

Output only the first sequences with the same sequence ID (up to first space)

=item B<-minLength> I<len>

Only output sequences that are at least as long as I<len>

=item B<-maxLength> I<len>

Only output sequences that are at most as long as I<len>

=item B<-count> I<rCount>

Stop after I<rCount> reads

=item B<-nametrim> I<pattern>

Remove text matching I<pattern> from the sequence names

=item B<-demultiplex> I<bcFileName>

Demultiplex reads using barcode assignments in I<bcFileName>. The file
format should be space, semicolon [;], or comma [,] separated, with
the barcode name (which musn't include spaces, semicolons, or commas)
as the first field.

=item B<-prefix> I<pString>

Prefix demultiplexed output file names with I<pString>

=item B<-chimeric>

Also output chimeric sequences (containing multiple barcodes).

=back

=head1 DESCRIPTION

Retrieve sequences from a fastq/fasta file.

=cut

my $idFileName = "";
my $bcFileName = "";
my $quiet = 0;
my $minLength = 1;
my $maxLength = 10 ** 12; # 1 Tbp
my $count = -1;
my $invert = 0; # invert logic
my $trim = 0;
my $trimString = "";
my $unique = 0;
my $prefix = "reads";
my $outputLengths = 0;
my $chimeric = 0;

GetOptions("idfile=s" => \$idFileName, "quiet!" => \$quiet,
           "reverse|v!" => \$invert, "trim=i" => \$trim,
           "unique!" => \$unique,
           "minLength=i" => \$minLength, "maxLength=i" => \$maxLength,
           "count=i" => \$count, "nametrim=s" => \$trimString,
           "demultiplex=s" => \$bcFileName, "prefix=s" => \$prefix,
           "lengths!" => \$outputLengths, "chimeric|x!" => \$chimeric)
  or pod2usage(1);

if($chimeric && !$bcFileName){
  pod2usage("Error: chimeric read output requires demultiplexing\n");
}

if($outputLengths && !$bcFileName){
  pod2usage("Error: read length output requires demultiplexing\n");
}

if($bcFileName && $invert){
  pod2usage("Error: demultiplexing is not compatible with logic inversion\n");
}

if($trim){
  $minLength = $minLength + $trim * 2;
  $maxLength = $maxLength + $trim * 2;
}

if($trimString){
  printf(STDERR "Will remove text matching '(%s)' from the sequence IDs",
         $trimString);
  $trimString =~ s/^\|//;
  $trimString = "($trimString)";
}

my %idsToGet = ();
my %outFiles = ();
my %idsSeen = ();
my %bcCounts = ();

# unknown commands are treated as identifiers
my @files = ();
while(@ARGV){
  my $arg = shift(@ARGV);
  if(-e $arg){
    push(@files, $arg);
  } else {
    $idsToGet{$arg} = 0;
  }
}
@ARGV = @files;

# use stdin if no files supplied
if(!@ARGV){
  @ARGV = '-' unless (-t STDIN);
}

if($idFileName){
  # read sequence IDs from input file
  if(!$quiet){
    printf(STDERR "Attempting to read from input file '$idFileName'\n");
  }
  my $idFile = new IO::Uncompress::Gunzip "$idFileName" or
    pod2usage("Unable to open $idFileName\n");
  while(<$idFile>){
    chomp;
    s/^[>@]//;
    s/[\s,].*$//;
    $idsToGet{$_} = 0;
  }
  close($idFile);
}

if($bcFileName){
  # read sequence IDs from input file
  if(!$quiet){
    printf(STDERR "Attempting to read barcodes from file '$bcFileName'\n");
  }
  my $bcFile = new IO::Uncompress::Gunzip "$bcFileName" or
    pod2usage("Unable to open $bcFileName\n");
  while(<$bcFile>){
    chomp;
    s/^[>@]//;
    my @F = split(/[\s,;]+/, $_, 2);
    $idsToGet{$F[1]} .= $F[0] . ";";
    $outFiles{$F[0]}{"present"} = 1;
  }
  close($bcFile);
}

if($bcFileName){
  ## set up out file lookup dictionary, and open files for writing
  foreach my $barcode (keys(%outFiles)){
    my $fileName = "${prefix}_${barcode}.fq.gz";
    my $lengthFileName = "${prefix}_${barcode}.txt.gz";
    $lengthFileName =~ s#/reads#/lengths#;
    if(!$lengthFileName =~ /lengths/){
      $lengthFileName =~ "${prefix}_lengths_${barcode}.txt.gz";
    }
    if(-e $fileName){
      pod2usage("Error: output file '$fileName' already exists;".
                " cannot continue\n");
    }
    my $outFile = new IO::Compress::Gzip($fileName);
    my $lengthFile = new IO::Compress::Gzip($lengthFileName);
    $outFiles{$barcode}{fileName} = $fileName;
    $outFiles{$barcode}{file} = $outFile;
    if($outputLengths){
      $outFiles{$barcode}{lengthFile} = $lengthFile;
    }
  }
}

if(!$quiet && (scalar(keys(%idsToGet)) > 0)){
  printf(STDERR "Read %d identifiers\n", scalar(keys(%idsToGet)));
}

if(!$quiet && $invert){
  printf(STDERR "Excluding IDs, rather than selecting\n");
} else {
  ## Stop when all IDs have been seen
  if(($count == -1) && (keys(%idsToGet))){
    $count = scalar(keys(%idsToGet));
  }
}

sub processSeq {
  my ($seqID, $seq, $qual, $containsFQ, $barcodes) = @_;
  if($trim > 0){
    $seq = substr($seq, $trim, length($seq)-($trim * 2));
    if($qual){
      $qual = substr($qual, $trim, length($qual)-($trim * 2));
    }
  }
  if(!$barcodes){
    if($qual){
      printf("@%s\n%s\n+\n%s\n", $seqID, $seq, $qual);
    } elsif(!$containsFQ) {
      $seq =~ s/(.{100})/$1\n/g;
      $seq =~ s/\n$//;
      printf(">%s\n%s\n", $seqID, $seq);
    }
  } else {
    while($barcodes =~ s/(.*?);//){
      my $barcodeName = $1;
      my $barcodeFile = $outFiles{$barcodeName}{file};
      $bcCounts{$barcodeName}++;
      if($qual){
        printf($barcodeFile "@%s\n%s\n+\n%s\n", $seqID, $seq, $qual);
      } elsif(!$containsFQ) {
        $seq =~ s/(.{100})/$1\n/g;
        $seq =~ s/\n$//;
        printf($barcodeFile ">%s\n%s\n", $seqID, $seq);
      }
      if($outputLengths){
        my $shortID = $seqID;
        $shortID =~ s/ .*$//;
        my $barcodeLengthFile = $outFiles{$barcodeName}{lengthFile};
        printf($barcodeLengthFile "%-8d %s\n", length($seq), $shortID);
      }
    }
  }
}

my $inQual = 0; # false
my $seqID = "";
my $qualID = "";
my $seq = "";
my $qual = "";
my $barcodes = "";
my $containsFQ = 0; # false
my $duplicateCount = 0;
foreach my $file (@ARGV) {
  # This little gunzip dance makes sure the script can handle both
  # gzip-compressed and uncompressed input, regardless of whether
  # or not it is piped
  my $z = new IO::Uncompress::Gunzip($file, "transparent", 1) or
    pod2usage("gunzip failed: $GunzipError\n");
  while(<$z>){
    chomp;
    chomp;
    if (!$inQual) {
      if(/\x00/){ ## detect corrupt files, wait for next good read
        $seqID =~ s/ .*$//;
        print(STDERR "Warning: corrupt sequence found at $seqID [NUL]\n");
        $seqID = "";
        $qualID = "";
        $seq = "";
        $qual = "";
        next;
      } elsif (/^(>|@)((.+?)( .*?\s*)?)$/) {
        my $newSeqID = $2;
        my $newShortID = $3;
        my $testSeqID = $newSeqID;
        my $testShortID = $newShortID;
        if($trimString){
          $testShortID =~ s/$trimString//;
          $testSeqID =~ s/$trimString//;
        }
        if ($seqID && (length($seq) >= $minLength) &&
            (length($seq) <= $maxLength)) {
          my $barcodeCount = ($barcodes =~ tr/;//);
          if($chimeric || $barcodeCount < 2){
            processSeq($seqID, $seq, $qual, $containsFQ, $barcodes);
          }
          if (--$count == 0) {
            $seqID = "";
            last;
          }
        }
        $seq = "";
        $qual = "";
        if ((!(keys(%idsToGet)) ||
             exists($idsToGet{$testSeqID}) ||
             exists($idsToGet{$testShortID})) xor $invert) {
          $barcodes = $idsToGet{$testShortID} unless
	      !exists($idsToGet{$testShortID});
          delete $idsToGet{$testSeqID};
          delete $idsToGet{$testShortID};
          if(exists($idsSeen{$newShortID})){
            $duplicateCount++;
            if($unique){
              $seqID = "";
            } else {
              if($duplicateCount < 20){
                print(STDERR "Warning: Sequence ID '${newShortID}' has ".
                      "already been seen. Use the '-u' option to exclude.\n");
              } elsif($duplicateCount == 20){
                print(STDERR "[Other duplicate sequences were seen]\n");
              }
              $seqID = $newSeqID;
            }
          } else {
            $seqID = $newSeqID;
          }
          $idsSeen{$newShortID} = 1;
        } else {
          $seqID = "";
        }
      } elsif (/^\+(.*)$/) {
        $inQual = 1;     # true
        $containsFQ = 1; # true
        $qualID = $1;
      } else {
        if(/@/){
          $seqID =~ s/ .*$//;
          print(STDERR "Warning: corrupt sequence found at $seqID ".
                "[header in sequence string]\n");
          $seqID = "";
          $qualID = "";
          $seq = "";
          $qual = "";
          next;
        } else {
          $seq .= $_;
        }
      }
    } else {
      if(/\x00/){ ## detect corrupt files, wait for next good read
        $seqID =~ s/ .*$//;
        print(STDERR "Warning: corrupt sequence found at $seqID [NUL]\n");
        $seqID = "";
        $qualID = "";
        $seq = "";
        $qual = "";
        $inQual = 0; # false
      } else {
        $qual .= $_;
        if (length($qual) > (length($seq) + 2)) {
          $seqID =~ s/ .*$//;
          print(STDERR "Warning: corrupt sequence found at $seqID ".
                "[quality string too long]\n");
          $seqID = "";
          $qualID = "";
          $seq = "";
          $qual = "";
          $inQual = 0; # false
        } elsif (length($qual) >= length($seq)) {
            $inQual = 0;            # false
        }
      }
    }
  }
}

if ($seqID && (length($seq) >= $minLength) &&
    (length($seq) <= $maxLength)) {
  my $barcodeCount = ($barcodes =~ tr/;//);
  if($chimeric || $barcodeCount < 2){
    processSeq($seqID, $seq, $qual, $containsFQ, $barcodes);
  }
}

if(!$quiet){
  if(!$unique && ($duplicateCount > 1)){
    print(STDERR "${duplicateCount} duplicate IDs were seen and not excluded. ".
          "Use the '-u' option to exclude duplicate sequences.\n");
  } elsif($unique) {
    print(STDERR "Number of duplicates excluded: $duplicateCount\n");
  }
}

if($bcFileName){
  ## close up files
  foreach my $barcode (sort(keys(%outFiles))){
    my $fileName = "${prefix}_${barcode}.fq.gz";
    close($outFiles{$barcode}{file});
    if($outputLengths){
      close($outFiles{$barcode}{lengthFile});
    }
    if(exists($bcCounts{$barcode})){
      printf("%9d %s\n", $bcCounts{$barcode}, $barcode);
    }
  }
}

