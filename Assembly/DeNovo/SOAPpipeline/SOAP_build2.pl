use Cwd;
use File::Basename;

my $rootdir = getcwd;

#### Get script arguments
$CPU = $ARGV[0];

#### Open projects file
system("mkdir data.out tmp tmp/fastq tmp/assemblies");

my @list = `ls reads/*.fq\_1`;

#### loop over all datasets
foreach my $genome (@list){

  # prepare infiles
  chomp($genome);
  my $bsn = basename($genome);
  $bsn =~ s/\.fq\_1//;
  my $r2 = $genome;
  $r2 =~ s/\.f\_1/\.fq\_2/;

  ### Get CPU load
  open PIPE, "uptime |";
  my $line = <PIPE>;
  close PIPE;
  $line =~ s/\s//g;
  my @lineArr =  split /:/, $line;
  my $times = $lineArr[@lineArr-1];
  my @timeArr = split /,/, $times;
  my $load = $timeArr[0] + 1;
  print "Current CPU load is $load \n\n\n";

  if($load < $CPU) {
    print "OK --- start $bsn\n";

    $command = "cp $genome tmp/fastq/R1.reads.fq";
    print "Running: $command\n\n";
    system("$command");   
   
    $command = "cp $r2 tmp/fastq/R2.reads.fq";
    print "Running: $command\n\n";
    system("$command");  

    # run SOAP
    $command = "./bin/SOAPdenovo-63mer all -K 63 -p 10 -M 2 -D 5 -F -s params/SOAPpe.config -o tmp/assemblies/$bsn";
#     $command = "./bin/SOAPdenovo-127mer all -K 61 -p 10 -M 3 -D 5 -s params/SOAP.config -o tmp/assemblies/$genome";
    print "Running: $command\n\n";
    system("$command");

    # move outputs
    $command = "cp tmp/assemblies/$bsn.scafS* tmp/assemblies/$bsn.contig data.out/";
    print "Running: $command\n\n";
    system("$command");

    # clean mess
    $command = "rm tmp/fastq/* tmp/assemblies/$bsn.*";
    print "\n\n\n##########\n##########\nRUNNING:\n$command\n##########\n";
    system("$command");

    } else {
    print "WAIT... CPU load is maximised\n\n\n";
    redo;
    }
  }

