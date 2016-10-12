package DeconSeqConfig;

use strict;

use constant DEBUG => 0;
use constant PRINT_STUFF => 1;
use constant VERSION => '0.4.1';
use constant VERSION_INFO => 'DeconSeq version '.VERSION;

use constant ALPHABET => 'ACGTN';

use constant DB_DIR => 'db/';
use constant TMP_DIR => 'tmp/';
use constant OUTPUT_DIR => 'output/';

use constant PROG_NAME => 'bwa64';  # XXX should be replaced by MAC or 64 based on your system architecture
use constant PROG_DIR => '';

use constant DBS => {conta => {name => 'Contaminants',  #database name used for dispaly
                               db => 'conta',             #database name as defined with -p for bwa index -p ...
                               parts => 1},                       #number of parts the database consists of (usually one)
                     bact => {name => 'Bacterial genomes',
                              db => 'bactDB',
                               parts => 1},
                     vir => {name => 'Viral genomes',
                             db => 'virDB',
                               parts => 1}};
use constant DB_DEFAULT => 'hsref';

#######################################################################

use base qw(Exporter);

use vars qw(@EXPORT);

@EXPORT = qw(
             DEBUG
             PRINT_STUFF
             VERSION
             VERSION_INFO
             ALPHABET
             PROG_NAME
             PROG_DIR
             DB_DIR
             TMP_DIR
             OUTPUT_DIR
             DBS
             DB_DEFAULT
             );

1;
