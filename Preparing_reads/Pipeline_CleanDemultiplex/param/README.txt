Important Notes about how to prepare your param files:

BARCODES:
Put Barcodes here, following tab-delimited spreadsheet format. DO NOT LEAVE BLANK LINES, will crash pipeline. No headers either

e.g.
ToB24	AATGATGC


PRIMERS:
Put all potential sequencing primer / adapter sequences here. Will be choped out from reads. Fasta format.


CONTAMINANTS:
put here only technical sequences such phiX or Ecoli contaminants. reads matching to these references will be tossed away. Fasta format
DO NOT PUT primer sequences in here!