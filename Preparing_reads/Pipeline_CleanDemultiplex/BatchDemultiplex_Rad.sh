for i in {1..12}
  do
  perl CleanFastq_Pipeline_queue3_Rad.pl /data/data/MAS-RAD-100B/OK/MAS$i param/Barcodes_MAS_RAD\_$i.txt MASRAD$i 10
  rm -rf data.out/MASRAD$i/tmp
  rm  data.out/MASRAD$i/stacks/*.fq\_2
  done
