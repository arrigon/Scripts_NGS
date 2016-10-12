#!/bin/bash
for i in {0..9}
  do
  wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.0$i.tar.gz
  tar zxvpf nt.0$i.tar.gz
  rm nt.0$i.tar.gz
  done

for i in {10..24}
  do
  wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.$i.tar.gz
  tar zxvpf nt.$i.tar.gz
  rm nt.$i.tar.gz
  done
