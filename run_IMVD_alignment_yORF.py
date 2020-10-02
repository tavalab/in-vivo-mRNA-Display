##      Panos Oikonomou, Columbia University
##      Saeed Tavazoie Lab
##
## 	(c) 2020 The Trustees of Columbia University in the City of New York. 
##	This work may be reproduced and distributed for academic 
##	non-commercial purposes only without further authorization, 
##	but rightsholder otherwise reserves all rights.
##
##	python2.7
##	requirements: cutadapt v1.8.3, bowtie2 v2.2.6, samtools v1.2, bedtools v2.17.0
##	

import subprocess
import shlex
import tempfile
import os
import re
import glob 
import string
import numpy as np 
import pandas as pd;
import pysam;

from multiprocessing import Pool

from Bio import SeqIO
import csv

DO_TRIM_ADAPTERS	= True;
DO_ALIGN_SEQUENCES	= True;
DO_PROCESS_SAMFILES	= True;
DO_MERGE_FILES		= False;
DO_NEW_MERGE_FILES	= True;

READDIR='data/'
SAVEDIR='data/out/'
TEMPDIR='data/out/intermediate/'
TRIMDIR='data/out/fastq_trimmed/'
CNTDIR ='data/out/counts/'

# adapter sequences to clip
PHRED_BASE=33
ADAP_SEQ_TAIL		= "GACGCTACGTCTCGTCGCGTGCGTAGT"
ADAP_SEQ_HEAD_5p	= "GTACAAAAAAGCAGGCTAC"
ADAP_SEQ_HEAD_3p	= "CACTTTGTACAAGAAAGC"

BARCODE_SEQ_FILE  	= "data/barcodes/myBarcodes.fa"
BARCODE_SEQ_FILE1 	= "data/barcodes/myBarcodes.anchor.fa"

LIBRARY_INDEX_ALL	= "data/library/yLinker-yORF-Library-5p-3p-index";
LIBRARY_INDEX_3P	= "data/library/yLinker-yORF-Library-3p-index";
LIBRARY_INDEX_5P	= "data/library/yLinker-yORF-Library-5p-index";

SAMDIR_ALL 		= TRIMDIR + '5p-3p' + '/'
SAMDIR_3P 		= TRIMDIR  + '3p' + '/'
SAMDIR_5P 		= TRIMDIR  + '5p' + '/'


pos_data 		= 'data/library/yLinker-yORF-Library-5p.positions'
dtcPos_5P = pd.read_csv(pos_data, header=0, index_col=0, delimiter="\t")
pos_data 		= 'data/library/yLinker-yORF-Library-3p.positions'
dtcPos_3P = pd.read_csv(pos_data, header=0, index_col=0, delimiter="\t")

LIBRARY_INDEX 		= LIBRARY_INDEX_ALL
SAMDIR 			= SAMDIR_ALL
DO_END 			= 'ALL';

MIN_LENGTH	 	= 20;

# target sequences
COMMON_PREFIX		= "R1_index_"
COMMON_SUFFIX_R1	= ""
COMMON_SUFFIX_R2 	= ""
COMMON_SUFFIX 		= ".fastq.gz" 


## for processing more than one index eg
## in_prefixes = ["ind2_CGATGTAT", "ind6_GCCAATAT", "ind8_CTTGTAAT"]; 
## change the barcode_map.tsv file accordingly!
## all internal indexes from myBarcodes.fa must be listed in barcode_map.tsv even if they were not used!
## there should be a single digit following ind !
in_prefixes 		= ["ind2_CGATGTAT"]; 
 
finalTag = '';


BARCODE_MAP	 	= "data/barcodes/barcode_map.tsv"

out_prefixes = [];

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)
if not os.path.exists(TEMPDIR):
    os.makedirs(TEMPDIR)
if not os.path.exists(TRIMDIR):
    os.makedirs(TRIMDIR)
if not os.path.exists(SAMDIR_ALL):
    os.makedirs(SAMDIR_ALL)
if not os.path.exists(SAMDIR_5P):
    os.makedirs(SAMDIR_5P)
if not os.path.exists(SAMDIR_3P):
    os.makedirs(SAMDIR_3P)
if not os.path.exists(SAMDIR):
    os.makedirs(SAMDIR)
if not os.path.exists(CNTDIR):
    os.makedirs(CNTDIR)

print;

SAMPLE_MAP = {};
PRIME_MAP = {};
SampleToIllumIndx = {}
SampleTo5pIndx = {}
SampleTo3pIndx = {}
SampleToAllIndx = {}
FileToPrimeEnd = {}

with open(BARCODE_MAP,'r') as tsvin:
    tsvin = csv.reader(tsvin, delimiter='\t')
    for row in tsvin:
     illumIndx = row[0]
     myIndx = row[1]
     primeEnd =  row[2]
     popLabel = row[3]

     print illumIndx + "\t" + myIndx + "\t" + primeEnd  + "\t" + popLabel;
     SAMPLE_MAP[(illumIndx, myIndx)] = popLabel;
     PRIME_MAP[(illumIndx, myIndx)] = primeEnd;

     SampleToIllumIndx[SAMPLE_MAP[(illumIndx, myIndx)]] = illumIndx;
     if (primeEnd == '5p'):
      SampleTo5pIndx[SAMPLE_MAP[(illumIndx, myIndx)]] = myIndx;
      if popLabel in SampleToAllIndx:
       SampleToAllIndx[popLabel] = myIndx + "." +SampleToAllIndx[popLabel];
      else:
       SampleToAllIndx[popLabel] = myIndx;   
     if (primeEnd == '3p'):
      SampleTo3pIndx[SAMPLE_MAP[(illumIndx, myIndx)]] = myIndx;
      if popLabel in SampleToAllIndx:
       SampleToAllIndx[popLabel] = myIndx + "." + SampleToAllIndx[popLabel];
      else:
       SampleToAllIndx[popLabel] = myIndx;   

print;

fasta_sequences = SeqIO.parse(open(BARCODE_SEQ_FILE),'fasta');
myIndex_dict={}
for fasta in fasta_sequences:
   name, sequence = fasta.id, fasta.seq.tostring()
   print name + "\t" + sequence; 
   myIndex_dict[name] = sequence;

print;

###################################################################################################################################################################################
###################################################################################################################################################################################

def get_prefix(reNAME, infilename):
  # return the prfix of the read information in a given filename
  mymatch = reNAME.match(infilename)
  if mymatch is None:
    return None
  else:
    return mymatch.groups()[0]

###################################################################################################################################################################################

def return_input_file(inprefix):
  return os.path.join(READDIR, COMMON_PREFIX + inprefix + COMMON_SUFFIX_R1 + COMMON_SUFFIX);

def return_temp_file(inprefix):
  return os.path.join(TEMPDIR, COMMON_PREFIX + inprefix + COMMON_SUFFIX_R1)

def return_copy_file(sample_name, prime_str, illumIdx, myIdx):
  return os.path.join(TRIMDIR, sample_name + "_orfEnd_" + prime_str + "_illumIdx_" + illumIdx + "_myIdx_" + myIdx + ".fully_trimmed" + COMMON_SUFFIX);

def return_count_file(sample_name, prime_str, illumIdx, myIdx):
  return os.path.join(CNTDIR, sample_name + "_orfEnd_" + prime_str + "_illumIdx_" + illumIdx + "_myIdx_" + myIdx + ".cnt");

def return_range_file(sample_name, prime_str, illumIdx, myIdx):
  return os.path.join(CNTDIR, sample_name + "_orfEnd_" + prime_str + "_illumIdx_" + illumIdx + "_myIdx_" + myIdx + ".range.cnt");


###################################################################################################################################################################################

def preprocess_gz_file(inprefix):
  # do some initial preprocessing of a gz file, including trimming and quality score filtering
  # give the sample-specific identifier and output prefix that should be used

  infile_fwd 	= return_input_file(inprefix);
  tmpfile_fwd 	= return_temp_file(inprefix);


  # initially:  clip the 3' universal adapter sequence
  cutadapt_cmd_0 = "cutadapt --quality-base=%i -q 20 -a %s -n 1 -o %s.cutadap_3p.fastq.gz %s > %s.cutadap_3p.log 2> %s.cutadap_3p.err" % (PHRED_BASE, ADAP_SEQ_TAIL, tmpfile_fwd, infile_fwd, tmpfile_fwd, tmpfile_fwd);

  # first: clip the 5' custom index sequence and demultiplex
  cutadapt_cmd_1 = "cutadapt --quality-base=%i -q 20 -g file:%s -O 3 --no-indels -e 0.1 --no-trim --untrimmed-o %s.myBarcodeUntrimmed.fastq.gz -o %s.trimmed-{name}.myBarcode.fastq.gz %s.cutadap_3p.fastq.gz  > %s.cutadap_myI.log 2> %s.cutadap_myI.err" % (PHRED_BASE, BARCODE_SEQ_FILE1, tmpfile_fwd, tmpfile_fwd, tmpfile_fwd, tmpfile_fwd, tmpfile_fwd);

  # second:  clip the 5' universal adapter sequence -- if the PCR came from the 5' end
  cutadapt_cmd_2 = "ls %s.trimmed-*.myBarcode.fastq.gz | while read line; do (nline=${line/myBarcode.fastq.gz/myBarcode.cut_5p.5p_PCR.fastq.gz}; uline=${line/myBarcode.fastq.gz/myBarcode.untrimmed_5p.5p_PCR.fastq.gz}; cutadapt --quality-base=%i -m %i -g %s -o $nline --untrimmed-o $uline $line); done > %s.cutadap_5p.5p_PCR.log 2> %s.cutadap_5p.5p_PCR.err" % (tmpfile_fwd, PHRED_BASE, MIN_LENGTH, ADAP_SEQ_HEAD_5p, tmpfile_fwd, tmpfile_fwd);

  # third:   clip the 5' universal adapter sequence -- if the PCR came from the 3' end
  cutadapt_cmd_3 = "ls %s.trimmed-*.myBarcode.fastq.gz | while read line; do (nline=${line/myBarcode.fastq.gz/myBarcode.cut_5p.3p_PCR.fastq.gz}; uline=${line/myBarcode.fastq.gz/myBarcode.untrimmed_5p.3p_PCR.fastq.gz}; cutadapt --quality-base=%i -m %i -g %s -o $nline --untrimmed-o $uline $line); done > %s.cutadap_5p.3p_PCR.log 2> %s.cutadap_5p.3p_PCR.err" % (tmpfile_fwd, PHRED_BASE, MIN_LENGTH, ADAP_SEQ_HEAD_3p, tmpfile_fwd, tmpfile_fwd);

  # here we do the forward and reverse reads separately; for the forward read, *only* keep reads with the correct adapter

  print "running... " +  cutadapt_cmd_0;
  print;
  subprocess.call(cutadapt_cmd_0,shell=True, executable="/bin/bash")

  print "running... " +  cutadapt_cmd_1;
  print;
  subprocess.call(cutadapt_cmd_1,shell=True, executable="/bin/bash")

  print "running... " +  cutadapt_cmd_2;
  print;
  subprocess.call(cutadapt_cmd_2,shell=True, executable="/bin/bash")

  print "running... " +  cutadapt_cmd_3;
  print;
  subprocess.call(cutadapt_cmd_3,shell=True, executable="/bin/bash")

  demultiplexed_files = glob.glob("%s.trimmed-*.myBarcode.cut_5p.*_PCR.fastq.gz" % tmpfile_fwd)
  print "%s.trimmed-*.myBarcode.cut_5p.fastq.gz" % tmpfile_fwd;
  print;

  fname_re_myIndx    = re.compile(r"%s.trimmed-(.*).myBarcode.cut_5p.(.*)_PCR.fastq.gz" % tmpfile_fwd)  

  tmpBASE = os.path.join(TEMPDIR, COMMON_PREFIX)
  fname_re_illumIndx = re.compile(r"%sind\d_(.*)%s.trimmed-(.*).myBarcode.cut_5p.(.*)_PCR.fastq.gz" % (tmpBASE, COMMON_SUFFIX_R1) )  
  # when ABOVE Index 9 (two digit indexes) use below...
  # fname_re_illumIndx = re.compile(r"%sind\d\d_(.*)%s.trimmed-(.*).myBarcode.cut_5p.(.*)_PCR.fastq.gz" % (tmpBASE, COMMON_SUFFIX_R1) )  

  print "%sind\d_(.*)%s.trimmed-(.*).myBarcode.cut_5p.(.*)_PCR.fastq.gz" % (tmpBASE, COMMON_SUFFIX_R1);
  print;
  print;

  fname_re_5p3p = re.compile(r"%s.trimmed-myIndex\d\d\d.myBarcode.cut_5p.(.*)_PCR.fastq.gz" % tmpfile_fwd )  

  prefixes = set()
  samples = set()
  for filename in demultiplexed_files: 
   print filename;

   myIndx    = myIndex_dict[("%s"   % get_prefix(fname_re_myIndx, filename))]
   print myIndx;
   illumIndx = ("%s"   % get_prefix(fname_re_illumIndx, filename))
   print illumIndx;

   copyfile = return_copy_file(SAMPLE_MAP[(illumIndx, myIndx)], PRIME_MAP[(illumIndx, myIndx)], illumIndx, myIndx);
   print filename+ "\tx" + myIndx + "x\tx" + illumIndx + "x\t That's that\t" + copyfile;   
   samples.add(SAMPLE_MAP[(illumIndx, myIndx)]);
 
   mv_cmd = "cp %s %s" %(filename, copyfile);
   if ((('5p_PCR' in filename) and ('_5p' in copyfile)) or (('3p_PCR' in filename) and ('_3p' in copyfile))):
     subprocess.call(mv_cmd,shell=True)


  print "\nDone with preprocessing iteration!";

  for smpl in samples:
   these_files_r1 = glob.glob("%s%s*.fully_trimmed.fastq.gz" % (TRIMDIR, smpl))
   these_files_r1.sort()
   if (smpl in SampleToIllumIndx) and (smpl in SampleTo5pIndx) and (smpl in SampleTo3pIndx):
    all_file = "%s%s_orfEnd_All_illumIdx_%s_myIdx_%s.%s.fully_trimmed.fastq.gz" %(TRIMDIR, smpl, SampleToIllumIndx[smpl], SampleTo3pIndx[smpl], SampleTo5pIndx[smpl]);
    if all_file in these_files_r1:
     these_files_r1.remove(all_file)
    file_str_1 = " ".join(these_files_r1)
    zcat_cmd = "zcat %s | gzip > %s%s_orfEnd_All_illumIdx_%s_myIdx_%s.%s.fully_trimmed.fastq.gz" % (file_str_1, TRIMDIR, smpl, SampleToIllumIndx[smpl], SampleTo3pIndx[smpl], SampleTo5pIndx[smpl])
    print zcat_cmd;
    subprocess.call(zcat_cmd, shell=True)

  print "\nDone with preprocessing iteration!";
  print;

#################################################################################################################################################################################################

def return_sam_file(sample_name, prime_str, illumIndx, myIndx):
  trimmedfile = return_copy_file(sample_name, prime_str, illumIndx, myIndx);
  print trimmedfile;
  if (prime_str == '3p'):
   strimmedfile = string.replace(trimmedfile, TRIMDIR, SAMDIR_3P);
  elif (prime_str == '5p'):
   strimmedfile = string.replace(trimmedfile, TRIMDIR, SAMDIR_5P);
  else:
   strimmedfile = string.replace(trimmedfile, TRIMDIR, SAMDIR);
  rp_suffix = ".fully_trimmed" + COMMON_SUFFIX;
  return string.replace(strimmedfile, rp_suffix, ".sam");


def run_bowtie(trimmedfile):
  rp_suffix = ".fully_trimmed" + COMMON_SUFFIX;

  strimmedfile = string.replace(trimmedfile, TRIMDIR, SAMDIR);
  print strimmedfile; 
  samoutput = string.replace(strimmedfile, rp_suffix, ".sam");
  logoutput = string.replace(strimmedfile, rp_suffix, ".sam.log");
  erroutput = string.replace(strimmedfile, rp_suffix, ".sam.err");

  cmdline = 'bowtie2 --very-sensitive -N 1 --norc -x %s -U %s -S %s > %s 2> %s' % (LIBRARY_INDEX,  trimmedfile, samoutput, logoutput, erroutput);
  ##  http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

  print cmdline
  print;
  subprocess.call(cmdline, shell=True)
  print cmdline


def run_postprocess_sam(trimmedfile):
  rp_suffix  = ".fully_trimmed" + COMMON_SUFFIX;
  strimmedfile = string.replace(trimmedfile, TRIMDIR, SAMDIR);
  samoutput  = string.replace(strimmedfile, rp_suffix, ".sam");
  bamoutput  = string.replace(strimmedfile, rp_suffix, ".bam");
  srtoutput  = string.replace(strimmedfile, rp_suffix, ".sorted");
  srtoutput2 = string.replace(strimmedfile, rp_suffix, ".sorted.bam");
  bedoutput  = string.replace(strimmedfile, rp_suffix, ".bed");
  cntoutput  = string.replace(strimmedfile, rp_suffix, ".cnt");
  rngoutput  = string.replace(strimmedfile, rp_suffix, ".range.cnt");

  cmd_1 = 'samtools view -F 0x04 -Sb %s > %s' % (samoutput, bamoutput);
  ## http://www.htslib.org/doc/samtools.html
  cmd_2 = 'samtools sort %s %s' % (bamoutput, srtoutput);
  cmd_3 = 'bamToBed -i %s > %s' % (srtoutput2, bedoutput);
  cmd_4 = 'cut -f 1 %s | sort | uniq -c | awk \'{print $2 \"\\t\" $1}\' | sort -k2nr > %s'  % (bedoutput, cntoutput);
  cmd_5 = 'cut -f 1,2,3 %s | sort | uniq -c | awk \'{print $2 $3 \"_\" $4  \"\\t\" $1}\' | sort -k2nr > %s' % (bedoutput, rngoutput);


  print cmd_1
  print;
  subprocess.call(cmd_1, shell=True)
  print cmd_2
  print;
  subprocess.call(cmd_2, shell=True)
  print cmd_3
  print;
  subprocess.call(cmd_3, shell=True)
  print cmd_4
  print;
  subprocess.call(cmd_4, shell=True)
  print cmd_5
  print;
  subprocess.call(cmd_5, shell=True)

  cp_cmd = "cp %s %s" %(cntoutput, CNTDIR);
  print cp_cmd
  print; 
  subprocess.call(cp_cmd,shell=True)

  cp_cmd = "cp %s %s" %(rngoutput, CNTDIR);
  print cp_cmd
  print; 
  subprocess.call(cp_cmd,shell=True)



#################################################################################################################################################################################################

def read_counts_file(cntfile):
  utrToCnt = {};
  with open(cntfile,'r') as tsvin:
    tsvin = csv.reader(tsvin, delimiter='\t')
    #headers = tsvin.next()
    for row in tsvin:
      utr = row[0];
      cnt = int(row[1]);
      utrToCnt[utr] = cnt;
  return utrToCnt;



def read_counts_sam_file(samffile, prime_str):
  utrToCnt = {};
  samfile = pysam.AlignmentFile(samffile,'r');
  a = 0;
  m = 0;
  u = 0;
  if (prime_str == '3p'):
   dtcPos = dtcPos_3P;
  elif (prime_str == '5p'):
   dtcPos = dtcPos_5P;

  for al in samfile.fetch():
   if (al.pos >=0):
    name = al.reference_name;
    pos1 = dtcPos['Pos1'][name];
    pos2 = dtcPos['Pos2'][name]-1;
    accept_cond = (abs(al.pos-pos1) < 2) or (abs(al.pos-pos2) < 2);

    qual = al.mapping_quality;
    if (not accept_cond):
     m = m+1;
    else:
     a = a+1;
     utr = name;
     if utr in utrToCnt:
      utrToCnt[utr] = utrToCnt[utr] + 1;
     else:
      utrToCnt[utr] = 1;
   else:
    u = u+1;
  print "In file %s there are:\n\t%d aligned\t (%.1f%%)\n\t%d mis-aligned\t (%.1f%%)\n\t%d not aligned to ORF ends\t (%.1f%%)\n\t%d ORFs found in the file\n" %(samffile, a, (100.*a/(a+u+m)), m, (100.*m/(a+u+m)), u, (100.*u/(a+u+m)), len(utrToCnt.keys()))

  return utrToCnt;



#################################################################################################################################################################################################
#################################################################################################################################################################################################




if (__name__ == "__main__"):

  ##############################
  ## preprocess each data set
  ##############################


  if(DO_TRIM_ADAPTERS):
   pool1=Pool(processes=24)
   for i in range(len(in_prefixes)):
     inpr = in_prefixes[i]
     print "beginning %s" % (inpr)
     apply(preprocess_gz_file, [inpr])
   pool1.close()
   pool1.join()

  files_for_alignemnt = set();
  for key, sample_name in SAMPLE_MAP.items():
   illumIndx = key[0];
   myIndx    = key[1];
   prime_str = PRIME_MAP[(illumIndx, myIndx)];
   fcp = return_copy_file(sample_name, prime_str, illumIndx, myIndx);

   print fcp;
   if os.path.exists(fcp):
    files_for_alignemnt.add(fcp);
    FileToPrimeEnd[fcp] = prime_str;
    print fcp;
   if sample_name in SampleToAllIndx:
    fcp1 = return_copy_file(sample_name, 'All', illumIndx, SampleToAllIndx[sample_name]);
    print "AAA\t" + fcp1;
    if os.path.exists(fcp1):
     files_for_alignemnt.add(fcp1);
     FileToPrimeEnd[fcp1] = 'All';

  for afile in files_for_alignemnt:
   if (FileToPrimeEnd[afile] == 'All'):
    DO_END = 'ALL';
    LIBRARY_INDEX = LIBRARY_INDEX_ALL
    SAMDIR = SAMDIR_ALL
   if (FileToPrimeEnd[afile] == '3p'):
    DO_END = '3P';
    LIBRARY_INDEX = LIBRARY_INDEX_3P
    SAMDIR = SAMDIR_3P
   if (FileToPrimeEnd[afile] == '5p'):
    DO_END = '5P';
    LIBRARY_INDEX = LIBRARY_INDEX_5P
    SAMDIR = SAMDIR_5P

   print afile;
   print FileToPrimeEnd[afile];

   if (DO_ALIGN_SEQUENCES):
    print "\n\nDO_ALIGN_SEQUENCES\n\n";
    run_bowtie(afile);
   if (DO_PROCESS_SAMFILES):
    print "\n\nDO_PROCESS_SAMFILES\n\n";
    run_postprocess_sam(afile);

#################################################################################################################################################################

  # count directly from the sam files and filter out reads that don't start exactly at the 5' and 3' ends
  if(DO_NEW_MERGE_FILES):
   files_for_counts = set();
   utrs_in_assay    = set();
   counts_dict = {};
   print "\n\nDO_NEW_MERGE_FILES\n\n";
   for key, sample_name in SAMPLE_MAP.items():
    illumIndx = key[0];
    myIndx    = key[1];
    sample       = SAMPLE_MAP[(illumIndx, myIndx)];
    prime_str    = PRIME_MAP[(illumIndx, myIndx)];
    fcp = return_sam_file(sample_name, prime_str, illumIndx, myIndx);
    if os.path.exists(fcp):
     files_for_counts.add(fcp);
     c_dict = {};
     #print "read file\t" +fcp;
     c_dict = read_counts_sam_file(fcp, prime_str);

     #print "process sample\t" + sample + "-" + prime_str;
     c_dict_items = c_dict.items();
     for utr, cnts in c_dict_items:
       naked_utr = string.replace(utr, "-3p", "");
       naked_utr = string.replace(naked_utr, "-5p", "");
       c_dict[naked_utr] = c_dict.pop(utr);
       utrs_in_assay.add(naked_utr);
     counts_dict[sample + "-" + prime_str] = c_dict; 
     #print c_dict;

   for sample in SampleToAllIndx:
    samp3 = sample + "-3p"
    samp5 = sample + "-5p"
    if ( (samp3 in counts_dict) & (samp5 in counts_dict)):
     x = counts_dict[sample + "-3p"];
     y = counts_dict[sample + "-5p"];
     counts_dict[sample + "-All"] = { k: x.get(k, 0) + y.get(k, 0) for k in set(x) | set(y) }
     print "%d ORFs found in sample %s\n" %(len(counts_dict[sample + "-All"].keys()), sample)

   writefile  =  'new_strict_counts' + finalTag +  '.tsv';
   writefile3 =  'new_strict_columns' + finalTag +  '.tsv';
   sum_c_dict = {};
   sum_r_dict = {};
   ic = 1;
   with open(writefile3, 'w') as outfile3:
     outfile3.write("%i\tORF_id\n" % ic);
     ckeys = counts_dict.keys();
     ckeys.sort();
     for sample_key in ckeys:
      ic = ic +1;
      outfile3.write("%i\t%s\n" % (ic, sample_key));
     outfile3.close();

   with open(writefile,  'w') as outfile:
      outfile.write("ORF_id");
      ckeys = counts_dict.keys();
      ckeys.sort();
      for sample_key in ckeys:
       outfile.write("\t" + sample_key);
       sum_c_dict[sample_key] = sum(counts_dict[sample_key].values())/1000000.0
      outfile.write("\n");
      for utr in utrs_in_assay:
       outfile.write(utr);
       for sample_key in ckeys:
        if utr in counts_dict[sample_key]:
         outfile.write("\t" + np.str(counts_dict[sample_key][utr]));
        else:
         outfile.write("\t0");
       outfile.write("\n");
      outfile.close(); 


  exit();

