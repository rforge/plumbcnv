#!/bin/bash
  
#######################################################################
## Bash Script for Extracting LRR/BAF data from Genome Studio Files  ##
#######################################################################

# By: Nicholas Cooper, JDRF/WT DIL, CIMR Cambridge, Nov 2012 #

### INSTRUCTIONS ###
# Firstly, in the directory this script is run from, make sure there is 
# the following file:  'file.spec.txt' and that the default directory structure has
# been copied to the 'outdir' location. set 'loco' to source of raw data and
# ensure bim/vcf file is specified and that the column numbers match.
# decide whether to run call rate QC in plink or R. Assumes plink installed and command is in path
# Check and modify the user options below to specify locations of columns of data.
# Finally, run this script in the terminal.
####################


### DEFAULT USER OPTIONS - MODIFY TO SUIT ###
#
# set output directory:
#
    outdir=/chiswick/data/ncooper/metabochipCNVanalysis2012/   #../../
#
# set directory of source files (typically read-only)
#
    #loco=/ipswich/data/Immunochip/FinalReports/
    loco=/chiswick/data/store/metabochip/FinalReports/
#
# Need to create a file called file.spec.txt which contains each file name with variable locations
# in the genome studio file (or equivalent), spec which columns are sample ids, and which snp ids:
# and if using plink for SNP-QC then which are the allele columns. Headers must be as below.
# See example of file format for 3 input files, mix of zipped and raw text, plus one file with 
# different column locations for the data. Raw files must be in the dir: 'loco'
#
#                FILE GRP TYPE SAMP SNP A1 A2 LRR BAF
# sanger-controls.txt  1  txt    1   2  6  7  13  12
#    uva-controls.txt  2  txt    1   2  5  6  12  11
#    t1d-cases.txt.gz  3  gzip   1   2  6  7  13  12
#
# set maketemplate to "yes" to make a template for this file.spec.txt which
# user can manually check modify against the data desired to import 
#
   maketemplate="no"
#
# Maximum number of SNPs expected per sample [the closer to the actual # makes it faster]
# Please err on the side of too many, e.g, 2000000 is a safe default in the year 2012
#
    MaxSnp=2000000
#
# which data to be extracted 
# e.g, 'BAF' or 'LRR'
#
    fn='LRR'
#
# set getSnpSub= "yes" or "no"; whether to obtain list of markers and subject ids
# [must do this at least the first time]
#
    getSnpSub="no"   #yes
#
# set lgen= "yes" or "no"; whether to extract an LGEN file for SNP-QC in plink
#
    lgen="no"   #yes
#
# set fakefam= "yes" or "no"; whether to use existing real family file, or create a fake one
# [only used when lgen="yes"]
    fakefam="no"  #yes 
#
# set dolen= "yes" or "no"; whether to obtain file lengths for all input files
# this helps when importing the SNPs into snpMatrix to do SNP-QC
#
    dolen="no"  #yes
#
# set combin= "yes" or "no"; whether to keep multiple data files separate or combine them
#
    combin="no"
#
# skip main extraction skipmain="yes"
#
    skipmain="no"
#
# set max number of processors to use (default 1)
#
    maxproc=1
#
# set delseps= "yes" or "no"; whether to delete the separate files if combining
#
    delseps="no"
#
# set mapfromafile= "yes" or "no"; whether to create the .map from a file here
#
    mapfromafile="no"   #yes
#
# file name with chromosome, position and snp data for making a map file
# file name with the snp,chr and pos data in cols 2,3,4 #
  #  mfnm=$(echo $loco | sed 's/[/]$//g')/support.txt
  # be wary of trailing spaces!

    mfnm=$"/chiswick/data/store/metabochip/PLINK/Metabo_20100426_58C.bim"
#
# Further parameters for reading data for the map file:
  mSampCol=1  # column with sample label
  mSnpCol=2  # column with snp label
  mchrCol=3  # column with chr number
  mchrColB=1  # column with chr number in 'B' = bim file
  mposCol=4  # column with snp position
  mnexttype="txt"  # file with snp,chr,pos data is gzip or txt?
  gs="no"  # 'yes' if genome studio file or 'no' if a proper snp support file (e.g, bim)
#
# set rmv23= "yes" or "no"; whether to remove chromosome 23, X,Y,MT from map file
# probably should set to "yes" if plinkQC="yes", or "no" if plinkQC="no".
#
    rmv23="no"  #yes
#
# set plinkQC= "yes" or "no"; whether to run sample and snp qc (call rate, hwe) with plink
#
    plinkQC="no" #yes
#
# set plink call rate and HWE thresholds
#
    sampcr=.05
    snpcr=.05
    hwe=.00000001
#
# set your email address for notifications of the script finishing
#
    mymail=nick.cooper@cimr.cam.ac.uk
#
######### END USER OPTIONS ##########

### HELP TEXT ### 
USAGE="Script to convert genome studio files into long format files for fast R,Plink import.
Usage: getDataGS [-h] [-B] [-T] [-N] [-S] [-L] [-F] [-l] [-C] [-D] [-M] [-m] [-R] [-P] [-E] [-x] [-y] [-z]
    -h      Shows this help
    -B      Set maximum possible number of unique SNPs to search for [eg 2 million]
    -T      Set to LRR or BAF (which to import)
    -O      Location (directory) where your CNV folder structure starts. Output will go to /...DATA/ColumnData/ 
    -F      Location (directory) where the raw genome studio datafiles are
    -t      Alongside -F and -O can use to create a template for the essential 'file.spec.txt'
    -N      Maximum number of threads/cores to use
    -S      Generate list of SNPs and SAMPLES
    -L      Generate plink lgen file(s)
    -f      Generate a fake family file for plink (as opposed to having a real one prepared)
    -l      Make file with file lengths of input data
    -C      Combine multiple input file into a single file for output
    -D      Delete the separate files after merging into one file (-C)
    -X      Skip main data extraction (e.g, because you only want to generate annotation)
    -M      Create plink .map file from one of these files
    -m      Specify file name to use for creating the map file in (-M)
    -G      Whether the -m file is a genome studio (long) file
    -R      Remove non-autosome SNPs before doing sample callrate QC
    -P      Do sample and SNP QC in plink - requires less memory,disk usage than R (but may be slower)
    -E      Email address to send notification emails
    -x      Sample call rate for plink (default 0.05 = 95%)
    -y      Snp call rate for plink (default 0.05 = 95%)
    -z      Hardy Weinberg Equilibrium p threshold (default 0.00000001)
 Examples:
    make file.spec.txt template:  ./getDataGS.sh -t -F '/Raw/GenomeStudio' -O '/MyCNV'
    generate subject,snp lists:   ./getDataGS.sh -XS -F '/Raw/GenomeStudio' -O '/MyCNV/'
    import BAF:                   ./getDataGS.sh -T 'BAF' -N 3 -F '/Raw/GenomeStudio' 
    import LRR for QC in R:       ./getDataGS.sh -T 'LRR' -S -l -M -F '/Raw/GenomeStudio/'
    import LRR and do plink QC:   ./getDataGS.sh -SLMRPlf -T 'LRR' -F '/Raw/GenomeStudio' -m '/Raw/GenomeStudio/support.vcf' "

### parse command line options ###

while getopts B:T:N:E:m:F:O:x:y:hzfltLSCDMGRPX OPT;
do
    case $OPT in
    h)  echo "$USAGE"
        exit 0 ;;
    B)  MaxSnp=$OPTARG ;;
    T)  fn=$OPTARG ;;
    N)  maxproc=$OPTARG ;;
    F)  loco=$OPTARG ;;
    O)  outdir=$OPTARG ;;
    m)  mfnm=$OPTARG ;;
    E)  mymail=$OPTARG ;;
    x)  sampcr=$OPTARG ;;
    y)  snpcr=$OPTARG ;;
    z)  hwe=$OPTARG ;;
    S)  getSnpSub="yes" ;;
    L)  lgen="yes" ;;
    f)  fakefam="yes" ;;
    t)  maketemplate="yes" ;;
    l)  dolen="yes" ;;
    C)  combin="yes" ;;
    D)  delseps="yes" ;;
    X)  skipmain="yes" ;;
    M)  mapfromafile="yes" ;;
    G)  gs="yes" ;;
    R)  rmv23="yes" ;;
    P)  plinkQC="yes" ;;
    *)  echo "$USAGE" >&2
        exit 1 ;;
    esac
done

## remove final forward-slash from directories ##
loco=$(echo $loco | sed 's/[/]$//g')
outdir=$(echo $outdir | sed 's/[/]$//g')

##### MAKE A TEMPLATE FOR file.spec.txt ######

if [ "$maketemplate" = "yes" ]
then
 tt="\t"
 ls $loco > allfiles.TEMPLATE
 grep '.txt$' allfiles.TEMPLATE > txt.TEMPLATE
 grep '.gz$' allfiles.TEMPLATE > gz.TEMPLATE
 cnt=1
 while IFS= read -r file
 do
   echo -e $file"$tt"$cnt"$tt"txt"$tt"1"$tt"2"$tt"6"$tt"7"$tt"13"$tt"12 > line.$cnt.TEMPLATE
   cnt=$(( $cnt + 1 ))
 done < "txt.TEMPLATE" 

 while IFS= read -r file
 do
   echo -e $file"$tt"$cnt"$tt"gzip"$tt"1"$tt"2"$tt"6"$tt"7"$tt"13"$tt"12 > line.$cnt.TEMPLATE
   cnt=$(( $cnt + 1 ))
 done < "gz.TEMPLATE" 

 echo -e FILE"$tt"GRP"$tt"TYPE"$tt"SAMP"$tt"SNP"$tt"A1"$tt"A2"$tt"LRR"$tt"BAF > hdr.TEMPLATE

 cat hdr.TEMPLATE line.*.TEMPLATE > file.spec.txt
 rm *.TEMPLATE
 mv file.spec.txt "$outdir"/LRRDATA/file.spec.txt
 echo Template: file.spec.txt written to "$outdir"/LRRDATA/file.spec.txt
 echo Please check this file and edit where necessary before running the main import job.
 echo The column number for each variable will vary across genome studio exports.
 echo "The default is to assume each file is a different cohort (GRP) but this may not be true."
 exit
fi

##############################################


echo
echo Options in use are:
echo -------------------
echo MaxSnp: $MaxSnp
echo Datatype: $fn
echo Max number of cores: $maxproc
echo Filename with Chr,Pos,Snp support information: $mfnm
echo Generate Plink .Map file using support filename above?: $mapfromafile
echo SNP location support file is in genome studio format?: $gs
echo "Location of main CNV directory (output then to /...DATA/ColumnData/): $outdir "
echo Directory containing raw data files: $loco
echo Email for notifications: $mymail
echo Get SNP/Sample IDs?: $getSnpSub
echo Make Plink LGEN file?: $lgen
echo Make fake Plink family file?: $fakefam
echo Write file lengths to file?: $dolen
echo Combine all into 1 file?: $combin
echo If combining, delete separate files?: $delseps
echo Skip main extraction?: $skipmain
echo Remove non-autosomes prior to sample QC in Plink?: $rmv23
echo Do QC in Plink?: $plinkQC
echo


### Warn if other options don't make sense when Plink-QC is selected ###

if [ "$plinkQC" = "yes" ]
then
  if [ "$mapfromafile" = "no" -o "$lgen" = "no" -o "$fakefam" = "no" ]
  then
    echo Warning! Plink QC option is set, but at least one of the options -L, -M and -f, are not.
    echo Please ctrl-C and modify command if you do not already have .lgen, .map and .fam files
  fi
  if [ "$mfnm" = "$loco/support.txt" -o "$rmv23" = "no" -o "$dolen" = "no" -o $getSnpSub = "no" ]
  then
    echo Warning! Plink QC option is set, but at least one of the options -m, -R, -l and -S, are not.
    echo Please ctrl-C and modify command if you are not sure that these options are appropriate.
  fi
else
  echo Plink settings seem ok
fi

################################################################
#################       MAIN PROGRAM       #####################
################################################################

## Runs the rest of this script from the directory where you want the resulting files to go

cd "$outdir/$fn"DATA/ColumnData/

# IF BAF, make a copy of the file.spec from the LRR folder in the BAF folder
if [ "$fn" = "BAF" ]
then
 cp "$outdir/LRRDATA/file.spec.txt" file.spec.txt
 echo "In BAF mode plink options are reset to 'no'"
 lgen="no"
 fakefam="no"
 rmv23="no"
 plinkQC="no"
fi


## EXTRACT SAMPLE AND MARKER LISTS

# parse spec file
sed 1d ../file.spec.txt | sort -b -k1 > file.spec.temp
cut -f 1 file.spec.temp > FILES.TMP
cut -f 2 file.spec.temp > GRPS.TMP
cut -f 3 file.spec.temp > TYPES.TMP
cut -f 4 file.spec.temp > SAMPS.TMP
cut -f 5 file.spec.temp > SNPS.TMP
cut -f 6 file.spec.temp > A1.TMP
cut -f 7 file.spec.temp > A2.TMP
cut -f 8 file.spec.temp > LRR.TMP
cut -f 9 file.spec.temp > BAF.TMP
rm file.spec.temp

if [ "$fn" = "BAF" ]
then
    echo "Extracting BAF data from genome studio file"
    cat BAF.TMP > DATCOL.TMP
else
    echo "Extracting LRR data from genome studio file"
    cat LRR.TMP > DATCOL.TMP
fi

numfiles=$(wc -l FILES.TMP | cut -d ' ' -f 1)
echo "$numfiles datafiles expected"

if [ "$getSnpSub" = "yes" ] ; 
then 
    #GET LIST OF SNPs#
    # get snp list from an unzipped file if available, else from zipped
    echo 'extracting snp list...'
    nexttype=$(head -1 TYPES.TMP)
    SampCol=$(head -1 SAMPS.TMP)
    SnpCol=$(head -1 SNPS.TMP)
    fnm=$(head -1 FILES.TMP)
    if [ "$nexttype" = "gzip" ];
    then 
      SUBID1=$(zcat $loco/$fnm | head -50 | tail -1 | cut -f $SampCol)
      SNPID1=$(zcat $loco/$fnm | head -50 | tail -1 | cut -f $SnpCol)
      zcat $loco/$fnm | head -"$MaxSnp" | grep $SUBID1 | cut -f $SnpCol  > snplist.txt ;
    else
      SUBID1=$(head -50 $loco/$fnm  | tail -1 | cut -f $SampCol)  
      SNPID1=$(head -50 $loco/$fnm  | tail -1 | cut -f $SnpCol)
      head -"$MaxSnp" $loco/$fnm  | grep $SUBID1 | cut -f $SnpCol  > snplist.txt ;
    fi
    
    #GET LIST OF IDS (in order)#
    echo 'extracting ID list...'
    
    for i in $(seq 1 $numfiles);
    do
       nexttype=$(head -$i TYPES.TMP | tail -1 )
       SampCol=$(head -$i SAMPS.TMP | tail -1 )
       file=$(head -$i FILES.TMP | tail -1 )
       if [ "$nexttype" = "gzip" ];
       then 
        zcat $loco/$file | grep $SNPID1 | cut -f $SampCol  > $file.ids ;
       else
        grep $SNPID1 $loco/$file | cut -f $SampCol  > $file.ids
       fi
    done

    cat *.ids > subIdsALL.txt ;
    numids=$(wc -l subIdsALL.txt | cut -d ' ' -f 1)
    echo "found $numids subject IDs"
    numsnps=$(wc -l snplist.txt | cut -d ' ' -f 1)
    echo "found $numsnps snp IDs"
    if [ "$numids" -eq "0" ] ; then exit ; fi
    if [ "$numsnps" -eq "0" ] ; then exit ; fi
#    cp *.ids "$outdir"/ANNOTATION/
## unnecessary duplication? think about it
    mv *.ids "$outdir"/LRRDATA/ColumnIds/
    cp subIdsALL.txt "$outdir"/LRRDATA/ColumnIds/subIdsALL.txt
    mv subIdsALL.txt "$outdir"/ANNOTATION/SAMPLE_SORT/subIdsALL.txt
    mv snplist.txt "$outdir"/ANNOTATION/snpNames.txt;
fi    


## MAIN DATA EXTRACTION

# GETs 1 Variable (col 'n') from file and puts into giant vector file (each datapoint on a new line). 
# Such a file is very fast to read into R
# order will be all SNPs for one SUBJ, then next subject and so forth
# SNPS.by.SUBJs.MATRIX <- matrix(GC.dat,nrow=num.SNPs)
# colnames = subIdsALL; rownames = snplist.txt
# use script: 'plainVecFileToBigMatrix.R'

# takes about 35 mins on 'wootton' for 200K snps x 6K samps. 

if [ "$skipmain" = "no" ] ; 
then 
  echo "Assuming genome studio format, which means a file header, then a line containing the text [Data]"
  echo "followed by row of headers, then the actual long format data"
  for i in $(seq 1 $numfiles)
  do
     nexttype=$(head -$i TYPES.TMP | tail -1 )
     file=$(head -$i FILES.TMP | tail -1 )
     coln=$(head -$i DATCOL.TMP | tail -1 )
     echo "extracting $fn from column $coln of $loco/$file"

     if [ "$nexttype" = "gzip" ]
     then 
      zcat $loco/$file | sed -e '1,/\[Data\]/d' | sed 1d | cut -f $coln > $file.$fn.dat &
     else
      sed -e '1,/\[Data\]/d' $loco/$file | sed 1d | cut -f $coln > $file.$fn.dat &
     fi

     NPROC=$(($NPROC+1))
     if [ "$NPROC" -ge $maxproc ]
     then
        wait
        NPROC=0
     fi

  done
fi

# data for each separate file should now be in files called:  $file.$fn.dat
# where $file are the original file names and $fn is the type of data, eg. LRR.


## FILE LENGTHS ##

# get file lengths for use in other scripts
if [ "$dolen" = "yes" -o "$lgen" = "yes" ] ; 
then 
 echo "getting lengths of vector format data files, writing to file.lengths.txt"
 wc -l *.$fn.dat | sed '$d' | sort -b -k2 > $outdir/LRRDATA/file.lengths.txt ;
fi

## GENERATE LGEN FILE FOR SNP-QC IN PLINK ##

if [ "$lgen" = "yes" ] 
then 
  echo "extracting allele data for plink LGEN file..."
  echo "WARNING, any existing LGEN files in this directory will be erased" 
  for i in $(seq 1 $numfiles)
  do
       nexttype=$(head -$i TYPES.TMP | tail -1 )
       file=$(head -$i FILES.TMP | tail -1 )
       SampCol=$(head -$i SAMPS.TMP | tail -1 )
       SnpCol=$(head -$i SNPS.TMP | tail -1 )
       A1=$(head -$i A1.TMP | tail -1 )
       A2=$(head -$i A2.TMP | tail -1 )
       coln=$(head -$i DATCOL.TMP | tail -1 )
       echo "extracting id and allele data from columns $SampCol,$SnpCol,$A1,$A2 of $loco/$file"
       if [ "$nexttype" = "gzip" ]
       then
        zcat $loco/$file | sed -e '1,/\[Data\]/d' | sed 1d | cut -f "$SampCol,$SnpCol,$A1,$A2" > $file.lgen
       else
        sed -e '1,/\[Data\]/d' $loco/$file | sed 1d | cut -f "$SampCol,$SnpCol,$A1,$A2" > $file.lgen
       fi
  done
  mailx -s DONE_MAINLGENS $mymail < "."
  #### UP  TO  HERE ####
  wc -l *.lgen > filelensNEGL.txt
  ## combine into a single file
  cat *.lgen > combined.lgen
  ## generate fake family file, which just codes for the different files 1,2,..,n
  if [ "$fakefam" = "yes" ] 
  then
     #lll=$(wc -l subIdsALL.txt | cut -d ' ' -f 1)
     #printf '1\n%.0s' {$(seq 1 $lll)} > 1s.txt
     #printf '0\n%.0s' {$(seq 1 $lll)} > 0s.txt
     awk '{print $1, $1, 0, 0, 1, 1}' "$outdir"/ANNOTATION/subIdsALL.txt > snpdata.fam
     #paste subIdsALL.txt subIdsALL.txt 0s.txt 0s.txt 1s.txt 1s.txt  > snpdata.fam
     cut -f 1 combined.lgen > fams.txt
  fi
  ## write in proper lgen format
  paste fams.txt combined.lgen > snpdata.temp
  rm *.lgen
  mv snpdata.temp snpdata.lgen 
  rm fams.txt 0s.txt 1s.txt
  mailx -s DONE_RESTLGEN $mymail < "." 
fi


## COMBINE SEPARATE FILES INTO 1 BIG FILE ##

if [ "$combin" = "yes" ]  
then 
  echo "combining separate files into 1: $fn.dat..."
  cat *.$fn.dat > $fn.combined.dat 
  if [ "$delseps" = "yes" ]  
  then 
    num1=$(wc -l $fnm.$fn.dat | cut -d ' ' -f 1)
    numA=$(wc -l $fn.combined.dat | cut -d ' ' -f 1)
    #run this next line only if sure it has worked!
    if [ "$numA" -gt "$num1" ]
    then 
      rm *.$fn.dat 
      echo "deleting partial files..."  
      echo 
    fi
    mv $fn.combined.dat $fn.dat 
  fi
fi

rm *.TMP 

echo 'complete!'

mailx -s DONE $mymail < "."

# END main extraction #

# need to prepare own snpdata.map file, might need annotation file#

## if data is present can create map file for immunochip using chr and pos in GS file ##

if [ "$mapfromafile" = "yes" ]  
then
  mfnm=$(echo "$mfnm" | sed 's/ *$//g') # in case of trailing space(s)
  if [ "$gs" = "yes" ]
  then
  echo "Extracting SNP info from genome studio type file"
    echo "parsing $mnexttype file."
    echo "[assuming sample-id in $mSampCol, chr in col $mchrCol, pos in col $mposCol, snp-label in col $mSnpCol]"
    echo "if any of these column numbers are wrong please modify mSampCol, mchrCol, mposCol, mSnpCol in this script"
  if [ "$mnexttype" = "gzip" ];
    then 
        SUBID1=$(zcat $mfnm | head -50 | tail -1 | cut -f $mSampCol)
        zcat $mfnm | head -"$MaxSnp" | grep $SUBID1 | cut -f $mchrCol  > snpdata1.temp 
        zcat $mfnm | head -"$MaxSnp" | grep $SUBID1 | cut -f $mSnpCol  > snpdata2.temp 
        zcat $mfnm | head -"$MaxSnp" | grep $SUBID1 | cut -f $mposCol  > snpdata3.temp 
    else
        SUBID1=$(head -50 $mfnm  | tail -1 | cut -f $mSampCol)  
        head -"$MaxSnp" $mfnm  | grep $SUBID1 | cut -f $mchrCol  > snpdata1.temp 
        head -"$MaxSnp" $mfnm  | grep $SUBID1 | cut -f $mSnpCol  > snpdata2.temp 
        head -"$MaxSnp" $mfnm  | grep $SUBID1 | cut -f $mposCol  > snpdata3.temp 
    fi
  else
    echo "Extracting SNP support from support file, e.g, bim, vcf, map, etc."
    if [[ $mfnm = *.bim ]]
    then
     echo $(basename $mfnm) seems to be a bim file so changing chr column from $mchrCol to $mchrColB
     mchrCol=$mchrColB
    else
     echo '(not a bim file so using default column numbers for vcf/map)'
    fi
    echo "[assuming chr in col $mchrCol, pos in col $mposCol, snp-label in col $mSnpCol]"
    echo "if any of these column numbers are wrong please modify mchrCol, mposCol, mSnpCol in this script"
    cut -f $mchrCol $mfnm > snpdata1.temp
    cut -f $mSnpCol $mfnm > snpdata2.temp
    cut -f $mposCol $mfnm > snpdata3.temp
  fi
  paste snpdata1.temp snpdata2.temp snpdata3.temp > snpdata.map
  echo snpdata.map file created successfully
  cp snpdata.map $outdir/ANNOTATION/rawdata.map
  rm snpdata*temp
fi

if [ "$rmv23" = "yes" ]  
then
  echo "removing chromosome 23, XY, X, Y, MT from snpdata.map file so they are not used in plink sample callrate calculations"
  mv snpdata.map snpdata.incl.XY.MT.map
  awk '{if ($1 <= 22) print $1, $2, $3}' snpdata.incl.XY.MT.map > snpdata.map
fi

## run plink to get call rates for snps and samples ##

if [ "$plinkQC" = "yes" ]  
then
  echo "Plink will require files: snpdata.lgen, snpdata.map and snpdata.fam to run."
  ## see whether data files are in current working directory, may already have been moved
  if [ -e ./snpdata.lgen ]
  then
    echo "Found snpdata.lgen file, running plink from current directory"
    whichmode="local"
  else
    if [ -e "$outdir/CALLRATES/Plink/snpdata.lgen" ]
    then
      echo "Did not find file in current directory, but found in /CALLRATES/Plink/"
      echo "Attempting to run plink from there"
      whichmode="callrate"
      cd $outdir/CALLRATES/Plink/
    else
      echo "did not find plink .lgen file in either expected location, exiting script."
      exit
    fi
  fi

  plink --lfile snpdata  --map3 --missing --hardy --missing-genotype '-' --out snpdataout --noweb --geno "$snpcr" --hwe "$hwe" --mind "$sampcr"

  mv snpdataout.hwe snpdataout.hwe.messy
  grep ALL snpdataout.hwe.messy  > snpdataout2.hwe
  head -1 snpdataout.hwe.messy > snpdataout1.hwe
  cat snpdataout1.hwe snpdataout2.hwe > snpdataout.hwe
  rm snpdataout1.hwe snpdataout2.hwe snpdataout.hwe.messy
  # plink codes missing as spaces which causes problems for R import
  # replace with '-' symbol
  sed -i 's/ALL          /ALL    -    -/g' snpdataout.hwe
  sed -i 's/ALL         A/ALL    -    A/g' snpdataout.hwe
  sed -i 's/ALL         C/ALL    -    C/g' snpdataout.hwe
  sed -i 's/ALL         G/ALL    -    G/g' snpdataout.hwe
  sed -i 's/ALL         T/ALL    -    T/g' snpdataout.hwe
  # sometimes other characters than [A,C,G,T] used for alleles, replace any like this with X
  sed -i 's/ALL         [A-Z]/ALL    -    X/g' snpdataout.hwe

  if [ "$whichmode" = "local" ]
  then
    echo "moving output files to /CALLRATES/Plink/"
    mkdir $outdir/CALLRATES/Plink/
    mv snpdata* $outdir/CALLRATES/Plink/
    echo "complete"
  else
    echo "Complete"
  fi

fi


##### script to get plate support files into annotation directory ####
# CTRL=/ipswich/data/Immunochip/support/controls
# CASS=/ipswich/data/Immunochip/support/t1d-cases
# cd /chiswick/data/ncooper/ImmunochipReplication/ANNOTATION/
# cut -f 1,3-8 $CTRL/control-sample-lookup-2012-05-21.tab > control-plate.txt
# sed '1d' $CASS/case-sample-lookup-2012-05-16.tab > case-plate.txt
# cat control-plate.txt case-plate.txt > plate.lookup.txt


# id <- readLines(paste(dir$ano,sample.fn,sep=""))
# well <- substr(IDs,8,10)
# plate <- substr(IDs,1,6)
# lookup <- (cbind(id,plate,well))
# write.table(lookup,file=paste(dir$ano,"plate.lookup.txt",sep=""),row.names=F,col.names=T,sep="\t")
