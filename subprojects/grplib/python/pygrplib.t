#! /bin/sh

# May 3, 2007

# Portions of the script which must be customized are marked with "!!",
# below.

# The complete list of tests must be placed in "alltests", !!4, below.
# The test[s] for the SHORTTEST must be placed in "shortlist", !!5 below.

# pygrplib.t (this file)
# driver script (unix sh) for myglue package testing

# pygrplib.t.py (external file)
# Python testing script to be activated by the driver above
#

# syntax:
# pygrplib.t [<testid> ... ]

# interface between python and unix-sh (this file) is done through
# "export INFILE" and  "export OUTFILE", see !!6, !!7, !!8

# Test suite activation may be made through  setting and exporting
# INFILE, see !!7, !!8

# Developer's comments:
# Limited to ascii file output at the time being

######################################################################
# Initialization

# !!3
gluename="pygrplib"

# set up list of tests
# !!4
alltests="test1 test2 test3 test4 test5 test6 test7 test8 test9 test10 test11 test12"

# "short" test to run
# !!5
shortlist="test1"


# compute date string for log file
DT=`date +'%d%b%Y_%T'`


# convenience definitions
OUTDIR=$TESTOUT/$gluename
SAVDIR=$TESTSAV/$gluename
INDIR=$TESTIN/$gluename
LOGDIR=$TESTLOG/$gluename

# set up log file name
LOGFILE=$LOGDIR/${gluename}_log.$DT

#get rid of old logs
rm -f $LOGDIR/${gluename}_log.*

# Any tests specified on command line?
if test $# -gt 0; then
  # yes, do those tests
  testlist=$*
else
  # No, see if we are to do "short" test
  if test "x$SHORTTEST" = "x" ; then
    # No, do everything
    testlist=$alltests
  else
    # yes, do short test
    testlist=$shortlist
  fi
fi

# Make sure we have a log directory
if test -d $LOGDIR ; then
 :
else
  mkdir -p $LOGDIR
  if test $? -ne 0 ; then
    error_exit ""
  fi
fi


# Make sure we have an output directory
if test -d $OUTDIR ; then
 :
else
  mkdir -p $OUTDIR >> $LOGFILE 2>&1
  if test $? -ne 0 ; then
    error_exit "can't create output directory $OUTDIR"
  fi
fi

# check for directory environment variables
if test "x${TESTIN}" = "x" -o "x${TESTOUT}" = "x" -o "x${TESTSAV}" = "x" \
   -o "x${TESTLOG}" = "x" ; then
  error_exit "one or more of TESTIN/TESTOUT/TESTSAV/TESTLOG not defined"
fi

# announce ourselves
echo ""
echo "${gluename} regression" | tee $LOGFILE
echo ""

script_succeeded=0

######################################################################
# Begin per-test loop

for testid in $testlist
do

  # delete old outputs
  rm -f $OUTDIR/${testid}*


  # Set up file names
  outfile=$OUTDIR/${testid}.txt
  savfile=$SAVDIR/${testid}.txt


# !!6
# export OUTFILE, INFILEs etc. as environmental variables
# which will be grabbed by python package called below
#
#
  OUTFILE=$outfile;
  export OUTFILE;


  echo "running $testid" >> $LOGFILE

  ####################################################################
  # run the tool
   case ${testid} in

    test1) command_str="python ${gluename}.t.py > /dev/null"
      infile_name="pygrplib_pha1.fits"
      INFILE=$INDIR/${infile_name};
      TESTID=${testid};
      export INFILE;
      export TESTID;
      echo $command_str  | tee -a $LOGFILE
      eval $command_str  | tee -a $LOGFILE  2>&1
             ;;
    test2) command_str="python ${gluename}.t.py > /dev/null"
      infile_name="pygrplib_pha1.fits"
      INFILE=$INDIR/${infile_name};
      TESTID=${testid};
      export INFILE;
      export TESTID;
      echo $command_str  | tee -a $LOGFILE
      eval $command_str  | tee -a $LOGFILE  2>&1
             ;;
    test3) command_str="python ${gluename}.t.py > /dev/null"
      infile_name="pygrplib_pha1.fits"
      INFILE=$INDIR/${infile_name};
      TESTID=${testid};
      export INFILE;
      export TESTID;
      echo $command_str  | tee -a $LOGFILE
      eval $command_str  | tee -a $LOGFILE  2>&1
             ;;
    test4) command_str="python ${gluename}.t.py > /dev/null"
      infile_name="pygrplib_pha1.fits"
      INFILE=$INDIR/${infile_name};
      TESTID=${testid};
      export INFILE;
      export TESTID;
      echo $command_str  | tee -a $LOGFILE
      eval $command_str  | tee -a $LOGFILE  2>&1
             ;;
    test5) command_str="python ${gluename}.t.py > /dev/null"
      infile_name="pygrplib_pha1.fits"
      INFILE=$INDIR/${infile_name};
      TESTID=${testid};
      export INFILE;
      export TESTID;
      echo $command_str  | tee -a $LOGFILE
      eval $command_str  | tee -a $LOGFILE  2>&1
             ;;
    test6) command_str="python ${gluename}.t.py > /dev/null"
      infile_name="pygrplib_pha1.fits"
      INFILE=$INDIR/${infile_name};
      TESTID=${testid};
      export INFILE;
      export TESTID;
      echo $command_str  | tee -a $LOGFILE
      eval $command_str  | tee -a $LOGFILE  2>&1
             ;;
    test7) command_str="python ${gluename}.t.py > /dev/null"
      infile_name="pygrplib_pha1.fits"
      INFILE=$INDIR/${infile_name};
      TESTID=${testid};
      export INFILE;
      export TESTID;
      echo $command_str  | tee -a $LOGFILE
      eval $command_str  | tee -a $LOGFILE  2>&1
             ;;
    test8) command_str="python ${gluename}.t.py > /dev/null"
      infile_name="pygrplib_pha1.fits"
      INFILE=$INDIR/${infile_name};
      TESTID=${testid};
      export INFILE;
      export TESTID;
      echo $command_str  | tee -a $LOGFILE
      eval $command_str  | tee -a $LOGFILE  2>&1
             ;;
    test9) command_str="python ${gluename}.t.py > /dev/null"
      infile_name="pygrplib_pha1.fits"
      INFILE=$INDIR/${infile_name};
      TESTID=${testid};
      export INFILE;
      export TESTID;
      echo $command_str  | tee -a $LOGFILE
      eval $command_str  | tee -a $LOGFILE  2>&1
             ;;
    test10) command_str="python ${gluename}.t.py > /dev/null"
      infile_name="pygrplib_pha1.fits"
      binfile_name="binfile.fits"
      INFILE=$INDIR/${infile_name};
      TESTID=${testid};
      BINFILE=$INDIR/${binfile_name};
      export INFILE;
      export TESTID;
      export BINFILE;
      echo $command_str  | tee -a $LOGFILE
      eval $command_str  | tee -a $LOGFILE  2>&1
             ;;
    test11) command_str="python ${gluename}.t.py > /dev/null"
      infile_name="testin12.fits"
      INFILE=$INDIR/${infile_name};
      TESTID=${testid};
      export INFILE;
      export TESTID;
      echo $command_str  | tee -a $LOGFILE
      eval $command_str  | tee -a $LOGFILE  2>&1
             ;;
    test12) command_str="python ${gluename}.t.py > /dev/null"
      infile_name="testin13.fits"
      INFILE=$INDIR/${infile_name};
      TESTID=${testid};
      export INFILE;
      export TESTID;
      echo $command_str  | tee -a $LOGFILE
      eval $command_str  | tee -a $LOGFILE  2>&1
             ;;
  esac


  ####################################################################
  # check the outputs

  # Init per-test error flag
  mismatch=1

  # if different tests need different kinds of comparisons, use a
  #  case ${testid} in...  here

  ######################################################################
  # ascii files
  #
   diff $OUTDIR/${testid}.txt $SAVDIR/${testid}.txt_std >/dev/null 2>>$LOGFILE
   if  test $? -ne 0 ; then
     echo "ERROR: TEXT MISMATCH in $OUTDIR/${testid}.txt" >> $LOGFILE
     mismatch=0
   fi

  ####################################################################
  # Did we get an error?
  if test $mismatch -eq 0 ; then
    # Yes
    echo "${testid} NOT-OK"
    script_succeeded=1
  else
    # No
    echo "${testid} OK"
  fi

done
# end pre-test loop
######################################################################


######################################################################
# report results

# blank line
echo ""

if test $script_succeeded -eq 0; then
    echo "${gluename} : PASS" | tee -a $LOGFILE
else
    echo "${gluename} : FAIL" | tee -a $LOGFILE
fi

echo "log file in ${LOGFILE}"

exit $script_succeeded
