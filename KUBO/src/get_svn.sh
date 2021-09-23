#!/bin/sh
# script to create subversion info file (Fortran code) for LSMS2

OUTFILE="test_svn.F90"
if [ $# -gt 0 ] ; then OUTFILE=$1; fi

UNIT="myunit"
MYUSER=`whoami 2> /dev/null`
SYSTEMNAME=`hostname -f 2> /dev/null`
#DATE=`date +"%D %T" 2> /dev/null`
DATE=`date +"%A, %B %d %Y %T" 2> /dev/null`

rm -rf $OUTFILE

svn info 2> /dev/null | awk \
-v OUTFILE="$OUTFILE" \
-v UNIT="$UNIT" \
-v SQ="'" \
-v MYUSER="$MYUSER" \
-v SYSTEMNAME="$SYSTEMNAME" \
-v DATE="$DATE" '
BEGIN{
  FMT1=sprintf("   write(%s,%s(\"\# \",t5,a,t30,\": \",a)%s)",UNIT,SQ,SQ)
  FMT2=sprintf("   write(%s,%s(\"\# \",t5,a)%s)",UNIT,SQ,SQ)
  FMT3=sprintf("   write(%s,%s(a)%s)",UNIT,SQ,SQ)
  FMT4=sprintf("   write(%s,%s(24x,a)%s)",UNIT,SQ,SQ)
  printf("! \n") >>OUTFILE
  printf("! this routine is automatically generated at compilation time\n") >>OUTFILE
  printf("! DO NOT EDIT \n") >>OUTFILE
  printf("! \n") >>OUTFILE
#  printf("subroutine print_SVN\n") >>OUTFILE
#  printf("implicit none\n") >>OUTFILE

#  printf("%s repeat(\"-\",80)\n",FMT3)>>OUTFILE
#  printf("%s trim(\"*********************************\")\n",FMT4)>>OUTFILE
#  printf("%s trim(\"*     Output from print_svn     *\")\n",FMT4)>>OUTFILE
#  printf("%s trim(\"*********************************\")\n",FMT4)>>OUTFILE
#  printf("%s repeat(\"=\",80)\n",FMT3)>>OUTFILE
 
  svnfound=0
}

{ # main part
   if(match($0,"^Repository Root:"))
   {printf ("%s \"Repository Root\", trim(\"%s\") \n",FMT1,$3) >> OUTFILE}


   if(match($0,"^URL:"))
   {printf ("%s \"Repository URL\", trim(\"%s\") \n",FMT1,$2) >> OUTFILE}

   if(match($0,"^Revision:"))
   {svnfound=1
    printf ("%s \"Revision\", trim(\"%s\") \n",FMT1,$2) >> OUTFILE}

   if(match($0,"^Last Changed Rev:"))
   {printf ("%s \"Last Changed Revision\", trim(\"%s\") \n",FMT1,$4) >> OUTFILE}

   if(match($0,"^Last Changed Date:"))
   { fx=""; for (i=4; i<=NF; i++) fx=fx$i" ";
     printf ("%s \"Last Changed Date\", trim(\"%s\") \n",FMT1,fx) >> OUTFILE}
}
END {

    if (svnfound = 0) 
    {msg="subversion information is not available"
     printf("%s trim(\"%s\")\n", FMT2,msg) >> OUTFILE} 
   
    if (MYUSER) 
     {printf("%s \"username\",trim(\"%s\")\n",FMT1,MYUSER) >> OUTFILE }
#    else  
#     {us="unknown"; printf("%s \"username\",trim(\"%s\")\n",FMT1,us) >> OUTFILE }

    if (SYSTEMNAME) 
     {printf("%s \"system name\",trim(\"%s\")\n",FMT1,SYSTEMNAME) >> OUTFILE }
#    else  
#     {sys="unknown"; printf("%s \"system name\",trim(\"%s\")\n",FMT1,sys) >> OUTFILE }

    if (DATE)
     {printf("%s \"code compiled\",trim(\"%s\")\n",FMT1,DATE) >> OUTFILE }

#   printf("%s \"\# \",repeat(\"-\",79)\n",FMT3)>>OUTFILE

#    printf("end subroutine print_SVN\n") >>OUTFILE
 
}'