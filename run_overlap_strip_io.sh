#/bin/sh -f
#Finding overlapping DEMs using meta.txt files
# example: ./run_overlap_strip_io.sh '/*/ArcticDEM/region*/strips/2m/'

#res=8;
#samples_cloudsAndWater/WV02_20150813_10300100473E1900_103001004613E200_500518774050_01_P003_500446924060_01_P003_8_meta.txt
#for line in `cat $dir/out`; do
#for line in `cat Antarcticxml`; do # wrong to use xml files
#for line in `cat Reg01metalist`; do
dir=$1
#find /*/ArcticDEM/region*/strips/2m/ -type f -name '*meta.txt' > Reglist
find $dir -type f -name '*meta.txt' > Reglist
for line in `cat Reglist`; do
#echo $line 
#infile=$dir/$line
infile=$line
xr=`awk '$1!="X:"{next};!i++{min=$2;max=$2;}{for(j=2;j<=NF;++j){min=(min<$j)?min:$j;max=(max>$j)?max:$j}}END{printf " %d\n", max}' $infile`
xl=`awk '$1!="X:"{next};!i++{min=$2;max=$2;}{for(j=2;j<=NF;++j){min=(min<$j)?min:$j;max=(max>$j)?max:$j}}END{printf " %d\n", min}' $infile`
yu=`awk '$1!="Y:"{next};!i++{min=$2;max=$2;}{for(j=2;j<=NF;++j){min=(min<$j)?min:$j;max=(max>$j)?max:$j}}END{printf " %d\n", max}' $infile`
yl=`awk '$1!="Y:"{next};!i++{min=$2;max=$2;}{for(j=2;j<=NF;++j){min=(min<$j)?min:$j;max=(max>$j)?max:$j}}END{printf " %d\n", min}' $infile`
#echo x range: $xl $xr '; y range:' $yl $yu
#if [ $xl -ge $xl0 -a $xl -le $xr0 -o $xr -ge $xl0 -a $xr -le $xr0 ] ; then
#  if [ $yu -ge $yl0 -a $yu -le $yu0 -o $yl -ge $yl0 -a $yl -le $yu0 ] ; then
#       echo overlapping map $infile
#	demfile=${infile/meta.txt/dem_browse.tif}
#	cp $demfile $odir/
#	echo $xl $xr $yl $yu $infile >> boundaries_Anta.dat 
#	echo $xl $xr $yl $yu $infile >> boundaries_reg01.dat 
#	echo $xl $xr $yl $yu $infile >> boundaries_reg31.dat 
	echo $xl $xr $yl $yu $infile >> boundaries_regall_strip.dat
#  fi
#fi

done

