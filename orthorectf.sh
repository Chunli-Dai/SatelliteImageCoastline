#/bin/sh -f

#multidir='/data1/pgc_projects/coastline/imagery/';
multidir='/data1/pgc_projects/dai_aleutians_multi_mono/imagery/';
odir='orthorectwork/'

for line in `ls $multidir/*/*.ntf `; do
ntffile=$line;
tiffile=`echo $line | sed -e "s/.ntf/.tif/g"`
file=$(basename "$tiffile")  #get filename without path
tiffile=$odir/$file
if [ -f "$tiffile" ]
then
	echo "$tiffile found."
else
time gdalwarp -rpc -et 0.01 -co tiled=yes -co compress=lzw -t_srs EPSG:3413 $ntffile $tiffile
fi
#exit
done
