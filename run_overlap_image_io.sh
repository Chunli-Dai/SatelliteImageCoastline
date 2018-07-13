#/bin/sh -f
#Finding overlapping DEMs using meta.txt files
# example: ./run_overlap_image_io.sh '/data1/pgc_projects/dai_aleutians_multi_mono/imagery/WV*/'
# example: ./run_overlap_image_io.sh '/data1/pgc_projects/dai_aleutians_multi_mono/imagery/[G-Q]*/'

#/data1/pgc_projects/dai_aleutians_multi_mono/imagery//WV02/WV02_20130923221601_1030010026BD6B00_13SEP23221601-M1BS-500127380110_01_P001.xml 
dir=$1
#find /data1/pgc_projects/dai_aleutians_multi_mono/imagery/WV*/ -type f -name '*.xml' > Reglist
find $dir -type f -name '*.xml' > Reglist

