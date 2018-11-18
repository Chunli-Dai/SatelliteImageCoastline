%directory
%General directory that contains the mosaic tile DEM files, such as /elev/dem/setsm/ArcticDEM/mosaic/v2.0/
macdir=[]; % Leave it blank for linux;Absolute main directory part.
%macdir='/Users/chunlidai/surge/';
tiledir=[macdir,'/data/chunli/coastline/'];%ArcticDEM mosaic tile directory. 
stripdir='/*/ArcticDEM/region*/strips/2m/'; %directory of strip files
%multidir=[macdir,'/data1/pgc_projects/dai_aleutians_multi_mono/imagery/WV*/']; % directory of mono multispectral images
%multidir=[macdir,'/data1/pgc_projects/coastline/imagery/']; %directory of multispectral images; cover all images in the targed zone
multidir=[macdir,'/data1/pgc_projects/coastline/results/iceland/multi_reykjavik/'];
probdir='./output/'; % the output directory that contains the output of CoastTileMonoMain.m (probability results and image boundary results).

%control parameters
resr=2;
width=10e3; %buffer width of the a priori coastline, e.g., 10 km (best).
widthstat=2e3; %buffer width of a priori coastline for calculating ndwi statistics, and the buffer of tile for calculation.
width2=100.;   %get rid of a long thin (100 m width) beach band for NDWI difference
cloudarea=100*100; % clouds over water area, minimum size 100m by 100m; if smaller, clouds get ignored. Also smallest island that would be kept.
lakearea=1000*500*4;%smallest water body (or lake) (m^2) that would be kept;
almt=1e3*1e3;  %2e3*2e3;%(lost of data) %1e3*1e3 (suggested value); %minimum areas for each piece
cloudflag=1; %1 apply cloud detection; 0 do not apply cloud detection (four times faster).
novlmt=3; %if number of repeats <= novlmt, set the area as edges/void.
novmax=60; % if number of repeats > novmax, only select novmax of them. Suggest value 60
cntmin=25*25; %unit:pixels. The size of a priori land/ocean area should be big enough to ensure reliable statistical analysis of the histogram of the region. (Liu and Jezek, 2004)

%control parameters for multispec.m
threshold=0.5; %Suggest value 0.5 or 0.3; % general threshold of NDWI for water classification.
probthre=50.;% threshold for water probability. %Best for migitating random coregistration offset: 50.
stdthres=0.5; % if NDWI STD > stdthres, discard the image. Suggest value 0.5
dmthres=0.6; % if mean_ocean - mean_land > dmthres, discard the image. Suggest value 0.6

% Revise the parent directory for the following two lines.
orthworkdir=[macdir,'/data/chunli/coastline/orthorectwork/']; %to store new orthorectifed images
tiledirnew=[macdir,'/data/chunli/coastline/tiledirnew/']; %to store new downloaded dems

