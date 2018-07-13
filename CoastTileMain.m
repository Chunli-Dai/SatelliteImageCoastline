% main program for getting coastline for each ArcticDEM tile
% Requirements: gdal software;ogr2ogr
% %%%% inputs needed
currentdir=pwd;
addpath(genpath(currentdir));

% %%%% control parameters
codedir=['/data/chunli/coastline/codec2/'];
addpath(genpath(codedir));
constant
width0=width+1e3; %buffer width of the a priori coastline, e.g., 2km. need to be slightly larger than width in Coastline.m
% Otherwise, one column of 0 at the edge leads to the unwanted erode in Me1=imerode(wm.z, ones(width/resrc)); (Coastline.m)

%Preparation: building folders
if ~exist(orthworkdir,'dir')
  mkdir(orthworkdir)
end
if ~exist(tiledirnew,'dir')
  mkdir(tiledirnew)
end
if ~exist('output','dir')
  mkdir('output')
end

if 0
%General directory that contains the mosaic tile DEM files, such as /elev/dem/setsm/ArcticDEM/mosaic/v2.0/
tiledir=[macdir,'/data/chunli/coastline/'];%ArcticDEM mosaic tile directory. 
stripdir='/*/ArcticDEM/region*/strips/2m/'; %directory of strip files
multidir=[macdir,'/data1/pgc_projects/coastline/imagery/']; %directory of multispectral images
end

%addpath(genpath(multidir));

%shpname='./GSHHS/GSHHS_f_L1.shp';% a priori coastline shapefile
[status , cmdout ]=system(['find ',codedir,' -name GSHHS_f_L1.shp']);
shpname=deblank(cmdout);

% %%% Preparation: get the list of strip files and boundries
[status , cmdout ]=system(['find ../ -name ','run_overlap_strip_io.sh']);
shname=deblank(cmdout);
filename='boundaries_regall_strip.dat'; %'boundaries_reg31.dat';
if ~exist(filename,'file')
    str=sprintf('%s  ''%s''',shname,deblank(stripdir));
  [status, cmdout]=system(str);
end

fprintf ('\n Step 0: geting the boundary for all files in the region.')
%READ INPUT PARAMETERS; getting the boundaries for all files
% filename='boundaries_reg31.dat';
fid = fopen(filename);
n = linecount(fid);
fid = fopen(filename);
% range=fscanf(fid, '%f', [4, n))';
range=zeros(n,4);
%exclude panchromatic bands; may be included later.
idd=[];
for i=1:n
   range(i,1:4)=fscanf(fid, '%f', [4, 1])';ifile=fgetl(fid);
   [demdir,name,ext] =fileparts([macdir,strtrim(ifile)]);
   f{i}=[name,ext];
   fdir{i}=[demdir,'/'];%working on two regions 
    satname=f{i}(1:4);
    if strcmp(satname,'WV01')||strcmp(satname,'GE01')
        idd=[idd;i];
    end
end
range(idd,:)=[];f(idd)=[];fdir(idd)=[];
display(['demdir=',demdir])


% ArcticDEM mosaic tile grids
if 0
xidg=5:8;
yidg=53:55;
[X,Y]=meshgrid(xidg,yidg);
[n,m]=size(X);
end % if 0 

filename='tilelist';
fid = fopen(filename);
n = linecount(fid);
fid = fopen(filename);
ftile=cell(n,1);
for i=1:n
   ftile{i}=fgetl(fid);
end

parpool(2)
%parfor xyid=1:n*m %5:8 %6:7  %1:80
parfor xyid=1:n
%   for yid=53:55 %1:80
%       for xids=1:2
%           for yids=1:2
                % 54_06_2_2_5m_v2 yid_xid_xids_yids %name convention
%   xid=7;yid=54;xids=1;yids=1;
%   xid=X(xyid);yid=Y(xyid);
%   tilefile=sprintf('%02d_%02d_%01d_%01d_5m_v2.0_reg_dem.tif',yid,xid,xids,yids);  %'54_06_2_2_5m_v2.0_reg_dem.tif';
    tilefile=ftile{xyid};
    [dir,ifile,ext] =fileparts(tilefile);
    r=1;
    xid= sscanf(ifile(r+3:(r+4)), '%g', 1);
    yid= sscanf(ifile(r:(r+1)), '%g', 1);
    xids= sscanf(ifile(r+6), '%g', 1);
    yids= sscanf(ifile(r+8), '%g', 1);

    % get the data boundary, rang0, of this DEM tile 
    dx=100e3;x0=-4000e3;y0=-4000e3;%xe=3400e3;ye=4000e3; %ArcticDEM Mosaic tiles coordinate reference;
    x=x0+(xid-1)*dx+(xids-1)*dx/2;y=y0+(yid-1)*dx+(yids-1)*dx/2;
%     rang0=[x x+dx/2 y y+dx/2];
    
    % Find whether this tile contains any coastline.
    % Buffer the tile boundary by width;
    rang0=[x-width0 x+dx/2+width0 y-width0 y+dx/2+width0];
    x0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];y0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];
    [lat0,lon0]=polarstereo_inv(x0,y0,[], [],70,-45);
    bb = geoshape(lat0,lon0,'Geometry','Polygon');
    tileshape=sprintf('output/%02d_%02d_%01d_%01d_tile.shp',yid,xid,xids,yids); 
    tilecoastname=sprintf('output/%02d_%02d_%01d_%01d_tilegshhs.shp',yid,xid,xids,yids);
    shapewrite(bb,tileshape);
    [status, cmdout]=system(['rm ',tilecoastname]);
    system(['time ogr2ogr -overwrite -clipsrc ',tileshape,' ',tilecoastname,' ',shpname]);
    %ogr2ogr -overwrite -clipsrc tile.shp tilegshhs.shp GSHHS/GSHHS_f_L1.shp
    S = shaperead(tilecoastname);
    cnt=length(S); %figure;mapshow(S); %InLand have one boundary box, i.e. cnt=1. For Ocean, cnt=0;
    if cnt==0; continue;end
    
    % find the dem tile file or download the data
    [status , cmdout ]=system(['find ',tiledir,' -name ',tilefile]); %status always 0, cmdout could be empty.
    if ~isempty(cmdout) && status ==0 % 
        tilefile=deblank(cmdout);
    else
        warning(['Tile file ',tilefile,' not found! Download it from website'])
        [dir,name,ext] =fileparts(tilefile);
    %     http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v2.0/43_59/43_59_2_1_5m_v2.0.tar
        tarfile=[name(1:17),'.tar'];
%       webtile=deblank(['http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v2.0/',name(1:5),'/',tarfile,'.gz']);
        webtile=deblank(['http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v2.0/',name(1:5),'/',tarfile]);
        system(['wget  ',webtile,' -o downloadlog'])
        system(['tar -xvf ',tarfile]);
	%collect all downloaded mosaic dems to tiledirnew/ % to do
	system(['mv *dem_meta.txt  *reg.txt *.tar *_reg_dem.tif *_reg_matchtag.tif ',tiledirnew]);
        [status , cmdout ]=system(['find .. -name ',tilefile]);
        tilefile=deblank(cmdout);
	if ~exist(tilefile,'file')
	  continue
	end
    end
fprintf (['\n Working on tile:',tilefile,'; \n'])

            [Co]=CoastTile(tilefile, S,range,f,fdir);
    tileshape=sprintf('output/%02d_%02d_%01d_%01d_tile*',yid,xid,xids,yids); 
	     system(['rm ',tileshape]);
%           end %yids
%       end%xids
%   end % yid
end %xid

