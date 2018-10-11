function [Co]=CoastTileMonoProb(tilefile,S);
% coastline from orthoimages;mosaicking of coastlines
% Requirements: 
%               
% Versions:
% CoastTilev2.m: to save space: do not save all groups, just save the data that has the maximum repeats, and disgrad others.
%  	July 2018: collect the edge line of each images, and set the a priori coastline buffer zone as edge.
%		   and the edge of overlapping groups.
%		   Crop the area of coastal band to save space.
% 

% Chunli Dai chunlidai1@gmail.com
% March 2018
% 

% close all
% clear 
% clc


%load and set up the parameter
constant

flagplot=0;

res='2';
% demdir=dir(['/data2/ArcticDEM/region_',regionnum,'*']);
macdir='/Users/chunlidai/surge/';
% macdir='/Users/chunli/surge/';
macdir=[];

% addpath(genpath([macdir,'/../Box Sync/ISSMCoastline/']));
Co=[];

if 0 %for testing
addpath(genpath([macdir,'/data/chunli/scripts/']));
addpath(genpath([macdir,'/data/chunli/coastline/']));
addpath(genpath(multidir));

tilefile=[macdir,'/data/chunli/coastline/orthorect1/54_06_2_2_5m_v2.0/54_06_2_2_5m_v2.0_reg_dem.tif'];
tilecoastname='tilegshhs.shp';
S = shaperead(tilecoastname);
cnt=length(S);
end

%resr=2.; %str2double(res)/dsr;
% resr=40.; 
resrc=40.; %for coregisteration
% res=mt.info.map_info.dx;
% demdir=[macdir,'/data2/ArcticDEM/region_08_canada_baffin/tif_results/8m/'];

%Grids of ArcticDEM Tiles
% 54_06_2_2_5m_v2 yid_xid_xids_yids
% 1,2 ; 2,2
% 1,1 ; 2,1
%x=x0+(xid-1)*dx; %left edge of a box;
dx=100e3;x0=-4000e3;y0=-4000e3;%xe=3400e3;ye=4000e3; %ArcticDEM Mosaic tiles coordinate reference;

% %%% Preparation: get rang0 from tile name
[dir,ifile,ext] =fileparts(tilefile);
r=1;
xid= sscanf(ifile(r+3:(r+4)), '%g', 1);
yid= sscanf(ifile(r:(r+1)), '%g', 1);
xids= sscanf(ifile(r+6), '%g', 1);
yids= sscanf(ifile(r+8), '%g', 1);
%     [xid,yid,xids,yids]
%    xid=14;yid=51;xids=2;yids=1;%51_14_2_1_5m_v2.0
x=x0+(xid-1)*dx+(xids-1)*dx/2;y=y0+(yid-1)*dx+(yids-1)*dx/2;
rang0b=[x x+dx/2 y y+dx/2]; %exact tile boundary;for calculation and final results griding.
ranget=round(rang0b/resr)*resr;rang0b=ranget;

rang0=[x-width x+dx/2+width y-width y+dx/2+width]; %tile boundary with buffer width for selecting files and grouping; 
				%need to be the buffer with of a priori incase of boundary being eroded.
ranget=round(rang0/resr)*resr;rang0=ranget;

rang1=[x-widthstat x+dx/2+widthstat y-widthstat y+dx/2+widthstat];ranget=round(rang1/resr)*resr;rang1=ranget;

%54_06_2_2_5m_v2.0_reg_dem.tif  54_06_2_2_coast_v1.0.tif 54_06_2_2_prob_v1.0.tif  54_06_2_2_lessdata_v1.0.tif 
ofile1=['output/',ifile(1:10),'coast',num2str(probthre),'_v1.0.shp'];
probfile=[ifile(1:10),'prob_v1.0.tif'];
[status , cmdout ]=system(['find ',probdir,' -name ',probfile]);
if  ~isempty(cmdout) %
    ofile2=deblank(cmdout);
else
    fprintf([probfile,' not found! \n'])
    return
end
% ofile4=['output/',ifile(1:10),'bound_v1.0.tif']; %image boundaries, no data areas
ofile4=strrep(ofile2,'prob','bound');

% rang0 = [  -3450000    -3400000     1350000     1400000];

xeq=(rang0(1)+rang0(2))/2;yeq=(rang0(3)+rang0(4))/2;
[lateq,loneq]=polarstereo_inv(xeq,yeq,[], [],70,-45);
formatSpec = '%6.1f';

% rang0=[-3467 -3455 110 124 ]*1e3;
x0=[rang1(1) rang1(2) rang1(2) rang1(1) rang1(1) ];y0=[rang1(4) rang1(4) rang1(3) rang1(3) rang1(4) ];
[lat0,lon0]=polarstereo_inv(x0,y0,[], [],70,-45);
tx=rang1(1):resr:rang1(2);ty=rang1(4):-resr:rang1(3);
xout=tx;yout=ty; %needs a buffer zone of the box avoiding excluding small water body at boundaries.

cnt=length(S);

% %%% Preparation: get water mask from a priori coastline
xw=rang0(1):resrc:rang0(2);yw=rang0(4):-resrc:rang0(3); % a priori water mask
nwx=length(xw);
% method 1 %poly2mask, fast 0.5 sec.
tic
smg=false(nwx,nwx);%water mask from a priori coastline shapefiles
for j=1:cnt
    [sx,sy]=polarstereo_fwd(S(j).Y,S(j).X,[], [],70,-45);
    dx=resrc;
    idx=round((sx-rang0(1))/dx)+1;idy=round((sy-rang0(4))/(-dx))+1;
    idt=isnan(idx)|isnan(idy);
    idx(idt)=[];idy(idt)=[];
    sm=poly2mask(idx,idy,nwx,nwx); % fast, apply to each polygon one by one.
    smg=smg|sm;
end
wm=[];wm.x=xw;wm.y=yw;wm.z=smg;clear smg; %1 land, 0 water
if(sum(wm.z(:))==0||sum(wm.z(:))==nwx*nwx);fprintf('This tile contain no coastline (all land or all ocean).');return;end % if all land or all water, i.e. no coastline.
% get the buffer region of coastline: +- width of the coastline;
Md1 = imdilate(wm.z, ones(width/resrc));
Me1=imerode(wm.z, ones(width/resrc)); % erode the box edge also; to fix: change width0 in CoastTileMain to be larger than width;
Mcb=Md1&~Me1; %coastal band


data=readGeotiff(ofile2);
[nsuby,nsubx]=size(data.z);
prob=data.z;xout=data.x;yout=data.y;

if ~exist(ofile4,'file')
fprintf([ofile4,' does not exist. \n'])
Medgs1=false(nsuby,nsubx); 
else
data=readGeotiff(ofile4);
Medgs1=logical(data.z); %1 edge; 0 non edge
end

%Based on probability only, prob=[0, 100]; if 255, it means no data.
jump=-1*ones(nsuby,nsubx,'int8'); %%-1,non value;1 water; 0 non-water
jump(prob>=probthre&prob~=255)=1;
jump(prob<probthre)=0;
Me1o=interp2(wm.x,wm.y,Me1,xout,yout','*nearest',1);%1 land, 0 water
jump(Me1o)=0;   %hi; fill inland with land;%for inland regions far from coast, set it as edges.
Md1o=interp2(wm.x,wm.y,Md1,xout,yout','*nearest',1);
jump(Md1o==0)=1;   %hi; fill the ocean with water
Medgsap=Me1o|~Md1o; %set the a priori buffer zone as edges.
clear Me1o Md1o
%jump(prob==255|novlpf<=1)=-1; %considear novlpf<=2 %no data/edge; Move after fill (Modfil) to mitigate the separation of water bodies (Bug 6).

%get coast from water mask.
%apply the coastline buffer band to the mask
%Medgs=isnan(jump);
Medgs=(jump==-1);%-1,non value;1 water; 0 non-water
Modj=jump;Modj(Medgs)=0;%edges as land,but would be filterd out again later.
Modj= bwareaopen(Modj, lakearea/resr/resr); %remove small clusters; smallest lake that would be kept;
%Around edge of image, add a thin band of water, so that the small area at edge be filled; then later remove it; 
Med=(imdilate(~Medgs,ones(4)))&Medgs;
Modj(Med)=1;%1 is water
%smallest island is 46 meters long by 16 meters (Bishop Rock).
Modfil = bwareaopen(~Modj, cloudarea/resr/resr); %fill small areas: 1e4*4m^2, smallest island that would be kept.
Modfil=~Modfil;
Modfil(Med)=0;%remove the artifical water around edge;
clear Med Modj jump


B = bwboundaries(~Modfil,'noholes'); %^_^ keep the unwanted lake islands.
n=length(B);xo=cell(n,1);yo=cell(n,1);
clear prob 
tic
peri=sqrt(4*pi*cloudarea/resr/resr); % shortest perimeter of a given area.
for k=1:n %can be slow if poor data quanlity, lots of scatterred points;workstereoprob2wocoregcrop/55_06_2_1_coast_v1.0.shp includes island of lakes!
    xid=B{k}(:,2); yid=B{k}(:,1);zid=ones(size(xid));
    if (length(xid)<peri);continue;end %to be immune to bug 11.
    %get the boundaries in terms of mask
    BWb=zeros(size(Modfil));BWb((xid-1)*nsuby+yid)=1;
    Mt=BWb&Medgs1;%filter out edges
    ne=sum(Mt(:));
    if ne~=0 %find the id and replace them with NaN;
       [idy,idx]=find(Mt==1);
       
       for j=1:ne
          Ml=xid==idx(j)&yid==idy(j);
          zid(Ml)=0;%fall into edges
       end
       
    end
    x=xout(xid);y=yout(yid);    
    [LAT,LON]=polarstereo_inv(x,y,[], [],70,-45);
    LAT(zid==0)=nan;LON(zid==0)=nan;%polygon but with only valid polylines displayed.
    zidi=find(zid==0);zidid=zidi(2:end)-zidi(1:end-1);id=find(zidid==1);%delete sequential nans
    idx=zidi(id+1);
    %figure;plot(x,y,'go');hold on;plot(x(idx),y(idx),'r.')
    LON(idx)=[];LAT(idx)=[];
    xo{k}=LON;yo{k}=LAT;
end
fprintf('Retrieve boundary')
toc
clear Mt BWb Medgs1
idd=find(cellfun(@isempty,xo)); %fix bug 12
xo(idd)=[];yo(idd)=[];
% shp1 = struct('Geometry', 'PolyGon', 'X', xo, 'Y', yo);
shp1 = struct('Geometry', 'PolyLine', 'X', xo, 'Y', yo);
%shapewrite(shp1, 'coastline.shp');
if ~isempty(shp1)
shapewrite(shp1, ofile1);
end

if flagplot==1
[X,Y]=meshgrid(xout,yout);
figure;imagesc(xout,yout,double(Modfil),'alphadata',~Medgs)
hold on;plot(X(M),Y(M),'ro')
xl=X(M);yl=Y(M);
if 0
hold on;plot(xl(probl>0.6)*1e-3,yl(probl>0.6)*1e-3,'r.');
hold on;plot(xl(probl<=0.6&probl>0.2)*1e-3,yl(probl<=0.6&probl>0.2)*1e-3,'m.');
hold on;plot(xl(probl<=0.2&probl>=0.1)*1e-3,yl(probl<=0.2&probl>=0.1)*1e-3,'k.');
hold on;plot(xl(probl<0.1)*1e-3,yl(probl<0.1)*1e-3,'c.');
print('-dpdf','-r300','Tilemask') 
% saveas(gcf,'Tilemask','fig')
savefig(gcf,'Tilemask','compact')

figure;hist(probl)
saveas(gcf,'probabilityw','fig')
end

% plot the chosen piece id for each output pixel
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]);
hold all
xt=xout*1e-3;yt=yout*1e-3;zt=iselop; %to check if the order is consistent with idun
imagesc(xt,yt,zt);
colormap jet;colorbar; shading interp; axis equal; view(0,90)
set(gca,'FontSize', 18);
plot(xeq*1e-3, yeq*1e-3,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
axis(rang0*1e-3)
title('Chosen piece isel')
print('-dpdf','-r300','Chosenisel') 
%saveas(gcf,'Chosenisel','fig')
%savefig(gcf,'Chosenisel','compact')

%Final count of repeats for the chosen piece 
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
imagesc(xt,yt,novlpf);
colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
hold on
plot(xeq*1e-3, yeq*1e-3,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
title('Final repeat, exclude clouds')
print('-dpdf','-r300','OverlappingCount_Final') 
%saveas(gcf,'OverlappingCount_Final','fig')
%savefig(gcf,'OverlappingCount_Final','compact')
end %if flagplot==1

end
