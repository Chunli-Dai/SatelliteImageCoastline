function [XYbi,rangei]=imagebd(ifile);
%get the Polygon boundary, and the rectangle range of a mono image from .xml file
% e.g. /data1/pgc_projects/dai_aleutians_multi_mono/imagery//WV02/WV02_20130923221601_1030010026BD6B00_13SEP23221601-M1BS-500127380110_01_P001.xml 

	vstr={'<ULLON>','<ULLAT>','<URLON>','<URLAT>','<LRLON>','<LRLAT>','<LLLON>','<LLLAT>'};

	metafile=ifile;
	c=textread(metafile,'%s','delimiter','\n');

% Get the Footprint Vertices X, Y, close the loop
	n=length(vstr)/2;
	lon=zeros(n+1,1);lat=zeros(n+1,1);
	for i=1:n*2
	str=vstr(i);
	r=find(~cellfun(@isempty,strfind(c,str)));
	if isempty(r)
		warning(['xml file is different as anticipated.',ifile])
        XYbi=[0 0]; 	rangei=[0 0 0 0];
		return;
	end %
	c2=c{r(1)};
	r1=strfind(c2,'>');r2=strfind(c2,'</');
    c2([1:r1(1),r2(1):end])='';
	z = sscanf(c2, '%g', 1);
	j=ceil(i/2);
	if mod(i,2)  %1 odd number, 0 even
           lon(j)=z;
	else
	   lat(j)=z;
	end
	end % if i
	lon(n+1)=lon(1);lat(n+1)=lat(1); %close the loop
	[Xb,Yb]=polarstereo_fwd(lat,lon,[], [],70,-45);

        XYbi=[Xb,Yb]; %n by 2
	rangei=[min(Xb) max(Xb) min(Yb) max(Yb)];

return
end
