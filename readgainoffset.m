function [GAIN,OFFSET,Esun]=readgainoffset(satname) %read Gain OFFset data, GainOffset.txt.
%readgainoffset.m
%
if strcmp(satname,'GE01') || strcmp(satname,'QB02') ||strcmp(satname,'IK01')
  iBlue=1;iGreen=2;iRed=3;iNIR1=4;iPan=5;
  ig=[iBlue,iGreen,iRed,iNIR1,iPan];
  GAIN=zeros(5,1);OFFSET=zeros(5,1);Esun=zeros(5,1);
elseif strcmp(satname,'WV02') ||strcmp(satname,'WV03') 
% refers to Maglione 2014
%Coastal 0.400 - 0.450; Blue 0.450 - 0.510; Green 0.510 - 0.580; Yellow 0.585 - 0.625
%Red 0.630 - 0.690; Red-Edge 0.705 - 0.745; NIR1 0.770 - 0.895; NIR2 0.860 - 0.940
  iNIR2=8;iCoastal=1;iRed=5;iGreen=3; iNIR1=7;
    iBlue=2;iYellow=4;iRE=6;iPan=9;
    ig=[iCoastal,iGreen,iRed,iNIR2,iNIR1,iBlue,iYellow,iRE,iPan];

  GAIN=zeros(9,1);OFFSET=zeros(9,1);Esun=zeros(9,1);
elseif strcmp(satname,'WV01')
    ig=[1];
  GAIN=zeros(1,1);OFFSET=zeros(1,1);Esun=zeros(1,1);
end

%read string refers to coastlineArcticDEMv2.m and ChangedeTolv42msv.m
infile='GainOffset.txt';
c=textread(infile,'%s','delimiter','\n');
r=find(~cellfun(@isempty,strfind(c,satname)));
for i=ig
str=c{r+i};
r1=strfind(str,':');
str(1:r1)='';
t1= sscanf(str, '%g %g %g', 3);
GAIN(i)=t1(1);OFFSET(i)=t1(2);Esun(i)=t1(3);
end
end
