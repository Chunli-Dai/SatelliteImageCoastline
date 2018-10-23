function [idx,idn]=separatelines(M,slmt)
%separate lines by nans. 
%Input: M, a vector that has nans, which can be used to separate M into multi segments.
%       slmt, the minimum number of points of a segment.
%Output: idx, the kept parts of M that have segment length > 5.
%        idn, index of M that is one element before each segment.

%debug at CoastTileMono.m
if 0 %test
idb=find(zidid~=1);
M=zidi;M(idb+1)=nan;
end

idx = any(isnan(M),2);
idy = 1+cumsum(idx);
idz = 1:size(M,1);
C = accumarray(idy(~idx),idz(~idx),[],@(r){M(r,:)});

[nc,~]=size(C);

idx=[];idn=[];
for j=1:nc
    [ncj,~]=size(C{j});
    if ncj <=slmt ; continue;end  %if a segment is less than 5 points, ignore it.
    idx=[idx;C{j}(:);];
    idni=min(find(idy==j)); %index of M that is nan and the segment after which is kept.
    idn=[idn;idni];
end
% hold all;plot(x(idx),y(idx),'m<')
end