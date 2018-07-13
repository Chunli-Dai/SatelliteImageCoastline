function BW2=icluster(BW,nlb,demgt,cloudth)
% ICLUSTER: 
%   Identify the cluster that contain a given point.
%   Update: return all clusters that are clouds over water.
% 
%   P : the coordinates of given points, in the form of the index of BW. 
%   BW can be a logical or numeric array of any dimension, and it must be
%   nonsparse.
%   nlb: width of the cluster to be eroded.
%   demgt: land cover classification of target strip, %-1,non value;1 water; 0 non-water
%   cloudth: threshold of cloud pixels
%   Output:
%   BW2: the image cluster of BW that contain the given points.
%   Update: only produce the first cloudth pixels for faster speed.

%   See also BWCONNCOMP, CONNDEF, REGIONPROPS.
%   Chunli Dai, April 2018, chunlidai1@gmail.com

    %get rid of a long thin (100 m width) beach band
    Mextra=BW;
    Me1=imerode(Mextra, ones(nlb));   
%   Md1 = imdilate(Me1, ones(width2/resr));
    %remove small clusters
    Modj= bwareaopen(Me1, 10*nlb);
    CCe1 = bwconncomp(Modj);

    %find the original clusters of Mextra that are not filtered during the erode, still appears in Me1 (eroded Mextra).
    CCe1st=zeros(CCe1.NumObjects,1);
    for k=1:CCe1.NumObjects
    CCe1st(k)=CCe1.PixelIdxList{k}(1); % first point of each cluster.
    end
    P=CCe1st;
    CC = bwconncomp(BW);
    np=length(P);
    ide1=zeros(np,1);
    BW2=BW;BW2(:)=0;    
    for k=1:CC.NumObjects
        %find the original clusters of Mextra that are not filtered during the erode
        Lia=ismember(P,CC.PixelIdxList{k}); %if the 1st is found in 2nd.
        idt=find(Lia==1);
        if ~isempty(idt)
            
            %check if this cluster is a long thin (100 m width) beach band 
            BW3=BW;BW3(:)=0;
            BW3(CC.PixelIdxList{k})=BW(CC.PixelIdxList{k});
            Me1=imerode(BW3, ones(nlb));  
            ratio=sum(Me1(:))/sum(BW3(:))*100; % in percentage of how many points left.
            if ratio<10 % either very small or very long shape
                continue
            end
            
            if 1
            % check whether this cluster is clouds or land with occasional water
            Me1=imdilate(BW3, ones(3))-BW3; % the location of immediate surrounding of the extra area. 
            M1e1=demgt(logical(Me1)); % 
            ratio=sum(M1e1==1)./sum(M1e1~=-1)*100.;
            if ratio > 50 % surrounding has some water (more than 50%): means clouds over water
                % do nothing
            else % no water surround the extra area, means it's land area, and the stacked NDWI is lake.
                continue
            end
            end
            
            ide1(idt)=k; %ide1 the id of BW clusters that contain given Points.
            BW2(CC.PixelIdxList{k})=BW(CC.PixelIdxList{k});
            
            % for faster speed
            cloudsum=sum(sum(BW2)); %pixels
            if cloudsum > cloudth 
                return
            end
       end
    end
    
return
end