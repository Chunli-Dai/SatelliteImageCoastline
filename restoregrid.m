function [un,idun]=restoregrid(un,idun,dx,idg,novlp) 
%Fix bug 16
%The problem with grouping: when images are too many, the size of each unique group is too small, -> delete too many pixels that actually have lots of data. 
%If change almt=0, then the computation time is too long, due to the large number of groups (calculation of clouds detection for each group).
%Solution: restore the deleted small groups.

almt2=2e3*2e3; % or 3e3*3e3 %grouping based on regular grids size of 2 km by 2 km.
nx2=round(sqrt(almt2)/dx);
M=double(idun==0);
Modj=bwareaopen(M, nx2*nx2); %ignore pixel groups that have area < nx2*nx2 *dx*dx

if sum(Modj(:))==0; return;end

% figure;plot(X(Modj)*1e-3,Y(Modj)*1e-3,'o-')
[ny,nx]=size(Modj);
ix=1:nx;iy=1:ny;
[IX,IY]=meshgrid(ix,iy);

ix1=min(IX(Modj));ix2=max(IX(Modj));
iy1=min(IY(Modj));iy2=max(IY(Modj));


% within the almt2 area, collect all image ids that fall within this area.
% -> Not gonna work, since the main program will do the overlapping area.
% select the id groups that 1) contain the image that has maximum appearance during all pixels; 
% 2) contain the image that occurs most in the rest of pixels  
% 3) contain the most repeats.
UL=zeros(size(Modj)); %Upper Left corner of each grid
SP=zeros(size(Modj)); %Selected pixel of each grid

idgmd=reshape(idg,ny,nx);

for ix=1:nx
    for iy=1:ny
    idgm{iy,ix}=num2str(idgmd{iy,ix}(:)');
    end
end

npc=length(un);
ipc=npc;
for ix=ix1:nx2:ix2 %row
    jj=ix:min((ix+nx2-1),ix2);
for iy= iy1:nx2:iy2 %column
    UL(iy,ix)=1;
    kk=iy:min((iy+nx2-1),iy2);
    Mj=zeros(size(Modj),'logical');idgmi=repmat({''},ny,nx);
    novlpj=novlp;
    Mj(kk,jj)=Modj(kk,jj);
    idgmi(kk,jj)=idgm(kk,jj);
    idgmi(~Mj)={''}; %mask out area
    novlpj(~Mj)=0;
    
    nt=sum(Mj(:));
    if (nt==0) continue;end
    
    idgmis=idgm(Mj);
    
    idc=[];
    for i1=1:nt
        idc=[idc;str2num(idgmis{i1})'];
    end
    [id1]=mode(idc);
    id1s=num2str(id1);
    
    %1 find the image that has maximum appearance during all pixels; 
%     idgmil=reshape(idgmi,ny*nx,1);
    r=find(~cellfun(@isempty,strfind(idgmi,id1s)));
%   r=find(contains(idgmi,id1s));
    M1=false(ny,nx);
    M1(r)=1;
    
    %2 the image that occurs most in the rest of pixels  
    M1r=Mj&~M1;
    
    if sum(M1r(:))~=0
    idgmis=idgm(M1r);
    idc=[];
    for i1=1:nt
        idc=[idc;str2num(idgmis{i1})'];
    end
    [id1]=mode(idc);
    id1s=num2str(id1);
    r=find(contains(idgmi,id1s));
    M2=false(ny,nx);
    M2(r)=1;
    else
        M2=false(ny,nx);
    end
    Mc=M1&M2;

    if sum(Mc(:))~=0 % both 1 and 2 satisfies;
        novlpj(~Mc)=0; %mask out pixels not satisfying both 1 and 2.
    end
    %search the maximum repeats
     [~,imr]=max(novlpj(:));

     imrc=find(Mj~=0);
    if ~ismember(imr,imrc)
        imr=imrc(1);
        fprintf(['Restore pixels problem: peak nov is not inside the grid!'])
    end
        
    %select idg
    SP(imr(1))=1;
    
    ipc=ipc+1;
    un{ipc}=num2str(idg{imr(1)}(:)');
    idun(logical(Mj))=ipc;
end
end

if 0
xo=X(Modj);yo=Y(Modj);
    figure;
    hold all;
    plot(xo*1e-3,yo*1e-3,'b.')
% hold all;plot(xo(1:3:end),yo(1:3:end),'k>-') %if simply every 3 points
plot(X(UL==1)*1e-3,Y(UL==1)*1e-3,'gs') %corner left upper
plot(X(SP==1)*1e-3,Y(SP==1)*1e-3,'ro') %%selected 

end

return
end

