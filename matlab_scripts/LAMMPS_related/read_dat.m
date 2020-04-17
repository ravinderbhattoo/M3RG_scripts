function [box,info_atoms,coors,id,label]=read_dat(fd,flag_sort,varargin)

% [box,info_atoms,coordinate,id,label]=read_dat(fd,flag_sort)
% 11/05/2015, Bu, fix skipping the pair section
% 10/30/2015, Bu, fix several bugs
% 10/23/2015, Bu, correct bug in reading tilts
% 10/10/2015, Bu, overhaul

title=fgetl(fd);
fgetl(fd);
na=cell2mat(textscan(fgetl(fd),'%uatoms'));
% fgetl(fd);
ntypes=cell2mat(textscan(fgetl(fd),'%uatom types'));
fgetl(fd);

box=zeros(3,3);
temp=textscan(fgetl(fd),'%f%fxlo xhi');
box(1,1,:)=temp{2}-temp{1};
temp=textscan(fgetl(fd),'%f%fylo yhi');
box(2,2,:)=temp{2}-temp{1};
temp=textscan(fgetl(fd),'%f%fzlo zhi');
box(3,3,:)=temp{2}-temp{1};
temps=fgetl(fd);
if ~isempty(temps)
    temp=textscan(temps,['%f%f%fxy xz yz']);
    box(2,1,:)=temp{1};
    box(3,1,:)=temp{2};
    box(3,2,:)=temp{3};
end

fgetl(fd);
fgetl(fd);
mass=reshape(cell2mat(textscan(fd,'%f%f')),[],2);

temps=fgetl(fd);
if length(temps)>4 && strcmp(temps(1:4),'Pair')
    for i=1:ntypes+2
        fgetl(fd);
    end
    temps=fgetl(fd);
end
atomstr=temp;
fgetl(fd);

frame=cell2mat(textscan(fd,'%f%f%f%f%f%f%f%f%f'));
id=frame(:,1);
label=frame(:,2);
charge=frame(:,3);
coors=frame(:,4:6);

info_atoms=zeros(ntypes,3);
info_atoms(:,1:2)=mass;
for i=1:ntypes
    index=find(label==i,1,'first');
    if ~isempty(index)
        info_atoms(i,3)=charge(index);
    end
end

if flag_sort==1
    [~,I]=sort(id);
    for i=1:size(id,2)
        coors(:,:,i)=coors(I(:,i),:,i);
        label(:,i)=label(I(:,i),i);
        id(:,i)=id(I(:,i),i);
    end
    %             if diff(id,1,2)~=0
    %                 error('check sorting!')
    %             end
end

