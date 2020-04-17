function [] = write_pd(fd,coordinate,type,vol)

na=size(coordinate,1);

for j=1:size(coordinate,3)
    for i=1:na
        fprintf(fd,'%20.12e %20.12e %20.12e %2.0f %20.12e\n',coordinate(i,:,j),type(i),vol(i));
    end
end
