clc
clear
xmin=0
xmax=12
ymin=0
ymax=8
dx=4
dy=4
kx=0.000001
ky=0.000001
hmax=13
hmin=8
 nodesX=((xmax-xmin)/dx+1);
 nodesY=((ymax-ymin)/dy+1);
 nodes=nodesX*nodesY;
 elements=(nodesX-1)*(nodesY-1)*2;
abc=zeros(nodes,3);
b=zeros(elements,3);
elementX=(nodesX-1)*2;
elementY=(nodesY-1)*2;
offset=zeros(elements,1);
type=zeros(elements,1);
kelem=zeros(3,3);
grauslibertade=zeros(nodes,elements);
kglobal=zeros(nodes,nodes);
hmed=(hmax+hmin)/2;
h=zeros(nodes,1);
Q=zeros(nodes,1);
Ad=zeros(elements,1);
velocity=zeros(elements,2);
vcoord=zeros(elements,2);
vresult=zeros(elements,1);
for i=1:elements
    offset(i)=3*i;
    type(i)=5;
end 

for k=1:nodesY    
    for l=1:nodesX    
        i= nodesX*(k-1)+l;
        for j=1:2
            if j==1
            abc(i,j)=xmin+(l-1)*dx;            
            else
            abc(i,j)=ymin+(k-1)*dy;
            end
        end
    end
end
for k=1:nodesY-1
    
    for l=1:nodesX-1        
        for m=0:1        
            i=2*(nodesX-1)*(k-1)+(2*l-1+m);            
                for j=1        
                 b(i,j)=(i+m+1)/2+(k-1);        
                end                 
                for j=2
                    if m==0
                        b(i,j)=(i+m+1)/2+(k);                       
                    else                        
                         b(i,j)=nodesX*(k)+l+1;
                        end                     
                end                 
                for j=3                    
                    b(i,j)=nodesX*(k)+l; 
                end
        end
    end
end




%abc

connectivity=b-1
%offset
%type
%nodes
%elements
%b

for j=1:elements
    for k=1:3
        grauslibertade(b(j,k),j)=k;
    end
end
grauslibertade

for i=1:elements
    Ad(i,1)=abc(b(i,2),1)*abc(b(i,3),2)-abc(b(i,3),1)*abc(b(i,2),2)-abc(b(i,1),1)*abc(b(i,3),2)+abc(b(i,1),1)*abc(b(i,2),2)+abc(b(i,3),1)*abc(b(i,1),2)-abc(b(i,2),1)*abc(b(i,1),2);
end

%Ad

for k=1:elements
  
            kelem(1,1)=((abc(b(k,2),2)-abc(b(k,3),2))^2*kx+(abc(b(k,3),1)-abc(b(k,2),1))^2*ky)/(2*Ad(k,1));
            kelem(1,2)=((abc(b(k,2),2)-abc(b(k,3),2))*(abc(b(k,3),2)-abc(b(k,1),2))*kx+(abc(b(k,3),1)-abc(b(k,2),1))*(abc(b(k,1),1)-abc(b(k,3),1))*ky)/(2*Ad(k,1));
            kelem(1,3)=((abc(b(k,2),2)-abc(b(k,3),2))*(abc(b(k,1),2)-abc(b(k,2),2))*kx+(abc(b(k,3),1)-abc(b(k,2),1))*(abc(b(k,2),1)-abc(b(k,1),1))*ky)/(2*Ad(k,1));
            kelem(2,1)=((abc(b(k,3),2)-abc(b(k,1),2))*(abc(b(k,2),2)-abc(b(k,3),2))*kx+(abc(b(k,1),1)-abc(b(k,3),1))*(abc(b(k,3),1)-abc(b(k,2),1))*ky)/(2*Ad(k,1));
            kelem(2,2)=((abc(b(k,3),2)-abc(b(k,1),2))^2*kx+(abc(b(k,1),1)-abc(b(k,3),1))^2*ky)/(2*Ad(k,1));
            kelem(2,3)=((abc(b(k,3),2)-abc(b(k,1),2))*(abc(b(k,1),2)-abc(b(k,2),2))*kx+(abc(b(k,1),1)-abc(b(k,3),1))*(abc(b(k,2),1)-abc(b(k,1),1))*ky)/(2*Ad(k,1));
            kelem(3,1)=((abc(b(k,1),2)-abc(b(k,2),2))*(abc(b(k,2),2)-abc(b(k,3),2))*kx+(abc(b(k,2),1)-abc(b(k,1),1))*(abc(b(k,3),1)-abc(b(k,2),1))*ky)/(2*Ad(k,1));
            kelem(3,2)=((abc(b(k,1),2)-abc(b(k,2),2))*(abc(b(k,3),2)-abc(b(k,1),2))*kx+(abc(b(k,2),1)-abc(b(k,1),1))*(abc(b(k,1),1)-abc(b(k,3),1))*ky)/(2*Ad(k,1));
            kelem(3,3)=((abc(b(k,1),2)-abc(b(k,2),2))^2*kx+(abc(b(k,2),1)-abc(b(k,1),1))^2*ky)/(2*Ad(k,1));
            
            for l=1:nodes
                for m=1:nodes
                
                    if grauslibertade(l,k)== 0   
                        kglobal(l,m)=kglobal(l,m);
                    elseif grauslibertade(m,k)== 0
                        kglobal(l,m)=kglobal(l,m);
                    else
                      kglobal(l,m)=kglobal(l,m)+kelem(grauslibertade(l,k),grauslibertade(m,k));
                    end
                end
            end
end




for i=1:nodes
    if abc(i,1)==xmax
        h(i,1)=hmed;
    elseif abc(i,1)==xmin
        h(i,1)=hmax;
    elseif abc(i,2)==ymax
       h(i,1)=hmax;
    else
        h(i,1)=0;
    end
end 

for i=1:nodes
    Q(i,1)=h(i,1);
end
Q



for i=1:nodes
   
    if Q(i,1)~=0  
        for k=1:nodes
        kglobal(i,k)=0;
        end
        kglobal(i,i)=1e0;
    end
    
    
end


for i=1:elements

        vcoord(i,1)=(abc(b(i,1),1)+abc(b(i,2),1)+abc(b(i,3),1))/3;
        vcoord(i,2)=(abc(b(i,1),2)+abc(b(i,2),2)+abc(b(i,3),2))/3;

end


    
kglobal
hp=zeros(nodes,1)
u=zeros(nodes,1)

S=kglobal\Q
for i=1:nodes
hp(i,1)=S(i,1)-abc(i,2);
end
u=hp*10 %kPa



for i=1:elements
        velocity(i,1)=kx*(S(b(i,1),1)*(abc(b(i,2),2)-abc(b(i,3),2))+S(b(i,2),1)*(abc(b(i,3),2)-abc(b(i,1),2))+S(b(i,3),1)*(abc(b(i,1),2)-abc(b(i,2),2)))/(Ad(i,1));
        velocity(i,2)=ky*(S(b(i,1),1)*(abc(b(i,3),1)-abc(b(i,2),1))+S(b(i,2),1)*(abc(b(i,1),1)-abc(b(i,3),1))+S(b(i,3),1)*(abc(b(i,2),1)-abc(b(i,1),1)))/(Ad(i,1));
end

for i=1:elements
        vresult(i,1)=(velocity(i,1)^2+velocity(i,2)^2)^0.5;
end

vcoord
velocity
vresult





nodes
elements
fid=fopen('MyFile.vtu','w');
fprintf(fid, '<?xml version="1.0"?>\n');
fprintf(fid, '<VTKFile type="UnstructuredGrid\" version=\"0.1\" byte_order="LittleEndian">\n');
fprintf(fid, '<UnstructuredGrid>\n');


fprintf(fid, '<Piece NumberOfPoints=\"');
fprintf(fid, '%d',nodes);
fprintf(fid, '\" NumberOfCells=\"');
fprintf(fid, '%d',elements);
fprintf(fid, '\"> \n');
fprintf(fid, '<Points> \n');
fprintf(fid, '<DataArray type=\"Float32\" NumberOfComponents=\"3" Format=\"ascii\">\n');
fprintf(fid, '%d %d %d \n',[abc].');
fprintf(fid, '</DataArray>\n');
fprintf(fid, '</Points>\n');
fprintf(fid, '<Cells>\n');
fprintf(fid, '<DataArray type="Int32" Name="connectivity" Format="ascii">\n');
fprintf(fid, '%d %d %d \n',[connectivity].');
fprintf(fid, '</DataArray>\n');
fprintf(fid, '<DataArray type="Int32" Name="offsets" Format="ascii">\n');
fprintf(fid, '%d \n',[offset]);
fprintf(fid, '</DataArray>\n');
fprintf(fid, '<DataArray type="UInt8" Name="types" Format="ascii">\n');
fprintf(fid, '%d \n',[type]);
fprintf(fid, '</DataArray>\n');
fprintf(fid, '</Cells>\n');
fprintf(fid, '<PointData Scalars="scalars">\n');
fprintf(fid, '<DataArray type="Float32" Name="hh" Format="ascii">\n');
fprintf(fid, '%d \n',[S]);
fprintf(fid, '</DataArray>\n');
fprintf(fid, '<DataArray type="Float32" Name="u" Format="ascii">\n');
fprintf(fid, '%d \n',[u]);
fprintf(fid, '</DataArray>\n');
fprintf(fid, '<DataArray type="Float32" Name="v" Format="ascii">\n');
fprintf(fid, '%d \n',[vresult]);
fprintf(fid, '</DataArray>\n');
fprintf(fid, '</PointData>\n');
fprintf(fid, '</Piece>\n');
fprintf(fid, '</UnstructuredGrid>\n');
fprintf(fid, '</VTKFile>\n');



fclose(fid);
