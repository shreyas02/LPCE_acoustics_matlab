%% gmshReader
%  This function reads the *.msh file generated from gmsh and creates *.crd
%  (coordinates) file, *.cnn (connectivity file) and *.nbc (nodal file) for
%  boundaries in both 2D and 3D.
%
clc;
clear all;

file = './mesh_v2.msh';
wrkDir = './' ;
problemString = 'cylinderbar' ;
elemType = '4Quad' ;

fileId = fopen(file);

while ~feof(fileId)
    tmp = fscanf(fileId, '%s',1);
    if strcmp(tmp,'$PhysicalNames');
        fgets(fileId);
        nPhysicalNames = fscanf(fileId,'%d',1); fgets(fileId);
        for i=1:nPhysicalNames
            pnId1(i) = fscanf(fileId,'%d',1);
            pnId2(i) = fscanf(fileId,'%d',1);
            pnStr = fscanf(fileId,'%s',1);
            if (i==1)
                B1 = strrep(pnStr,'"','');
            elseif (i==2)
                B2 = strrep(pnStr,'"','');
            elseif (i==3)
                B3 = strrep(pnStr,'"','');
            elseif (i==4)
                B4 = strrep(pnStr,'"','');
            elseif (i==5)
                B5 = strrep(pnStr,'"','');
            elseif (i==6)
                B6 = strrep(pnStr,'"','');
            elseif (i==7)
                B7 = strrep(pnStr,'"','');
            elseif (i==8)
                B8 = strrep(pnStr,'"','');
            elseif (i==9)
                B9 = strrep(pnStr,'"','');
            elseif (i==10)
                B10 = strrep(pnStr,'"','');
            else
                warning('Too many boundaries!!');
            end
        end

    elseif strcmp(tmp,'$Nodes');
        fgets(fileId);
        nNodes = fscanf(fileId,'%d',1); fgets(fileId);
        tmp = fscanf(fileId,'%d %e %e %e',[4 nNodes]);
        nodeId = tmp(1,:)';
        Crds = tmp(2:4,:)';
    elseif strcmp(tmp,'$Elements');
        fgets(fileId);
        nEl = fscanf(fileId,'%d',1); fgets(fileId);
        
        for i=1:nEl
            Elem.Id(i) = fscanf(fileId,'%d',1);
            Elem.Type(i) = fscanf(fileId,'%d',1);
            pT = fscanf(fileId,'%d',1);
            for j=1:pT
                Elem.Tag(i,j) = fscanf(fileId,'%d',1);
            end
            if (Elem.Type(i) == 1)
                conn2Line(i,:) = fscanf(fileId,'%d %d',[1 2]);
            elseif (Elem.Type(i) == 2)
                conn3Tri(i,:) = fscanf(fileId,'%d %d %d',[1 3]);
            elseif (Elem.Type(i) == 3)
                conn4Quad(i,:) = fscanf(fileId,'%d %d %d %d',[1 4]);
            elseif (Elem.Type(i) == 4)
                conn4Tet(i,:) = fscanf(fileId,'%d %d %d %d',[1 4]);
            elseif (Elem.Type(i) == 5)
                conn8Hex(i,:) = fscanf(fileId,'%d %d %d %d %d %d %d %d'...
                    ,[1 8]);
            elseif (Elem.Type(i) == 6)
                conn6Prism(i,:) = fscanf(fileId,'%d %d %d %d %d %d'...
                    ,[1 6]);
            elseif (Elem.Type(i) == 7)
                conn5Pyr(i,:) = fscanf(fileId,'%d %d %d %d %d',[1 5]);
            end
        end
    elseif strcmp(tmp,'');
        
    end
    
end

fclose(fileId);

% Write *.crd file
% filename = sprintf('%s/%s.crd',wrkDir,problemString);
% fileId = fopen(filename,'w');
% for i=1:nNodes
%     fprintf(fileId,'%d %e %e %e\n',nodeId(i),Crds(i,1),Crds(i,2),Crds(i,3));
% end
% fclose(fileId);
crd = [nodeId Crds] ;

% Write *.cnn file
% filename = sprintf('%s/%s.fluid.cnn',wrkDir,problemString);
% fileId = fopen(filename,'w');
% ElemFound = 0 ;
% a = 1;
% for i=1:nEl
% if (strcmp(elemType,'3Tri') && conn3Tri(i,1) ~= 0)
%     ElemFound = 1;
%         fprintf(fileId,'%d %d %d %d\n',a,conn3Tri(i,1),conn3Tri(i,2),...
%             conn3Tri(i,3)) ;
%         a = a+1 ;
% elseif (strcmp(elemType,'4Quad') && conn4Quad(i,1) ~= 0)
%     ElemFound = 1;
%         fprintf(fileId,'%d %d %d %d %d\n',a,conn4Quad(i,1),conn4Quad(i,2),...
%             conn4Quad(i,3),conn4Quad(i,4)) ;
%         a = a+1 ;
% elseif (strcmp(elemType,'4Tet') && conn4Tet(i,1) ~= 0)
%     ElemFound = 1;
%         fprintf(fileId,'%d %d %d %d %d\n',a,conn4Tet(i,1),conn4Tet(i,2),...
%             conn4Tet(i,3),conn4Tet(i,4)) ;
%         a = a+1 ;
% elseif (strcmp(elemType,'8Hex') && conn8Hex(i,1) ~= 0)
%     ElemFound = 1;
%         fprintf(fileId,'%d %d %d %d %d %d %d %d %d\n',a,conn8Hex(i,1),...
%             conn8Hex(i,2),conn8Hex(i,3),conn8Hex(i,4),conn8Hex(i,5),...
%             conn8Hex(i,6),conn8Hex(i,7),conn8Hex(i,8)) ;
%         a = a+1 ;
% elseif (strcmp(elemType,'6Prism') && conn6Prism(i,1) ~= 0)
%     ElemFound = 1;
%         fprintf(fileId,'%d %d %d %d %d %d %d\n',a,conn6Prism(i,1),...
%             conn6Prism(i,2),conn6Prism(i,3),conn6Prism(i,4),...
%             conn6Prism(i,5),conn6Prism(i,6)) ;
%         a = a+1 ;
% elseif (strcmp(elemType,'5Pyr') && conn5Pyr(i,1) ~= 0)
%     ElemFound = 1;
%         fprintf(fileId,'%d %d %d %d %d %d\n',a,conn5Pyr(i,1),...
%             conn5Pyr(i,2),conn5Pyr(i,3),conn5Pyr(i,4),...
%             conn5Pyr(i,5)) ;
%         a = a+1 ;
% end
% end
% fclose(fileId);
if (strcmp(elemType,'3Tri'))
    conn = conn3Tri ;
    conn(conn==0) = [];
    conn = reshape(conn,[],3);
elseif (strcmp(elemType,'4Quad'))
    conn = conn4Quad ;
    conn(conn==0) = [];
    conn = reshape(conn,[],4);
elseif (strcmp(elemType,'4Tet'))
    conn = conn4Tet ;
    conn(conn==0) = [];
    conn = reshape(conn,[],4);
elseif (strcmp(elemType,'8Hex'))
    conn = conn8Hex ;
    conn(conn==0) = [];
    conn = reshape(conn,[],8);
elseif (strcmp(elemType,'6Prism'))
    conn = conn6Prism ;
    conn(conn==0) = [];
    conn = reshape(conn,[],6);
end

% Write *.nbc files
BcNodes1 = [];
BcNodes2 = [];
BcNodes3 = [];
BcNodes4 = [];
BcNodes5 = [];
BcNodes6 = [];
BcNodes7 = [];
BcNodes8 = [];
BcNodes9 = [];
BcNodes10 = [];
for i=1:nEl
    if (Elem.Type(i) == 1)
        if (Elem.Tag(i,1) == pnId2(1))
            BcNodes1 = [BcNodes1; conn2Line(i,:)];
            nBC1 = size(unique(BcNodes1),1) ;
        elseif (Elem.Tag(i,1) == pnId2(2))
            BcNodes2 = [BcNodes2; conn2Line(i,:)];
            nBC2 = size(unique(BcNodes2),1) ;
        elseif (Elem.Tag(i,1) == pnId2(3))
            BcNodes3 = [BcNodes3; conn2Line(i,:)];
            nBC3 = size(unique(BcNodes3),1) ;
        elseif (Elem.Tag(i,1) == pnId2(4))
            BcNodes4 = [BcNodes4; conn2Line(i,:)];
            nBC4 = size(unique(BcNodes4),1) ;
        elseif (Elem.Tag(i,1) == pnId2(5))
            BcNodes5 = [BcNodes5; conn2Line(i,:)];
            nBC5 = size(unique(BcNodes5),1) ;
        elseif (Elem.Tag(i,1) == pnId2(6))
            BcNodes6 = [BcNodes6; conn2Line(i,:)];
            nBC(i) = size(unique(BcNodes6),1) ;
        elseif (Elem.Tag(i,1) == pnId2(7))
            BcNodes7 = [BcNodes7; conn2Line(i,:)];
            nBC(i) = size(unique(BcNodes7),1) ;
        elseif (Elem.Tag(i,1) == pnId2(8))
            BcNodes8 = [BcNodes8; conn2Line(i,:)];
            nBC(i) = size(unique(BcNodes8),1) ;
        elseif (Elem.Tag(i,1) == pnId2(9))
            BcNodes9 = [BcNodes9; conn2Line(i,:)];
            nBC(i) = size(unique(BcNodes9),1) ;
        elseif (Elem.Tag(i,1) == pnId2(10))
            BcNodes10 = [BcNodes10; conn2Line(i,:)];
            nBC(i) = size(unique(BcNodes10),1) ;
        end 
    elseif (Elem.Type(i) == 2)
        if (Elem.Tag(i,1) == pnId2(1))
            BcNodes1 = [BcNodes1; conn3Tri(i,:)];
            nBC1 = size(unique(BcNodes1),1) ;
        elseif (Elem.Tag(i,1) == pnId2(2))
            BcNodes2 = [BcNodes2; conn3Tri(i,:)];
            nBC2 = size(unique(BcNodes2),1) ;
        elseif (Elem.Tag(i,1) == pnId2(3))
            BcNodes3 = [BcNodes3; conn3Tri(i,:)];
            nBC3 = size(unique(BcNodes3),1) ;
        elseif (Elem.Tag(i,1) == pnId2(4))
            BcNodes4 = [BcNodes4; conn3Tri(i,:)];
            nBC4 = size(unique(BcNodes4),1) ;
        elseif (Elem.Tag(i,1) == pnId2(5))
            BcNodes5 = [BcNodes5; conn3Tri(i,:)];
            nBC5 = size(unique(BcNodes5),1) ;
        elseif (Elem.Tag(i,1) == pnId2(6))
            BcNodes6 = [BcNodes6; conn3Tri(i,:)];
            nBC6 = size(unique(BcNodes6),1) ;
        elseif (Elem.Tag(i,1) == pnId2(7))
            BcNodes7 = [BcNodes7; conn3Tri(i,:)];
            nBC7 = size(unique(BcNodes7),1) ;
        elseif (Elem.Tag(i,1) == pnId2(8))
            BcNodes8 = [BcNodes8; conn3Tri(i,:)];
            nBC8 = size(unique(BcNodes8),1) ;
        elseif (Elem.Tag(i,1) == pnId2(9))
            BcNodes9 = [BcNodes9; conn3Tri(i,:)];
            nBC9 = size(unique(BcNodes9),1) ;
        elseif (Elem.Tag(i,1) == pnId2(10))
            BcNodes10 = [BcNodes10; conn3Tri(i,:)];
            nBC10 = size(unique(BcNodes10),1) ;
        end
    elseif (Elem.Type(i) == 3)
        if (Elem.Tag(i,1) == pnId2(1))
            BcNodes1 = [BcNodes1; conn4Quad(i,:)];
            nBC1 = size(unique(BcNodes1),1) ;
        elseif (Elem.Tag(i,1) == pnId2(2))
            BcNodes2 = [BcNodes2; conn4Quad(i,:)];
            nBC2 = size(unique(BcNodes2),1) ;
        elseif (Elem.Tag(i,1) == pnId2(3))
            BcNodes3 = [BcNodes3; conn4Quad(i,:)];
            nBC3 = size(unique(BcNodes3),1) ;
        elseif (Elem.Tag(i,1) == pnId2(4))
            BcNodes4 = [BcNodes4; conn4Quad(i,:)];
            nBC4 = size(unique(BcNodes4),1) ;
        elseif (Elem.Tag(i,1) == pnId2(5))
            BcNodes5 = [BcNodes5; conn4Quad(i,:)];
            nBC5 = size(unique(BcNodes5),1) ;
        elseif (Elem.Tag(i,1) == pnId2(6))
            BcNodes6 = [BcNodes6; conn4Quad(i,:)];
            nBC6 = size(unique(BcNodes6),1) ;
        elseif (Elem.Tag(i,1) == pnId2(7))
            BcNodes7 = [BcNodes7; conn4Quad(i,:)];
            nBC7 = size(unique(BcNodes7),1) ;
        elseif (Elem.Tag(i,1) == pnId2(8))
            BcNodes8 = [BcNodes8; conn4Quad(i,:)];
            nBC8 = size(unique(BcNodes8),1) ;
        elseif (Elem.Tag(i,1) == pnId2(9))
            BcNodes9 = [BcNodes9; conn4Quad(i,:)];
            nBC9 = size(unique(BcNodes9),1) ;
        elseif (Elem.Tag(i,1) == pnId2(10))
            BcNodes10 = [BcNodes10; conn4Quad(i,:)];
            nBC10 = size(unique(BcNodes10),1) ;
        end
    end
end
BCCyl = BcNodes1 ;
BCOutercyl = BcNodes2;
BCradial = BcNodes3;
BCFluid = BcNodes4;

% save('BCLeft1');
% save('BCLeft2');
% save('BCRight');
% save('BCTop');
% save('BCBottom');
% save('BCAll');
    
% for i=1:nPhysicalNames
%     if (i==1 && i~=nPhysicalNames)
%         filename = sprintf('%s/%s.%s.nbc',wrkDir,problemString,B1);
%         fileId = fopen(filename,'w');
%         if(Elem.Type(i)==1)
%             fprintf(fileId,'%d %d\n',BcNodes1);
%         elseif(Elem.Type(i)==2)
%             fprintf(fileId,'%d %d %d\n',BcNodes1);
%         elseif(Elem.Type(i)==3)
%             fprintf(fileId,'%d %d %d %d\n',BcNodes1);
%         end
%         fclose(fileId);
%     elseif (i==2 && i~=nPhysicalNames)
%         filename = sprintf('%s/%s.%s.nbc',wrkDir,problemString,B2);
%         fileId = fopen(filename,'w');
%         if(Elem.Type(i)==1)
%             fprintf(fileId,'%d %d\n',BcNodes2);
%         elseif(Elem.Type(i)==2)
%             fprintf(fileId,'%d %d %d\n',BcNodes2);
%         elseif(Elem.Type(i)==3)
%             fprintf(fileId,'%d %d %d %d\n',BcNodes2);
%         end
%         fclose(fileId);
%     elseif (i==3 && i~=nPhysicalNames)
%         filename = sprintf('%s/%s.%s.nbc',wrkDir,problemString,B3);
%         fileId = fopen(filename,'w');
%         if(Elem.Type(i)==1)
%             fprintf(fileId,'%d %d\n',BcNodes3);
%         elseif(Elem.Type(i)==2)
%             fprintf(fileId,'%d %d %d\n',BcNodes3);
%         elseif(Elem.Type(i)==3)
%             fprintf(fileId,'%d %d %d %d\n',BcNodes3);
%         end
%         fclose(fileId);
%     elseif (i==4 && i~=nPhysicalNames)
%         filename = sprintf('%s/%s.%s.nbc',wrkDir,problemString,B4);
%         fileId = fopen(filename,'w');
%         if(Elem.Type(i)==1)
%             fprintf(fileId,'%d %d\n',BcNodes4);
%         elseif(Elem.Type(i)==2)
%             fprintf(fileId,'%d %d %d\n',BcNodes4);
%         elseif(Elem.Type(i)==3)
%             fprintf(fileId,'%d %d %d %d\n',BcNodes4);
%         end
%         fclose(fileId);
%     elseif (i==5 && i~=nPhysicalNames)
%         filename = sprintf('%s/%s.%s.nbc',wrkDir,problemString,B5);
%         fileId = fopen(filename,'w');
%         if(Elem.Type(i)==1)
%             fprintf(fileId,'%d %d\n',BcNodes5);
%         elseif(Elem.Type(i)==2)
%             fprintf(fileId,'%d %d %d\n',BcNodes5);
%         elseif(Elem.Type(i)==3)
%             fprintf(fileId,'%d %d %d %d\n',BcNodes5);
%         end
%         fclose(fileId);
%     elseif (i==6 && i~=nPhysicalNames)
%         filename = sprintf('%s/%s.%s.nbc',wrkDir,problemString,B6);
%         fileId = fopen(filename,'w');
%         if(Elem.Type(i)==1)
%             fprintf(fileId,'%d %d\n',BcNodes6);
%         elseif(Elem.Type(i)==2)
%             fprintf(fileId,'%d %d %d\n',BcNodes6);
%         elseif(Elem.Type(i)==3)
%             fprintf(fileId,'%d %d %d %d\n',BcNodes6);
%         end
%         fclose(fileId);
%     elseif (i==7 && i~=nPhysicalNames)
%         filename = sprintf('%s/%s.%s.nbc',wrkDir,problemString,B7);
%         fileId = fopen(filename,'w');
%         if(Elem.Type(i)==1)
%             fprintf(fileId,'%d %d\n',BcNodes7);
%         elseif(Elem.Type(i)==2)
%             fprintf(fileId,'%d %d %d\n',BcNodes7);
%         elseif(Elem.Type(i)==3)
%             fprintf(fileId,'%d %d %d %d\n',BcNodes7);
%         end
%         fclose(fileId);
%     elseif (i==8 && i~=nPhysicalNames)
%         filename = sprintf('%s/%s.%s.nbc',wrkDir,problemString,B8);
%         fileId = fopen(filename,'w');
%         if(Elem.Type(i)==1)
%             fprintf(fileId,'%d %d\n',BcNodes8);
%         elseif(Elem.Type(i)==2)
%             fprintf(fileId,'%d %d %d\n',BcNodes8);
%         elseif(Elem.Type(i)==3)
%             fprintf(fileId,'%d %d %d %d\n',BcNodes8);
%         end
%         fclose(fileId);
%     elseif (i==9 && i~=nPhysicalNames)
%         filename = sprintf('%s/%s.%s.nbc',wrkDir,problemString,B9);
%         fileId = fopen(filename,'w');
%         if(Elem.Type(i)==1)
%             fprintf(fileId,'%d %d\n',BcNodes9);
%         elseif(Elem.Type(i)==2)
%             fprintf(fileId,'%d %d %d\n',BcNodes9);
%         elseif(Elem.Type(i)==3)
%             fprintf(fileId,'%d %d %d %d\n',BcNodes9);
%         end
%         fclose(fileId);
%     elseif (i==10 && i~=nPhysicalNames)
%         filename = sprintf('%s/%s.%s.nbc',wrkDir,problemString,B10);
%         fileId = fopen(filename,'w');
%         if(Elem.Type(i)==1)
%             fprintf(fileId,'%d %d\n',BcNodes10);
%         elseif(Elem.Type(i)==2)
%             fprintf(fileId,'%d %d %d\n',BcNodes10);
%         elseif(Elem.Type(i)==3)
%             fprintf(fileId,'%d %d %d %d\n',BcNodes10);
%         end
%         fclose(fileId);
%     elseif ( i==nPhysicalNames )
%         filename = sprintf('%s/%s.all.nbc',wrkDir,problemString);
%         fileId = fopen(filename,'w');
%         for j=1:nNodes
%             fprintf(fileId,'%d\n',j);
%         end
%         fclose(fileId);
%     end
% end

% filename = sprintf('%s/%s.meshData',wrkDir,problemString);
% fileId = fopen(filename,'w');
% fprintf(fileId,'nCrds %d\n',nNodes);
% fprintf(fileId,'nCnns %d\n',a-1);
% fprintf(fileId,'nPhysicalNames %d\n',nPhysicalNames);
% for i=1:nPhysicalNames
%     if (i==1&& i~=nPhysicalNames) fprintf(fileId,'%s %d\n',B1,nBC1);
%     elseif (i==2&& i~=nPhysicalNames) fprintf(fileId,'%s %d\n',B2,nBC2);
%     elseif (i==3&& i~=nPhysicalNames) fprintf(fileId,'%s %d\n',B3,nBC3);
%     elseif (i==4&& i~=nPhysicalNames) fprintf(fileId,'%s %d\n',B4,nBC4);
%     elseif (i==5&& i~=nPhysicalNames) fprintf(fileId,'%s %d\n',B5,nBC5);
%     elseif (i==6&& i~=nPhysicalNames) fprintf(fileId,'%s %d\n',B6,nBC6);
%     elseif (i==7&& i~=nPhysicalNames) fprintf(fileId,'%s %d\n',B7,nBC7);
%     elseif (i==8&& i~=nPhysicalNames) fprintf(fileId,'%s %d\n',B8,nBC8);
%     elseif (i==9&& i~=nPhysicalNames) fprintf(fileId,'%s %d\n',B9,nBC9);
%     elseif (i==10&& i~=nPhysicalNames) fprintf(fileId,'%s %d\n',B10,nBC10);
%     end
% end
% fclose(fileId);

save Mesh.mat BCCyl  crd BCFluid BCOutercyl BCradial
