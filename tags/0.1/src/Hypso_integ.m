function  HI_rout()
% Hypsometric Integral calculator
% By Sagy Cohen (sagy.cohen@uon.edu.au)
% Calculate the Hypsometric Integral for each pixel at catchment.  
% Each pixel is considered as an outlet for its contributing area.
% Requires two input files (ESRI ASCII format):
%   (1) DEM - defined in line 19
%   (2) Flow direction grid - defined in line 20
% Produce 5 output files (ESRI ASCII format) names can be modified at the 
% end of this file:
%   (1) HI.txt - pixel scale hypsometric integral;
%   (2) max_elev.txt - the maximum elevation of the catchment flowing
%   thorough each pixel;
%   (3) Elev_Acc.txt - the sum of the elevation (m) of all the pixels
%   flowing thorough each pixel;
%   (4) flowacc.txt - Contributing area in pixels;
%   (5) junctions.txt - how many of a pixel's 8 neighbors flow into it;
 
% To run with different dataset change the file name here:  
newData1 = importdata('DEM_demo.txt', ' ', 6); %DEM ASCII FILE
flow_dir=int16(arcgridread('FlowDir_demo.txt'));   %Flow-Direction ASCII file
tmpdir=flow_dir;

% change flow-direction coding from ArcGIS to TauDEM
if max(max(flow_dir,[],1))==128
     tmpdir(flow_dir==2)=8;
     tmpdir(flow_dir==4)=7;
     tmpdir(flow_dir==8)=6;
     tmpdir(flow_dir==16)=5;
     tmpdir(flow_dir==32)=4;
     tmpdir(flow_dir==64)=3;
     tmpdir(flow_dir==128)=2;
end
flow_dir=tmpdir;
% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end
elev=newData1.data;
%FlowDirectionCalc(elev);
header=char(newData1.textdata);
N= [0 -1;1 -1;1 0;1 1;0 1;-1 1;-1 0;-1 -1;0 0];
D= [0 1;-1 1;-1 0;-1 -1;0 -1;1 -1;1 0;1 1];
[r c]=size(elev);
area=ones(r,c);
ElevAcc=elev;
junctions=zeros(r,c);
flow_dir(flow_dir(:,1)==4)=0;flow_dir(flow_dir(:,1)==5)=0;...
    flow_dir(flow_dir(:,1)==6)=0;
tmp=flow_dir(:,end);
tmp(flow_dir(:,end)==8)=0;tmp(flow_dir(:,end)==2)=0;...
    tmp(flow_dir(:,end)==1)=0;
flow_dir(:,end)=tmp;

tmp=flow_dir(1,:);
tmp(flow_dir(1,:)==2)=0;tmp(flow_dir(1,:)==3)=0;tmp(flow_dir(1,:)==4)=0;
flow_dir(1,:)=tmp;

tmp=flow_dir(end,:);
tmp(flow_dir(end,:)==6)=0;tmp(flow_dir(end,:)==7)=0;tmp(flow_dir(end,:)==8)=0;
flow_dir(end,:)=tmp;

tmpJunc=junctions;

for i=1:8
    Flows8(:,:,i)=circshift(flow_dir,[N(i,1) N(i,2)]);
end

%claculate a raster which contain the number of pixels flowing in to each
%pixel (junctions)
junctions=(Flows8(:,:,1)==5)+(Flows8(:,:,2)==6)+(Flows8(:,:,3)==7)...
    +(Flows8(:,:,4)==8)+(Flows8(:,:,5)==1)+ (Flows8(:,:,6)==2)+...
    (Flows8(:,:,7)==3)+ (Flows8(:,:,8)==4);

tmpJunc=junctions;
MaxElev=elev;
%Create an array with the index of source pixels (junction 0 and
%elevation larger then zero
JuncN0=find(junctions==0);%&elev>0);

%Create a sorted array of source pixels index and thier elevation 
SourceIZ=flipdim(sortrows([JuncN0,elev(JuncN0)],2),1);
SourceI=SourceIZ(1:size(SourceIZ))';
i=1;
while sum(SourceI)>0
sum(SourceI)
    %for i=1:size(SourceI)

if SourceI(i)>0 
       ext=1;
       dir=flow_dir(SourceI(i));
       if dir==0
           SourceI(i)=0;
           SourceI=sort(SourceI,'descend');
           i=i+1;
           continue
       end
       %extract the x,y index from linear index
       [tr tc]=ind2sub([r c],SourceI(i)); 
       RCindex=[tr tc];
       %Calculate the location index of the downflow pixel
       FlowIndx=sub2ind([r c],RCindex(1)+D(dir,1), RCindex(2)+ D(dir,2));
         
       if FlowIndx>0 && FlowIndx<(r*c)&&SourceI(i)>0&&SourceI(i)<(r*c)
           %while tmpJunc(FlowIndx)<3&&ext==1&&tmpJunc(SourceI(i))>=0%%%%%
           while ext==1&&tmpJunc(SourceI(i))==0
               dir=flow_dir(SourceI(i));
               if dir==0
                  SourceI(i)=0;
                  SourceI=sort(SourceI,'descend');
                  break
               end
               [tr tc]=ind2sub([r c],SourceI(i)); 
               RCindex=[tr tc];
               FlowIndx=sub2ind([r c],RCindex(1)+D(dir,1), RCindex(2)+ D(dir,2));

               if MaxElev(SourceI(i))> MaxElev(FlowIndx)
                   MaxElev(FlowIndx)=MaxElev(SourceI(i));
               end
               
               area(FlowIndx)= area(FlowIndx)+area(SourceI(i));
               ElevAcc(FlowIndx)=ElevAcc(FlowIndx)+ElevAcc(SourceI(i));
               tmpJunc(FlowIndx)=tmpJunc(FlowIndx)-1;%%%
               SourceI(i)=FlowIndx;

               if sum(SourceI==FlowIndx)>1%%%
                   SourceI(i)=0;
                   ext=2;
               end
               if tmpJunc(FlowIndx)>0 %||tmpJunc(SourceI(i))>0
                   ext=2;
               end
           end
       end
   %    waitbar(sum(SourceI)\sum(SourceIZ(:,1)))
end
    SourceI=sort(SourceI,'descend');
    if i==size(SourceI,1)||SourceI(i)==0
       i=1;
   else
    i=i+1;
   end
end
%Calculate the Hypsometric Integral
MeanElev=ElevAcc./area;
HI=(MeanElev-elev)./(MaxElev-elev);

HI(isnan(HI))=-9999; % replace 'NaN' with -9999

%Write output files
dlmwrite('HI.txt', header,'%F');
dlmwrite('HI.txt', HI,'-append','roffset', 0, 'delimiter', ' ');
%HIOut=arcgridread('HI.txt');
%assignin('base', 'HIOut', HIOut);

dlmwrite('max_elev.txt', header,'%F');
dlmwrite('max_elev.txt', MaxElev,'-append','roffset', 0, 'delimiter', ' ');
%MaxElevOut=arcgridread('max_elev.txt');
%assignin('base', 'MaxElevOut', MaxElevOut);

dlmwrite('Elev_Acc.txt', header,'%F');
dlmwrite('Elev_Acc.txt', ElevAcc,'-append','roffset', 0, 'delimiter', ' ');
%Elev_Acc=arcgridread('Elev_Acc.txt');
%assignin('base', 'Elev_Acc', Elev_Acc);

dlmwrite('flowacc.txt', header,'%F');
dlmwrite('flowacc.txt', area,'-append','roffset', 0, 'delimiter', ' ');
%flow_acc=arcgridread('flowacc.txt');
%assignin('base', 'flow_acc', flow_acc);

dlmwrite('junctions.txt', header,'%F');
dlmwrite('junctions.txt', junctions,'-append','roffset', 0, 'delimiter', ' ');
%junctions=arcgridread('junctions.txt');
%assignin('base', 'junctions', junctions);
end

