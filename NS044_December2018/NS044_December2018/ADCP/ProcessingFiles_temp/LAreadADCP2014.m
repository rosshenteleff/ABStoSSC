
% Read and plot wave and current data from vector

%% add Matlab function yearday to directory

%addpath('D:/data');
%addpath('D:\Data\Matlab_Functions');
%disp('Path Generated');

%% Variables to change
% filename =('D:\2012\ADCP\May\');
% day =('AMAY6A01');

%% CHANGE THIS DIRECTORY!
filename =('G:\MATLAB\ADCP\');

day =('AY25A01');
numberofbins = 17;
disp('Variables Uploaded');
%colormapoutward = importdata('G:\colormapoutward.txt');

% KP.Velo = load ('KP_C_Velo_10m_18L_J8_12_2013','-mat');
% KP.WaterLevel = load ('KP_C_WL_10m_18L_J8_12_2013','-mat');
julianday = 125.8125;
% julianday2= 164.43;
% julianday3= 166;
% h_interval = 5;
% layers= 18;
% fakelayers=300;
%% Files and importing



                  %REMEMBER TO CHANGE FILENAMES
sensor =[filename,day,'.sen'];
header2 =[filename,day,'.hr2'];
depth =[filename,day,'.d1'];
velocity1 =[filename,(day),'.v1'];
velocity2 =[filename,(day),'.v2'];
velocity3 =[filename,(day),'.v3'];
amplitude1 =[filename,(day),'.a1'];
amplitude2 =[filename,(day),'.a2'];
amplitude3 =[filename,(day),'.a3'];
correlation1 =[filename,(day),'.c1'];
correlation2 =[filename,(day),'.c2'];
correlation3 =[filename,(day),'.c3'];
disp('Filenames Located');


%% Import data

data.sensor = importdata(sensor);
data.header2 = importdata(header2);
data.depthfile = importdata(depth);
data.v1 = importdata(velocity1);
data.v2 = importdata(velocity2);
data.v3 = importdata(velocity3);
data.a1 = importdata(amplitude1);
data.a2 = importdata(amplitude2);
data.a3 = importdata(amplitude3);
data.c1 = importdata(correlation1);
data.c2 = importdata(correlation2);
data.c3 = importdata(correlation3);
disp('Data Imported');
%% sensor data (date, orientation, pressure, temperature and water depth)
data.mo =  data.sensor(:,1);
data.dy =  data.sensor(:,2);
data.yr =  data.sensor(:,3);
data.hr =  data.sensor(:,4);
data.mn =  data.sensor(:,5);
data.sc =  data.sensor(:,6);
data.yd =  yearday(data.yr,data.mo,data.dy + data.hr/24 + data.mn/(24*60)+data.sc/(24*3600));
data.heading =data.sensor(:,13);
data.pitch =data.sensor(:,14);
data.roll =data.sensor(:,15);
data.pre = data.sensor(:,16);
data.temp = data.sensor(:,17);
data.dep = data.pre*10000/(9810);
disp('Sensor Array Complete');
%% header data (average velocity, amplitude and correlation)
data.ave = data.header2 (:,7);
data.avn = data.header2 (:,8);
data.avv = data.header2 (:,9);

data.aae = data.header2 (:,10);
data.aan = data.header2 (:,11);
data.aav = data.header2 (:,12);

data.ace = data.header2 (:,13);
data.acn = data.header2 (:,14);
data.acv = data.header2 (:,15);

data.fakemag = (data.ave.^2+ data.avn.^2).^(.5);
disp('Header Array Complete');
%% depth of bins (vertical and beam bins)
%NOTE: This data is taken from the hdr file and made into a seperate file
%txt file named "filename".dl  (ie Aug_11.d1).
data.dbe  = data.depthfile (:,2);
data.dbn  = data.depthfile (:,2);
data.dbv  = data.depthfile (:,3);
disp('Depth of Bins Array Complete');
%% Velocities, amplitudes and correlation files have no need to be made into
%structures.
data.mag = sqrt(data.v1(:,3:99).^2 + data.v2(:,3:99).^2);
%data.mag3d = sqrt(data.v1.^2 + data.v2.^2+ data.v3.^2);
disp('Velocity Array Complete');

%% filtered data for a data set above 80%
data.fv1=data.v1;
data.fv2=data.v2;
data.fv3=data.v3;
i=1;
j=1;
count1=0;
L = length (data.sensor);
onestep = ones(L,1);
data.count = onestep-1;
disp('Start first while loop');
while j <= L
    row =j;
    %disp(j)
    while i <= numberofbins
        column =i;
        
        % %%%%%%%%%%%%%% filtered data for a data set above 80%
        if data.c1((row),((column)+2)) < 80
            data.fv1((row),((column)+2))= NaN;         
        end    
        if data.c2((row),((column)+2)) < 80
            data.fv2((row),((column)+2))= NaN;
        end    
        if data.c3((row),((column)+2)) < 80
            data.fv3((row),((column)+2))= NaN;
        end
         
        % %%%%%%%%%%%%%% Direction 2D 
        % North taken as 0 degrees; positive is clockwise
        %north east quadrant
        if data.v1((row),((column)+2))>0 && data.v2((row),((column)+2))>0            
            data.angle2d((row),(column)) = atand(data.v2((row),((column)+2))/data.v1((row),((column)+2)));   
        end
        %Southeast quadrant
        if data.v1((row),((column)+2))<0 && data.v2((row),((column)+2))>0
            data.angle2d((row),(column)) = (atand(data.v2((row),((column)+2))/data.v1((row),((column)+2)))+90);  
        end
        %southwest quadrant
        if data.v1((row),((column)+2))<0 && data.v2((row),((column)+2))<0
            data.angle2d((row),(column)) = (atand(data.v2((row),((column)+2))/data.v1((row),((column)+2)))+180); 
        end
        %northwest quadrant
        if data.v1((row),((column)+2))>0 && data.v2((row),((column)+2))<0
            data.angle2d((row),(column)) = (atand(data.v2((row),((column)+2))/data.v1((row),((column)+2)))+270);
        end
        
        % %%%%%%%%%%%%% Direction 3D
        %up is positive
        if data.v3((row),((column)+2))>0
           data.angle3d((row),(column)) = atand(data.v3((row),((column)+2))/data.mag((row),((column))));  
        end
        if data.v3((row),((column)+2))<0
           data.angle3d((row),(column)) = -(atand(data.v3((row),((column)+2))/data.mag((row),((column)))));  
        end
        
        % # of bins for calculation for Depth averaged velocity
        if data.dep((row),:) > data.dbe((column),:)
            count1=count1+1;
            data.count(row,1) = count1;
        elseif data.dep((row),:) < data.dbe((column),:)
            data.count(row,1) = data.count(row,1)+0;
        end     
        
        
        % %%%%%%%%%%%%% Cutting data above water level
       if data.dep((row),:) < data.dbe((column),:); 
            data.v1((row),((column)+2))=NaN;
            data.a1((row),((column)+2))=NaN;
            data.c1((row),((column)+2))=NaN;
            data.fv1((row),((column)+2))=NaN;
            data.mag((row),((column)))=NaN;
            
        end
        if data.dep((row),:) < data.dbn((column),:);
            data.v2((row),((column)+2))=NaN;
            data.a2((row),((column)+2))=NaN;
            data.c2((row),((column)+2))=NaN;
            data.fv2((row),((column)+2))=NaN;
        end
        if data.dep((row),:) < data.dbv((column),:); 
            data.v3((row),((column)+2))=NaN;
            data.a3((row),((column)+2))=NaN;
            data.c3((row),((column)+2))=NaN;
            data.fv3((row),((column)+2))=NaN;
            data.angle2d((row),((column)))= NaN;
            data.angle3d((row),((column)))= NaN;
        end
        i=i+1;
    end
    % Calculation of Depth Average Velocity
    data.dav(row,1) = (sum(data.mag(row,1:(data.count(row,1))))/data.count(row,1));
    count1=0;
    i=1;
    j=j+1;
    if j == round (L/4) 
        disp('25%...')
    end
    if j == round (L/2) 
        disp('50%...')
    end
     if j == round (3*L/4) 
        disp('75%...')
    end
end 
disp('While Loop Finished');

%% Magnitude of Filtered Velocity
data.fmag = (data.fv1(:,3:99).^2+ data.fv2(:,3:99).^2).^(0.5);
%data.fmag=[];
i=1;
count1=0;
j=1;
summagnitude=0;
disp('Start Second while loop')
while j <= L
    row =j;
    %disp(j)
    while i <= numberofbins
        column =i;        
        if data.dep((row),:) > data.dbe((column),:) && data.fmag(row,column) > 0 
            count1=count1+1;
            data.count(row,1) = count1;
            summagnitude = summagnitude + data.fmag(row,column);
        elseif data.dep((row),:) < data.dbe((column),:)
            data.count(row,1) = data.count(row,1)+0;
        end    
    i=i+1;    
    end
    % Calculation of Depth Average Velocity
    data.fdav(row,1) = summagnitude/data.count(row,1);
    count1=0;
    summagnitude =0;
    i=1;
    j=j+1;
    if j == round (L/4) 
        disp('25%...')
    end
    if j == round (L/2) 
        disp('50%...')
    end
     if j == round (3*L/4) 
        disp('75%...')
    end
end
disp('While Loop Finished');
%% Colormap for filtered map

fig = colormap;
fig(1,:) = 1;
%fig(32,:) = 1;
colormap(fig);

colormap1= colormap;
colormap2= colormap;
colormap1(1:31,3)= colormap2(1:31,1);
colormap1(1:31,1)= colormap2(1:31,3);
disp('Colormap Array Complete');
%% Model calculations starts

% Add a time length to model run
% KP.Velo.data.Time = [];
% KP.Velo.data.Time(1,1)= 0;
% KP.WaterLevel.data.Time(1,1)= 0;
% aa=1;
% while aa <= length(KP.Velo.data.XComp)
%     i=aa-1;    
%     KP.Velo.data.Time(aa,:)= h_interval*i/60/24;
%     KP.WaterLevel.data.Time(aa,:)= h_interval*i/60/24;
%     aa=aa+1;
% end
% KP.Velo.data.Time = KP.Velo.data.Time+julianday;
% KP.WaterLevel.data.Time = KP.WaterLevel.data.Time+julianday;
% % KP.ISNAN(:,:) = isnan(KP.Velo.data.Z(:,1,:));
% % Calculate magnitude/ Add depth to sigma layers
% KP.Velo.data.S = (KP.Velo.data.XComp.^2 + KP.Velo.data.YComp.^2).^(1/2);
% KP.Velo.data.Z = KP.Velo.data.Z + 0.5379;
% 
% %% Rewrite files to include NaN Value
% for i=1:length(KP.WaterLevel.data.Val)
%     if KP.Velo.data.Z(i,1,1)<.5379
%         KP.Velo.data.Z(i,1,:)=0;
%     end     
%     for j=1:layers
%         
%         KP.Z(i,j) = KP.Velo.data.Z(i,1,j);
%         KP.X(i,j) = KP.Velo.data.XComp(i,1,j);
%         KP.Y(i,j) = KP.Velo.data.YComp(i,1,j);
%         KP.S(i,j) = KP.Velo.data.S(i,1,j); 
% %         if KP.ISNAN(i,j) == 1
% %             KP.S(i,j)=0;
% %         end        
%         
% %         KP.XComp(i,j+2) = KP.Velo.data.XComp(i,1,j);
% %         KP.YComp(i,j+2) = KP.Velo.data.YComp(i,1,j);
% %         KP.S(i,j+2) = KP.Velo.data.S(i,1,j); 
% %         if j==layers
% % %             KP.Z(i,j+1)=(KP.Z(i,j)-KP.Z(i,j-1))+KP.Z(i,j);
% %             KP.Z(i,j+1)=0;
% %             
% %             KP.Z(i,j+2)=0;
% %             KP.S(i,j+2)=0;
% %             
% % %             KP.Z(i,1)=0;
% %             KP.S(i,1)=0;
% %             KP.S(i,2)=0;
% %         end
%     end
%     
%     if i == round (length(KP.WaterLevel.data.Val)/4) 
%         disp('25%...')
%     end
%     
%     if i == round (length(KP.WaterLevel.data.Val)/2) 
%         disp('50%...')
%     end
%     
%      if i == round (3*length(KP.WaterLevel.data.Val)/4) 
%         disp('75%...')
%         disp('First Loop Complete');
%     end
% end
% KP.FakeZ=KP.Z;
% KP.T = KP.WaterLevel.data.Time;
% % KP.Z(1+length(KP.WaterLevel.data.Val),:)= KP.Z(length(KP.WaterLevel.data.Val),:);
% % KP.S(1+length(KP.WaterLevel.data.Val),:)= 0;
% % KP.T(1+length(KP.WaterLevel.data.Val),:) = KP.T(length(KP.WaterLevel.data.Val),:)+h_interval/60/24;
% 
% %% Create fake layers
% 
% KPL.Z= KP.Z;
% KPL.T= KP.T;
% size = length(KP.S);
% KPL.S= ones(size,fakelayers);
% KPL.Z= ones(size,fakelayers);
% KPL.X= ones(size,fakelayers);
% KPL.Y= ones(size,fakelayers);
% 
% 
% for i=1:size
%    for j=1:fakelayers
%        KPL.Z(i,j)=9-(j-1)*.03;
%         for c=1:(layers-1)
%            
%             if KPL.Z(i,j)< KP.Z(i,c) && KPL.Z(i,j)>KP.Z(i,c+1);
%                 KPL.S(i,j) =KP.S(i,c);
%                 KPL.X(i,j) =KP.X(i,c);
%                 KPL.Y(i,j) =KP.Y(i,c);
%             end
%             if c== (layers-1)
%                 if KPL.Z(i,j)<= KP.Z(i,c) && KPL.Z(i,j)>=0;
%                     KPL.S(i,j) =KP.S(i,c);
%                     KPL.X(i,j) =KP.X(i,c);
%                     KPL.Y(i,j) =KP.Y(i,c);
%                 end
%             end
%            
%             if KPL.Z(i,j)> KP.WaterLevel.data.Val(i,1)
%                 KPL.S(i,j) =0;
%                 KPL.X(i,j) =-5;
%                 KPL.Y(i,j) =-5;
%             end
%             
%         end
%         
%    end
%     if KPL.S(i,1) ==1 && KPL.S(i,fakelayers)==1
%         KPL.S(i,:) = 0;
%         KPL.X(i,j) = -5;
%         KPL.Y(i,j) = -5;
%     end
% end

%% Write Imagesc First way
% clims = [ 0 .4 ];
% figure(2)

%  for i=5000:length(KP.T)
%     hold on
%         imagesc(KP.T(i,1),KP.Z(i,:),KP.S(i,:)',clims)
%     hold on
%     if i == round (length(KP.T)/4) 
%         disp('25%...')
%     end
%     if i == round (length(KP.T)/2) 
%         disp('50%...')
%     end
%      if i == round (3*length(KP.T)/4) 
%         disp('75%...')
%         disp('PLotting Loop Complete');
%     end
% end



%% Write Imagesc 


% clims = [ 0 .4 ];
% figure(1)
% colormap(fig);
% set(gca,'YDir','Normal')
% imagesc(KPL.T,KPL.Z(1,:),KPL.S',clims)
% set(gca,'YDir','Normal')
% hold on
% %% plot waterlevel line/ rewrite time scale
% KP.FakeTime = KP.WaterLevel.data.Time;
% KP.FakeVal  = KP.WaterLevel.data.Val;
% plot (KP.FakeTime,KP.FakeVal,'-k')
% hold on
% shading interp
% set(gca,'YDir','Normal')
% %set(gca,'YLim',[0 5])
% %set(gca,'XLim','tight')
% axis ([159 164 0 8])
% ylabel('Water Level (m)')
% set(gca,'CLim',[-0 .5])
% title('Filtered Velocity Easting component ');
% colorbar
% % colormap(fig);
% text(186.18,9.5,'u (m/s)');





















% %% Figure 1 Modelled
% 
% figure(1)
% subplot(311)
%     hold on; box on;
%     imagesc(KPL.T,KPL.Z(1,:),KPL.X',clims)
%     plot (KP.FakeTime,KP.FakeVal,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
% %     axis ([161.975 162.275 0 9])
%     ylabel('z (m)')
%     set(gca,'CLim',[-.5 .5])
%     title('Modelled Velocity E component ');
%     colorbar
%     colormap(fig);
%     
%  subplot(312)
%     hold on; box on;
%     imagesc(KPL.T,KPL.Z(1,:),KPL.Y',clims)
%     plot (KP.FakeTime,KP.FakeVal,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
% %     axis ([161.975 162.275 0 9])
%     ylabel('z (m)')
%     set(gca,'CLim',[-.5 .5])
%     title('Modelled Velocity N component ');
%     colorbar
%    colormap(fig);
%  
% subplot(313)
%     hold on; box on;
%     imagesc(KPL.T,KPL.Z(1,:),KPL.S',clims)
%     plot (KP.FakeTime,KP.FakeVal,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
% %    axis ([161.975 162.275 0 9])
%     ylabel('z (m)')
%     set(gca,'CLim',[0 .5])
%     title('Modelled Speed component ');
%     colorbar    
%     colormap(fig);
%  %  print -dpdf PCS_validation
% 
%  %% Figure 2 Combination
% 
% figure(2)
% subplot(321)
%     hold on; box on;
%     imagesc(KPL.T,KPL.Z(1,:),KPL.X',clims)
%     plot (KP.FakeTime,KP.FakeVal,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     axis ([161.975 162.275 0 9])
%     ylabel('z (m)')
%     set(gca,'CLim',[-.5 .5])
%     title('Modelled Velocity E component ');
%     colorbar
%     colormap(fig);
%     
%  subplot(323)
%     hold on; box on;
%     imagesc(KPL.T,KPL.Z(1,:),KPL.Y',clims)
%     plot (KP.FakeTime,KP.FakeVal,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     axis ([161.975 162.275 0 9])
%     ylabel('z (m)')
%     set(gca,'CLim',[-.5 .5])
%     title('Modelled Velocity N component ');
%     colorbar
%    colormap(fig);
%  
% subplot(325)
%     hold on; box on;
%     imagesc(KPL.T,KPL.Z(1,:),KPL.S',clims)
%     plot (KP.FakeTime,KP.FakeVal,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     axis ([161.975 162.275 0 9])
%     ylabel('z (m)')
%     set(gca,'CLim',[0 .5])
%     title('Modelled Speed component ');
%     colorbar    
%     colormap(fig);
% 
% subplot(322)
%     hold on; box on;
%     imagesc(data.yd, data.dbe,data.v1(:,3:99)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     axis ([161.975 162.275 0 9])
%     ylabel('z (m)')
%     set(gca,'CLim',[-.5 .5])
%     title('Non Filtered data Velocity E component ');
%     colorbar
%     colormap(fig);
%     
%  subplot(324)
%     hold on; box on;
%     imagesc(data.yd, data.dbn,data.v2(:,3:99)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     axis ([161.975 162.275 0 9])
%     ylabel('z (m)')
%     set(gca,'CLim',[-.5 .5])
%     title('Non Filtered data Velocity N component ');
%     colorbar
%    colormap(fig);
%  
% subplot(326)
%     hold on; box on;
%     imagesc(data.yd, data.dbv,data.mag(:,1:97)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     axis ([161.975 162.275 0 9])
%     ylabel('z (m)')
%     set(gca,'CLim',[0 .5])
%     title('Non Filtered data Speed component ');
%     colorbar    
%     colormap(fig);
%  %  print -dpdf PCS_validation
% [ax,h3]=suplabel('Modelled vs Data Velocity Arc Analysis Plot' ,'t');

%% Figure 3 Non Filtered data Velocity

% figure(3)
% subplot(311)
%     hold on; box on;
%     imagesc(data.yd, data.dbe,data.v1(:,3:99)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([126.9 127.2 0 9])
%     ylabel('z (m)')
%     set(gca,'CLim',[-.3 .3])
%     title('Non Filtered data Velocity E component ');
%     colorbar
%     colormap(fig);
%     
%  subplot(312)
%     hold on; box on;
%     imagesc(data.yd, data.dbn,data.v2(:,3:99)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([126.9 127.2 0 9])
%     ylabel('z (m)')
%     set(gca,'CLim',[-.3 .3])
%     title('Non Filtered data Velocity N component ');
%     colorbar
%    colormap(fig);
%  
% subplot(313)
%     hold on; box on;
%     imagesc(data.yd, data.dbv,data.v3(:,3:99)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([126.9 127.2 0 9])
%     ylabel('z (m)')
%     set(gca,'CLim',[-.1 .1])
%     title('Non Filtered data Velocity V component ');
%     colorbar    
%     colormap(fig);
%  %  print -dpdf PCS_validation
% 
% %% Figure 4 Filtered data Velocity
% 
% figure(4)
% subplot(311)
%     
%     hold on; box on;
%     imagesc(data.yd, data.dbe,data.fv1(:,3:99)')
%     plot(data.yd,data.dep,'-k') 
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([221 127.2 0 9])
%     ylabel('Water Level (m)')
%     set(gca,'CLim',[-.3 .3])
%     title('Filtered Velocity Easting component ');
%     colorbar
%     colormap(fig);
%     text(186.18,3.5,'u (m/s)');
%     
%  subplot(312)
%     hold on; box on;
%     imagesc(data.yd, data.dbn,data.fv2(:,3:99)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([126.9 127.2 0 9])
%     ylabel('Water Level (m)')
%     set(gca,'CLim',[-.3 .3])
%     title('Northing component ');
%     colorbar
%     colormap(fig);
%     text(186.18,3.5,'v (m/s)');
%  
% subplot(313)
%     hold on; box on;
%     imagesc(data.yd, data.dbv,data.fv3(:,3:99)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([126.9 127.2 0 9])
%     ylabel('Water Level (m)')
%     set(gca,'CLim',[-0.1 .1])
%     title('Vertical component ');
%     colorbar
%     colormap(fig);
%     text(186.18,3.5,'w (m/s)');
 %  print -dpdf PCS_validation
%% figure 5 Amplitude

%figure(5)
% subplot(311)
%     hold on; box on;
%     imagesc(data.yd, data.dbe,data.a1(:,3:99)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([223.55 223.85 0 5])
%     ylabel('z (m)')
%     set(gca,'CLim',[0 200])
%     title('Amplitude E component ');
%     colorbar
%     colormap(fig);
%     
%  subplot(312)
%     hold on; box on;
%     imagesc(data.yd, data.dbn,data.a2(:,3:99)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([223.55 223.85 0 5])
%     ylabel('z (m)')
%     set(gca,'CLim',[0 200])
%     title('Amplitude N component ');
%     colorbar
%     colormap(fig);
%  
% subplot(313)
%     hold on; box on;
%     imagesc(data.yd, data.dbv,data.a3(:,3:99)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([223.55 223.85 0 5])
%     ylabel('z (m)')
%     set(gca,'CLim',[0 200])
%     title('Amplitude V component ');
%     colorbar
%     colormap(fig);
%     
%     
%     %% figure 6 Correlation
% 
% figure(6)
% subplot(311)
%     hold on; box on;
%     imagesc(data.yd, data.dbe,data.c1(:,3:99)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([223.55 223.85 0 5])
%     ylabel('z (m)')
%     set(gca,'CLim',[0 100])
%     title('Correlation E component ');
%     colorbar
%     colormap(fig);
%     
%  subplot(312)
%     hold on; box on;
%     imagesc(data.yd, data.dbn,data.c2(:,3:99)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([223.55 223.85 0 5])
%     ylabel('z (m)')
%     set(gca,'CLim',[0 100])
%     title('Correlation N component ');
%     colorbar
%     colormap(fig);
%  
% subplot(313)
%     hold on; box on;
%     imagesc(data.yd, data.dbv,data.c3(:,3:99)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([223.55 223.85 0 5])
%     ylabel('z (m)')
%     set(gca,'CLim',[0 100])
%     title('Correlation V component ');
%     colorbar 
%     colormap(fig);
%     
    %% figure 7 Magnitude of velocity

% figure(7)
% 
%     hold on; box on;
%     imagesc(data.yd, data.dbe,data.mag(:,1:97)')
%     %imagesc(data.depthfile(:,3)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     axis ([223.55 223.85 0 5])
%     ylabel('z (m)')
%     set(gca,'CLim',[0 .2])
%     title('Magnitude of velocity');
%     colorbar
%     colormap(jet);
%     
%   
%     
%     %% figure 8 ploting Velocities
%  figure (8)   
%  hold on;  box on;
%  plot(data.yd,data.mag,'-b') 
%  title('Magnitude of velocities ');
%   xlabel('Yearday')
%   ylabel('Velocity (m/s')
%   axis ([223.6 223.8 0 .5])
%  
%     %% figure 9 Direction
% 
% %figure(9)
% %subplot(111)
% %    hold on; box on;
%     
%     %imagesc(data.yd, data.dbv,data.c3(:,3:99)')
%     %plot(data.yd,data.dep,'-k')
%     
% %    imagesc(data.yd, data.angle2d(:,1:97)',data.mag(:,1:97)')
% %    plot(data.yd,data.dep,'-k')
% %    shading interp
% %    set(gca,'YDir','normal')
% %    %set(gca,'YLim',[0 5])
% %    %set(gca,'XLim','tight')
% %    axis ([223.6 223.8 0 360])
% %    ylabel('z (m)')
% %    set(gca,'CLim',[0 .2])
% %    title('E velocity ');
% %    colorbar
% %    colormap(fig);
%   
% %% figure 10
% 
% figure (10)
% subplot (311)
% hold on;  box on;
%  plot(data.yd,data.dav,'-b') 
%  title('Non filtered depth averaged magnitude of velocity ');
%   xlabel('Yearday')
%   ylabel('Velocity (m/s')
%   axis ([223.6 223.8 0 .5])
% subplot (312)
% hold on;  box on;
%  plot(data.yd,data.fdav,'-r') 
%  title('Filtered depth averaged magnitude of velocity ');
%   xlabel('Yearday')
%   ylabel('Velocity (m/s')
%   axis ([223.6 223.8 0 .5])  
% subplot (313)
% hold on;  box on;
%  plot(data.yd,data.dav,'-b')
%  plot(data.yd,data.fdav,'-r')
%  legend('data.dav =black','data.fdav =red')
%  title('Depth averaged magnitude of velocity comparison ');
%   xlabel('Yearday')
%   ylabel('Velocity (m/s')
%   axis ([223.6 223.8 0 .5])  
% %% 
% figure (11)
% hold on;  box on;
% plot(data.yd,data.fdav,'-b')
% plot(data.yd,data.fmag,'-r')
%  legend('data.dav =black','data.fdav =red')
%  title('Filtered depth averaged magnitude of velocities ');
%   xlabel('Yearday')
%   ylabel('Velocity (m/s')
%   axis ([223.6 223.8 0 .5])  
%COPIED FROM 2012 SCRIPT
%% Figure 1 Water Depth

%figure(1) 
    
 %hold on;  box on;
 %plot(data.yd,data.dep,'-b') 
 %title('Pressure');
 %xlabel('Yearday')
  %ylabel('Water Depth')
 % legend ('Observed')
 %axis ([223.6 223.8 0 5])
 

  
  %% Figure 2 Temperature

%figure(2) 
    
% hold on;  box on;
 %plot(data.yd,data.temp,'-b') 
%title('Temperature');
 % xlabel('Yearday')
 % ylabel('Temperature')
  %legend ('Observed')
  %axis ([223.6 223.8 15 30])
 
 

%% Figure 3 Non Filtered data Velocity

% figure(5)
% subplot(311)
%     hold on; box on;
%     imagesc(data.yd, data.dbe,data.v1(:,3:numberofbins)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([223.6 223.8 0 5])
%     ylabel('z (m)')
%     set(gca,'CLim',[-.2 .2])
%     title('Non Filtered data Velocity E component ');
%     colorbar
%     colormap(fig);
%     
%  subplot(312)
%     hold on; box on;
%     imagesc(data.yd, data.dbn,data.v2(:,3:numberofbins)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([223.6 223.8 0 5])
%     ylabel('z (m)')
%     set(gca,'CLim',[-.2 .2])
%     title('Non Filtered data Velocity N component ');
%     colorbar
%    colormap(fig);
%  
% subplot(313)
%     hold on; box on;
%     imagesc(data.yd, data.dbv,data.v3(:,3:numberofbins)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([223.6 223.8 0 5])
%     ylabel('z (m)')
%     set(gca,'CLim',[-.02 .02])
%     title('Non Filtered data Velocity V component ');
%     colorbar    
%     colormap(fig);
%  %  print -dpdf PCS_validation
% 
%  %% Figure 4 Filtered data Velocity
% 
% figure(6)
% subplot(311)
%     hold on; box on;
%     imagesc(data.yd, data.dbe,data.fv1(:,3:numberofbins)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([223.6 223.8 0 5])
%     ylabel('z (m)')
%     set(gca,'CLim',[-.2 .2])
%     title('Filtered data Velocity E component ');
%     colorbar
%     colormap(fig);
%     
%  subplot(312)
%     hold on; box on;
%     imagesc(data.yd, data.dbn,data.fv2(:,3:numberofbins)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([223.6 223.8 0 5])
%     ylabel('z (m)')
%     set(gca,'CLim',[-.2 .2])
%     title('N component ');
%     colorbar
%     colormap(colormap1);
%  
% subplot(313)
%     hold on; box on;
%     imagesc(data.yd, data.dbv,data.fv3(:,3:numberofbins)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([223.6 223.8 0 5])
%     ylabel('z (m)')
%     set(gca,'CLim',[-.02 .02])
%     title('V component ');
%     colorbar
%     colormap(colormap1);
 %  print -dpdf PCS_validation
%% figure 5 Amplitude

% figure(7)
% subplot(311)
%     hold on; box on;
%     imagesc(data.yd, data.dbe,data.a1(:,3:numberofbins)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([223.6 223.8 0 5])
%     ylabel('z (m)')
%     set(gca,'CLim',[0 200])
%     title('Amplitude E component ');
%     colorbar
%     colormap(fig);
%     
%  subplot(312)
%     hold on; box on;
%     imagesc(data.yd, data.dbn,data.a2(:,3:numberofbins)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([223.6 223.8 0 5])
%     ylabel('z (m)')
%     set(gca,'CLim',[0 200])
%     title('Amplitude N component ');
%     colorbar
%     colormap(fig);
%  
% subplot(313)
%     hold on; box on;
%     imagesc(data.yd, data.dbv,data.a3(:,3:numberofbins)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([223.6 223.8 0 5])
%     ylabel('z (m)')
%     set(gca,'CLim',[0 200])
%     title('Amplitude V component ');
%     colorbar
%     colormap(fig);
    
    
    %% figure 6 Correlation

figure(8)
subplot(311)
    hold on; box on;
    imagesc(data.yd, data.dbe,data.c1(:,3:numberofbins)')
    plot(data.yd,data.dep,'-k')
    shading interp
    set(gca,'YDir','normal')
    %set(gca,'YLim',[0 5])
    %set(gca,'XLim','tight')
    %axis ([223.6 223.8 0 5])
    ylabel('z (m)')
    set(gca,'CLim',[0 100])
    title('Correlation E component ');
    colorbar
    colormap(fig);
    
 subplot(312)
    hold on; box on;
    imagesc(data.yd, data.dbn,data.c2(:,3:numberofbins)')
    plot(data.yd,data.dep,'-k')
    shading interp
    set(gca,'YDir','normal')
    %set(gca,'YLim',[0 5])
    %set(gca,'XLim','tight')
    %axis ([223.6 223.8 0 5])
    ylabel('z (m)')
    set(gca,'CLim',[0 100])
    title('Correlation N component ');
    colorbar
    colormap(fig);
 
subplot(313)
    hold on; box on;
    imagesc(data.yd, data.dbv,data.c3(:,3:numberofbins)')
    plot(data.yd,data.dep,'-k')
    shading interp
    set(gca,'YDir','normal')
    %set(gca,'YLim',[0 5])
    %set(gca,'XLim','tight')
    %axis ([223.6 223.8 0 5])
    ylabel('z (m)')
    set(gca,'CLim',[0 100])
    title('Correlation V component ');
    colorbar 
    colormap(fig);
    
    %% figure 7 Magnitude of velocity
% 
% figure(9)
% 
%     hold on; box on;
%     imagesc(data.yd, data.dbe,data.mag(:,1:97)')
%     %imagesc(data.depthfile(:,3)')
%     plot(data.yd,data.dep,'-k')
%     shading interp
%     set(gca,'YDir','normal')
%     %set(gca,'YLim',[0 5])
%     %set(gca,'XLim','tight')
%     %axis ([223.6 223.8 0 5])
%     ylabel('z (m)')
%     set(gca,'CLim',[0 .2])
%     title('Magnitude of velocity');
%     colorbar
%     colormap(jet);
    
  
    
    %% figure 8 ploting Velocities
%  figure (10)   
%  hold on;  box on;
%  plot(data.yd,data.mag,'-b') 
%  title('Magnitude of velocities ');
%   xlabel('Yearday')
%   ylabel('Velocity (m/s')
%   %axis ([223.6 223.8 0 .5])
 
    %% figure 9 Direction

%figure(9)
%subplot(111)
%    hold on; box on;
    
    %imagesc(data.yd, data.dbv,data.c3(:,3:numberofbins)')
    %plot(data.yd,data.dep,'-k')
    
%    imagesc(data.yd, data.angle2d(:,1:97)',data.mag(:,1:97)')
%    plot(data.yd,data.dep,'-k')
%    shading interp
%    set(gca,'YDir','normal')
%    %set(gca,'YLim',[0 5])
%    %set(gca,'XLim','tight')
%    axis ([223.6 223.8 0 360])
%    ylabel('z (m)')
%    set(gca,'CLim',[0 .2])
%    title('E velocity ');
%    colorbar
%    colormap(fig);
  
%% figure 10

figure (11)
subplot (311)
hold on;  box on;
 plot(data.yd,data.dav,'-b') 
 title('Non filtered depth averaged magnitude of velocity ');
  xlabel('Yearday')
  ylabel('Velocity (m/s')
  %axis ([223.6 223.8 0 .5])
subplot (312)
hold on;  box on;
 plot(data.yd,data.fdav,'-r') 
 title('Filtered depth averaged magnitude of velocity ');
  xlabel('Yearday')
  ylabel('Velocity (m/s')
  %axis ([223.6 223.8 0 .5])  
subplot (313)
hold on;  box on;
 plot(data.yd,data.dav,'-b')
 plot(data.yd,data.fdav,'-r')
 legend('data.dav =black','data.fdav =red')
 title('Depth averaged magnitude of velocity comparison ');
  xlabel('Yearday')
  ylabel('Velocity (m/s')
  %axis ([223.6 223.8 0 .5])  
%% 
figure (12)
hold on;  box on;
 plot(data.yd,data.fdav,'-b')
 plot(data.yd,data.fmag,'-r')
 legend('data.dav =black','data.fdav =red')
 title('Filtered depth averaged magnitude of velocities ');
  xlabel('Yearday')
  ylabel('Velocity (m/s')
  %axis ([223.6 223.8 0 .5])  
   
   %   %% end
% disp('End'); 
% disp('   ');   
%   %% wind rose
%   
%   
% 
%   
%   
%   
%   
%   