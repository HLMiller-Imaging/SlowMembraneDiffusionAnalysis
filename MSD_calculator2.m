function [MSD] = MSD_calculator2(tracks,TruncLength,apply_filter)
% MSD_calculator
% Calculates MSDs
% It also makes MSD plots for a set of 2D tracks data. You can
% choose the size to truncate tracks to, and it finds intervals of 
% consecutive steps of that length. 
% Make sure you set the  parameters at the end of this section 
% appropriately for your data.
%
% INPUTS
% tracks      output from Andreas, will get turned into:(columns: x (pixels), y(pixels), frame number, trajectory
%             number)
% TruncLength is the length you want to truncate all the tracks to. Should
%             be at least 4*TimeInt+1 for decent statistics 
% apply_filter 1 to apply a filter at MSD(1)=1 to get rid of free halotag before
%             optimising the number of points to fit
%
% OUTPUTS
% MSD         Each row is the mean MSD values for a trajectory. First
%             column is one step, increasing up to TruncLength-1
%
% EXAMPLE OF USE
% 1. I dragged and dropped the 'compileddata' you sent into the command window
% on Matlab. Then select 'Numeric Matrix' from OutputType dropdown menu in
% the import wizard. Click the green tick (Import selection).
% 2. Check (and adjust if needed) the parameters 'pixel', 'dT',
% Ctrl+S to save your changes
% 3. Type (or copy-paste) the following to the command line:%
% [MSD] = MSD_calculator(compiledData,25,1);
%
% Helen Miller Sept 2021

%%PARAMETERS TO SET FOR A DATA SET
pixel =0.117; % length per pixel in um (0.096 = 96nm)
dT = 0.08; % time per frame in seconds (1/frames per sec, should be the cycle time, not the exposure time)

%%
%quickly change the shape of the input 
trackstemp=tracks;
tracks=horzcat(trackstemp(:,3:4),trackstemp(:,2),trackstemp(:,11)-min(trackstemp(:,11))+1);

%determine the maximun track number
MaxTraj=max(tracks(:,4));
%%Truncate the tracks; chop up longer tracks
%make all tracks TruncLength frames long
%loop over track numbers
Excess=0;
for uu=1:max(tracks(:,4))
   rows=find(tracks(:,4)==uu);
   minitraj2=tracks(rows,:);
   count3r=1;
   flag=0;
   if length(rows)<TruncLength
       %it's too short, set traj no to zero
      tracks(rows,4)=0; 
   elseif length(rows)==TruncLength
       %it's the right length, see if it's consecutive
       rowsConseq=find((tracks(rows(2:end),3)-tracks(rows(1:end-1),3))==1);
       if length(rowsConseq)==TruncLength-1
           %it's consecutive, leave it be
       else %it's not consecutive; delete it
           tracks(rows,4)=0;
       end  
   else %it's more than the truncated length
       %loop through testing every set of 25 frames to find a consecutive
       %one
       qq=1;
       count3r=0;
       flag=1;
       while qq<(length(rows)-TruncLength+1)
           rowsConseq=find((minitraj2(qq+1:qq+TruncLength-1,3)-minitraj2(qq:qq+TruncLength-2,3))==1);
           if length(rowsConseq)==TruncLength-1
                %it's consecutive; 
                if count3r==0 %first consec interval
                    %skip to next interval
                    qq=qq+TruncLength-1;
                    count3r=count3r+1;
                    flag=0;
                else
                    tracks(rows(qq:qq+TruncLength-1),4)=MaxTraj+count3r+Excess;%create a new traj no for this new consec. interval
                    count3r=count3r+1;
                    qq=qq+TruncLength-1;
                end
           else %if it's not consecutive, set the first traj no to zero
               tracks(rows(qq),4)=0;               
           end 
           qq=qq+1;
       end
       %set any remaining track nos to zero
       tracks(rows(qq:end),4)=0;
   end 
   Excess=Excess+count3r-1+flag;
end
%now delete all the rows where you set the traj no to zero
tracks(tracks(:,4)==0,:)=[];

%Now calculate the MSDs you need
%preallocate variable for MSDs
MSD=-1*ones(max(tracks(:,4)),TruncLength-1);
for ww=1:max(tracks(:,4))
       rows2=find(tracks(:,4)==ww);
        minitraj=tracks(rows2,:);
        for TimeInt=1:TruncLength-1
            %calculate the number of steps you will get from each trajectory
            N=TruncLength-TimeInt;
%                 minitraj((1+TimeInt):end,1:2) - minitraj(1:end-TimeInt,1:2)
%                 sum((minitraj((1+TimeInt):end,1:2) - minitraj(1:end-TimeInt,1:2)).^2,2)
%                 sum((sum((minitraj((1+TimeInt):end,1:2) - minitraj(1:end-TimeInt,1:2)).^2,2)))./N
            MSD(ww,TimeInt)=sum((sum((minitraj((1+TimeInt):end,1:2) - minitraj(1:end-TimeInt,1:2)).^2,2)))./N; %HM this would be if all steps used, but we want just the independent ones
        end
        %AllIntervals=(sum((minitraj((1+TimeInt):end,1:2) - minitraj(1:end-TimeInt,1:2)).^2,2));
        %now average only the independent ones
        %AllIntervals(1:TimeInt:end)
        %MSD(ww,1)=mean(AllIntervals(1:TimeInt:end));
end

%Now all the tracks are the right length and consecutive
%get rid of any MSD that correspond to nonconsecutive rows
for pp=length(MSD(:,1)):-1:1 %loop backwards to not mess up the numbering
   if MSD(pp,1)<=0
       MSD(pp,:)=[];
   end
end

%convert into um
MSD = MSD* pixel^2; % convert from pixel to length units
CondMean=mean(MSD,1);

% calculate D from MSD, assuming 2D
% D = MSD./(4*TimeLags*dT);
%temporary testing filter
if apply_filter==1
MSD(MSD(:,1)>1,:)=[];
end

%Make an MSDs plot
MSDsPlot = figure;
axes2 = axes('Parent',MSDsPlot,'LineWidth',3,'FontSize',16);
hold(axes2,'all');
TimeLags=dT*(1:1:TruncLength-3);
for ii=1:length(MSD(:,1))
    plot(TimeLags,MSD(ii,1:length(TimeLags(1,:))),'Color',[0.6 0.6 0.6]);
end
CondMean=mean(MSD,1);
plot(TimeLags,CondMean(1:length(TimeLags(1,:))),'r');

xlim([0, max(TimeLags(1,:))]);
xlabel('Time Lag (seconds)','FontSize',16);
ylabel('MSD \mum^2','FontSize',16);
ylim([0, 1.05*max(max(MSD(:,1:length(TimeLags(1,:)))))]);
set(axes2,'FontSize',16,'LineWidth',3,'YMinorTick','on','YScale','log');
box on

%display values to screen
disp(['Number of trajectories used (after truncation) = ' num2str(length(MSD))]);
               
end

