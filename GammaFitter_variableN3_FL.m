function D = GammaFitter_variableN3_FL(tracks,TimeInt,TruncLength)
% GAMMMAFITTER_VARIABLEN
% This code fits a gamma distribution to a set of 2D track data. You can
% choose the step length and size to truncate tracks to and it will determine
% the independent steps and fit the distribution. Make sure you set the
% parameters at the end of this section appropriately for your data and the
% fit you want to do.
% E.g.1; TimeInt=1, TruncLength=25 ; the average of 24 independent 1 steps
%                                    is used to represent the diffusion
%                                    coefficient
% E.g.2; TimeInt=4, TruncLength=25 ; the average of 6 independant 4 frame
%                                    steps is used
%
% INPUTS
% tracks      output from Andreas, will get turned into:(columns: x (pixels), y(pixels), frame number, trajectory
%             number)
% TimeInt     is the time interval to use. E.g for to use 1 frame steps '1'
% TruncLength is the length you want to truncate all the tracks to. Should
%             be at least 4*TimeInt+1 for decent statistics 
%
% OUTPUTS
% D           column of msd values
%
% EXAMPLE OF USE
% 1. I dragged and dropped the 'compileddata' you sent into the command window
% on Matlab. Then select 'Numeric Matrix' from OutputType dropdown menu in
% the import wizard. Click the green tick (Import selection).
% 2. Check (and adjust if needed) the parameters 'pixel', 'dT',
% 'DhistFitType', 'InitGuess' and 'rangeD' just beneath this code
% explanation. Ctrl+S to save your changes
% 3. Type (or copy-paste) the following to the command line:%
% D = GammaFitter_variableN3_FL(Input,4,20);
%

% To change Input format
% Input = vertcat(CompiledTracks{:});


% Helen Miller August 2021. With a little help from Matt Stracey's fitting
% functions

%%PARAMETERS TO SET FOR A DATA SET
pixel = 0.1169; % length per pixel in um (0.096 = 96nm) -> Fred: set to 1, since script 3 is already converting to um
dT = 0.35; % time per frame in seconds (1/frames per sec, should be the cycle time, not the exposure time)
% DhistFitType= 'two_gamma_1constrained'; %choose one of _1constrained
DhistFitType= 'two_free';
% DhistFitType= 'one_free';
% 'three_free', 'two_free', 'one_free', 'two_gamma_1constrained';
%'three_gamma_1constrained_fit'; 
% InitGuess = [0.5 0.0005 0.4 0.005 0.1 0.02]; % 0.5 0.03 0.25 0.1]; % for 1 species fit need [fraction1 D1];
InitGuess = [0.7 0.0005 0.3 0.005];
% InitGuess = [0.7 0.000478 0.3 0.005]; %for constrained fit
% InitGuess = [1 0.0005];
% InitGuess for twogamma_1constrained  [percent D(fixed) percent D(guess)]
% for 1 species fit need [fraction1 D1];
% for two species fit need [fraction1 D1 fraction2 D2])
% for three species fit need [fraction1 D1 fraction2 D2 fraction3 D3])
rangeD = [0 0.0005 0.02]; % D range for the histogram in um^2/s
%%
%quickly change the shape of the input 
trackstemp=tracks;
tracks=horzcat(trackstemp(:,2:3),trackstemp(:,1),trackstemp(:,5));
% tracks=horzcat(trackstemp(:,3:4),trackstemp(:,1),trackstemp(:,11));

%calculate the number of independent steps you will get from each
%trajectory
N=floor((TruncLength-1)./TimeInt);
%check the relative sizes of TimeInt and TruncLength won't cause too much
%uncertainty (due to too little averaging)
if (4*TimeInt+1)>TruncLength
    disp('Warning: you are truncating your tracks to be short relative to the step size');
    disp('This means there is less averaging of the stochastic nature of the diffusion,')
    disp('and degrades the accuracy of your measurements. TruncLength> 4*TimeInt+1 is recommended');
    disp('press Ctrl+C to exit, or any other key to continue');
    pause
else %do nothing, this is fine
end    

%determine the maximun track number
MaxTraj=max(tracks(:,4));
%%Truncate the tracks; chop up longer tracks
%make all tracks 25 frames long
for uu=1:max(tracks(:,4))
   rows=find(tracks(:,4)==uu);
   if length(rows)<TruncLength
       %it's too short, get rid of it
      tracks(rows,:)=[]; 
   elseif length(rows)==TruncLength
       %it's the right length, leave it be
   else %it's more than the truncated length
       %work out how many times the truncated length it is
       NoMinis=floor(length(rows)./TruncLength); 
       %renumber truncframes length chunks into new trajectories
       for vv=2:NoMinis
           MaxTraj=MaxTraj+1;
           tracks(rows(((vv-1)*TruncLength)+1:(vv*TruncLength)),4)=MaxTraj;       
       end
       tracks(rows((NoMinis*TruncLength)+1:end),:)=[];
   end 
end
%all tracks should now be truncframes long and have unique frame numbers

%Now check every trajectory is consecutive, delete it if it's not; if it
%is, calculate the MSDs you need
%preallocate variable for MSDs
MSD=-1*ones(max(tracks(:,4)),1);
for ww=1:max(tracks(:,4))
       rows2=find(tracks(:,4)==ww);
       rowsConseq=find((tracks(rows2(2:end),3)-tracks(rows2(1:end-1),3))==1);
       if length(rowsConseq)==TruncLength-1
           %it's consecutive, leave it be, but calculate the MSDS
            minitraj=tracks(rows2,:);
         %  MSD(ww,1)=sum((sum((minitraj((1+TimeInt):end,1:2) - minitraj(1:end-TimeInt,1:2)).^2,2)))./N; HM this would be if all steps used, but we want just the independent ones
            AllIntervals=(sum((minitraj((1+TimeInt):end,1:2) - minitraj(1:end-TimeInt,1:2)).^2,2));
            %now average only the independent ones
            %AllIntervals(1:TimeInt:end)
            MSD(ww,1)=mean(AllIntervals(1:TimeInt:end));
       else %it's not consecutive; delete it
           tracks(rows2,:)=[];
       end
end
%Now all the tracks are the right length and consecutive
%get rid of any MSD that correspond to nonconsecutive rows
for pp=length(MSD(:,1)):-1:1 %loop backwards to not mess up the numbering
   if MSD(pp)==-1
       MSD(pp)=[];
   end
end

%convert into um 
MSD = MSD* pixel^2; % convert from pixel to length units

% calculate D from MSD, assuming 2D
D = MSD/(4*TimeInt*dT);

%plot the histogram
binSpacing = rangeD(2)-rangeD(1);
bins=rangeD(1):rangeD(2):rangeD(3);
histDCount =histcounts(D,bins); %(and these should be plotted at rangeD+binSpacing/2
ncounts=sum(histDCount)*binSpacing;
histDCount = histDCount./(sum(histDCount)*binSpacing); %normalize area

Dhist_plot2 = figure;
axes2 = axes('Parent',Dhist_plot2,'LineWidth',3,'FontSize',16);
box(axes2,'off');
hold(axes2,'all');
bar1 = bar(bins(1:end-1)+(binSpacing./2),histDCount,'BarWidth',1,'FaceColor',[.7 .7 .7],'EdgeColor',[.4 .4 .4],'LineWidth',4);
baseline1 = get(bar1,'BaseLine');
set(baseline1,'LineWidth',3);
xlim([min(bins), max(bins)]);
xlabel('Measured diffusion coefficients in \mum^{2}s^{-1}','FontSize',16);
ylabel('Probability density','FontSize',16);
box on

%display values to screen
fprintf('\n')
disp(['Number of trajectories used (after chopping them up) = ' num2str(length(MSD))]);

DhistMinSteps = N; % minimum number of steps for a track to be analyzed
x = bins+binSpacing/2;
[nparam,Expected]=Dhist_fitting(histDCount,x',DhistMinSteps, InitGuess ,DhistFitType,D,binSpacing,bins);
rchi2=ncounts*(sum(((histDCount-Expected').^2)./Expected'))./(length(histDCount)-nparam);
disp(['Reduced chi-squared = ' num2str(round(rchi2,3))]);
               
end


function [nparam,Expected]=Dhist_fitting(histDCount,x,N, initGuess, fit_type,D,binSpacing, binsD)
%fitting D histograms
% x =  D histogram range
% N = number of independent steps

%check initGuess is the correct length
if strcmp(fit_type,'one_free')
    if length(initGuess) ~= 2
        error('2 initial values neeeded for a one species fit: [fraction1 D1]')
    end
elseif strcmp(fit_type,'two_free') || strcmp(fit_type,'two_gamma_2constrained') || strcmp(fit_type,'two_gamma_1constrained')
    if length(initGuess) ~= 4
        error('4 initial values neeeded for a two species fit:  [fraction1 D1 fraction2 D2]')
    end
elseif strcmp(fit_type,'three_free') || strcmp(fit_type,'three_gamma_3constrained') || strcmp(fit_type,'three_gamma_2constrained') || strcmp(fit_type,'three_gamma_1constrained')
    if length(initGuess) ~= 6
        error('6 initial values neeeded for a three species fit:  [fraction1 D1 fraction2 D2 fraction3 D3]')
    end
end

%set least sq curve fitting options
curvefitoptions = optimset( 'lsqcurvefit');
curvefitoptions = optimset( curvefitoptions, 'Display', 'off', 'MaxFunEvals', 100000, 'MaxIter', 100000);

%new x values array for smooth plotting
x2 = linspace(0,max(x),1000);

switch fit_type 
    case 'one_free'
        %upper and lower bounds
        lb = 0 ;
        ub = inf;
        initGuess = initGuess(2);
        
        %least squares curve fit
        [c resnorm] = lsqcurvefit(@(parameters,x) ...
            one_gamma_fit_no_amplitude(...
            parameters,x,N),...
            initGuess, x(1:end-1)', histDCount, lb, ub, curvefitoptions);
        
        %plot the fit
        hold on;
        cFit = one_gamma_fit_no_amplitude(c,x2,N);
       % plot(x2, cFit, 'b-','LineWidth',3)
        
        %start bootstrapping
        NewD=zeros(10,1);
        for ii=1:10
            %choose random 80% of data with replacement
            RandVector=randi(length(D),round(0.8*length(D)),1);
            Dboot=D(RandVector);
            %make xboot and histDcountboot
            histDCountboot =histcounts(Dboot,binsD); %HM (and these should be plotted at rangeD+binSpacing/2
            histDCountboot = histDCountboot./(sum(histDCountboot)*binSpacing);
            %fit it using mostly the same params
            [cboot, ~] = lsqcurvefit(@(parameters,x) one_gamma_fit_no_amplitude(parameters,x,N),initGuess, x(1:end-1)', histDCountboot, lb, ub, curvefitoptions);
            NewD(ii)=cboot;  
            clear Dboot RandVector jj histDCountboot
        end
        %calc the std of the NewD
        Derror=std(NewD,1);
        CFitRange(1,:)=one_gamma_fit_no_amplitude(c+Derror,x2,N);
        CFitRange(2,:)=one_gamma_fit_no_amplitude(c-Derror,x2,N);
        [hl1, hb1]=boundedline(x2,cFit, [ -1*(max([CFitRange(1:2,:); cFit(1,:)])-cFit(1,:)) ; min([CFitRange(1:2,:);cFit(1,:)])-cFit(1,:) ]','b-','alpha');
        hl1.LineWidth=3;
        %boundedline(x2,cFit, [ -1*(max([CFitRange(1:2,:); cFit(1,:)])-cFit(1,:)) ; zeros(1,size(cFit,2)) ]','r-','alpha');
        % end bootstrapping
        
        disp(['  D = ', num2str(c,3),' \pm ', num2str(Derror,3)]);
        title({'Single species fit'; [ 'D = ', num2str(c,3),' \pm ', num2str(Derror,2)] },'FontSize',14);
        set(gca, 'Layer', 'Top')
        
        Expected = one_gamma_fit_no_amplitude(c,x(1:end-1),N);
        nparam=1;
        
    case 'two_free'
        
        initGuess = initGuess([1 2 4]);
        lb = [0 0 0];
        ub = [inf inf inf];
        
        c = lsqcurvefit(@(parameters,x) ...
            two_gamma_fit(...
            parameters,x,N),...
            initGuess, x(1:end-1)', histDCount, lb, ub, curvefitoptions);
              
        cFit = two_gamma_fit(c,x2,N);
        hold on;
        plot(x2, cFit, 'k--','LineWidth',3)
        
        cFit1 = one_gamma_fit(c(1:2),x2,N);
        hold on;
        plot(x2, cFit1, 'r-','LineWidth',3)
        
        cFit2 = one_gamma_fit([1-c(1), c(3)],x2,N);
        hold on;
       plot(x2, cFit2, 'b-','LineWidth',3)
        
        
        %start bootstrapping
        NewD=zeros(10,3);
        for ii=1:10
            %choose random 80% of data with replacement
            RandVector=randi(length(D),round(0.8*length(D)),1);
            Dboot=D(RandVector);
            %make xboot and histDcountboot
            histDCountboot =histcounts(Dboot,binsD); %HM (and these should be plotted at rangeD+binSpacing/2
            histDCountboot = histDCountboot./(sum(histDCountboot)*binSpacing);
            %fit it using mostly the same params
            cboot = lsqcurvefit(@(parameters,x) two_gamma_fit(parameters,x,N),initGuess, x(1:end-1)', histDCountboot, lb, ub, curvefitoptions);
            NewD(ii,:)=cboot(:);  
            clear Dboot RandVector jj histDCountboot
        end
        %calc the std of the NewD
        Derror=std(NewD,1); %should calc std down columns
        %return this as an error
        CFitRange(1,:)=one_gamma_fit([c(1) c(2)+Derror(2)],x2,N); %make the max and min with fitted best amplitude
        CFitRange(2,:)=one_gamma_fit([c(1) c(2)-Derror(2)],x2,N);
        CFitRange(3,:)=one_gamma_fit([1-c(1), c(3)+Derror(3)],x2,N);
        CFitRange(4,:)=one_gamma_fit([1-c(1), c(3)-Derror(3)],x2,N); 
%         plot(x2,((c(1)+Derror(1))./c(1))* CFitRange(1,:), 'm-','LineWidth',3)
%         plot(x2, ((c(1)+Derror(1))./c(1))*CFitRange(2,:), 'c-','LineWidth',3)
%         plot(x2, ((c(1)-Derror(1))./c(1))*CFitRange(1,:), 'y-','LineWidth',3)
%         plot(x2, ((c(1)-Derror(1))./c(1))*CFitRange(2,:), 'g-','LineWidth',3)
%         plot(x2,((1-c(1)+Derror(1))./(1-c(1)))* CFitRange(3,:), 'm-','LineWidth',3)
%         plot(x2, ((1-c(1)+Derror(1))./(1-c(1)))*CFitRange(4,:), 'c-','LineWidth',3)
%         plot(x2, ((1-c(1)-Derror(1))./(1-c(1)))*CFitRange(3,:), 'y-','LineWidth',3)
%         plot(x2, ((1-c(1)-Derror(1))./(1-c(1)))*CFitRange(4,:), 'g-','LineWidth',3)
%        pause
        [hlk, hpk]=boundedline(x2,cFit, [ -1*((( ((c(1)+Derror(1))./c(1)) *max([CFitRange(1:2,:);cFit1(1,:)]))-cFit1(1,:))+((((1-c(1)+Derror(1))./(1-c(1)))*max([CFitRange(3:4,:);cFit2(1,:)]))-cFit2(1,:))) ; ((((c(1)-Derror(1))./c(1))*min([CFitRange(1:2,:);cFit1(1,:)]))-cFit1(1,:))+((((1-c(1)-Derror(1))./(1-c(1)))*min([CFitRange(3:4,:);cFit2(1,:)]))-cFit2(1,:)) ]','k--','alpha');
        [hlr, hpr]=boundedline(x2,cFit1, [ -1*(( ((c(1)+Derror(1))./c(1)) *max([CFitRange(1:2,:);cFit1(1,:)]))-cFit1(1,:)); (((c(1)-Derror(1))./c(1))*min([CFitRange(1:2,:);cFit1(1,:)]))-cFit1(1,:) ]','r-','alpha');
        [hlb, hpb]=boundedline(x2,cFit2, [ -1*((((1-c(1)+Derror(1))./(1-c(1)))*max([CFitRange(3:4,:);cFit2(1,:)]))-cFit2(1,:)) ;(((1-c(1)-Derror(1))./(1-c(1)))*min([CFitRange(3:4,:);cFit2(1,:)]))-cFit2(1,:) ]','b-','alpha');
        hlk.LineWidth=3;
        hlr.LineWidth=3;
        hlb.LineWidth=3;
        %boundedline(x2,cFit2, [ min(CFitRange(7:8,:))-cFit2(1,:) ; -1*(max(CFitRange(7:8,:))-cFit2(1,:)) ]','b-','alpha')
        % end bootstrapping
        
        frac1 = c(1)*100;
        frac2 = 100 - frac1;
        
        disp(['  D = ', num2str(c(2),3),' \pm ', num2str(Derror(1,2),3), '  percentage = ', num2str(frac1,3),' \pm ', num2str(100*Derror(1,1),3)]);
        disp(['  D = ', num2str(c(3),3),' \pm ', num2str(Derror(1,3),3), '  percentage = ', num2str(frac2,3),' \pm ', num2str(100*Derror(1,1),3)]);
                
        title({'Two species fit;', [num2str(frac1,3),' \pm ', num2str(100*Derror(1,1),2), '% at D = ', num2str(c(2),3),' \pm ', num2str(Derror(1,2),2)],[num2str(frac2,3),' \pm ', num2str(100*Derror(1,1),2), '% at D = ', num2str(c(3),3),' \pm ', num2str(Derror(1,3),2)] },'FontSize',14)
        set(gca, 'Layer', 'Top')
        
        Expected = two_gamma_fit(c,x(1:end-1),N);
        nparam=3;
        
    case 'two_gamma_1constrained'
        
        initGuess2 = [initGuess(1) initGuess(4)];
        constrained_fit_vals = initGuess(2);
        lb = [0 0];
        ub = [inf inf];
        
        c  = lsqcurvefit(@(parameters,x) ...
            two_gamma_1constrained_fit(...
            parameters,x,N,constrained_fit_vals),...
            initGuess2, x(1:end-1)', histDCount, lb, ub, curvefitoptions);
        
        cFit = two_gamma_1constrained_fit(c,x2,N, constrained_fit_vals);
        hold on;
        plot(x2, cFit, 'k--','LineWidth',3)
        
        cFit1 = one_gamma_fit([c(1),constrained_fit_vals(1)],x2,N);
        hold on;
        plot(x2, cFit1, 'r-','LineWidth',3)
        
        cFit2 = one_gamma_fit([(1-c(1)),c(2)],x2,N);
        hold on;
        plot(x2, cFit2, 'b-','LineWidth',3)
% HM 25/08/2022 need to add bootstrapping here     
        %start bootstrapping
        NewD=zeros(10,2);
        for ii=1:10
            %choose random 80% of data with replacement
            RandVector=randi(length(D),round(0.8*length(D)),1);
            Dboot=D(RandVector);
            %make xboot and histDcountboot
            histDCountboot =histcounts(Dboot,binsD); %HM (and these should be plotted at rangeD+binSpacing/2
            histDCountboot = histDCountboot./(sum(histDCountboot)*binSpacing);
            %fit it using mostly the same params
            cboot = lsqcurvefit(@(parameters,x) two_gamma_1constrained_fit(parameters,x,N,constrained_fit_vals),initGuess2, x(1:end-1)', histDCountboot, lb, ub, curvefitoptions);
            NewD(ii,:)=cboot(:);  
            clear Dboot RandVector jj histDCountboot
        end
        %calc the std of the NewD
        Derror=std(NewD,1); %should calc std down columns
        %return this as an error
        % 25/08/2022 up to here...
        CFitRange(1,:)=one_gamma_fit([1-c(1) c(2)+Derror(2)],x2,N); %make the max and min with fitted best amplitude
        CFitRange(2,:)=one_gamma_fit([1-c(1) c(2)-Derror(2)],x2,N);
        CFitRange(3,:)=one_gamma_fit([c(1), constrained_fit_vals],x2,N);
        CFitRange(4,:)=one_gamma_fit([c(1), constrained_fit_vals],x2,N); 
%         plot(x2,((c(1)+Derror(1))./c(1))* CFitRange(3,:), 'm-','LineWidth',3)
%         plot(x2, ((c(1)+Derror(1))./c(1))*CFitRange(4,:), 'c-','LineWidth',3) 
%         plot(x2, ((c(1)-Derror(1))./c(1))*CFitRange(3,:), 'y-','LineWidth',3)
%         plot(x2, ((c(1)-Derror(1))./c(1))*CFitRange(4,:), 'g-','LineWidth',3)
%         plot(x2,((1-c(1)+Derror(1))./(1-c(1)))*CFitRange(1,:), 'm-','LineWidth',3)
%         plot(x2, ((1-c(1)+Derror(1))./(1-c(1)))*CFitRange(2,:), 'c-','LineWidth',3)
%         plot(x2, ((1-c(1)-Derror(1))./(1-c(1)))*CFitRange(1,:), 'y-','LineWidth',3)
%         plot(x2, ((1-c(1)-Derror(1))./(1-c(1)))*CFitRange(2,:), 'g-','LineWidth',3)
                                                            
         [hlk, hpk]=boundedline(x2,cFit, [ -1*((( ((c(1)+Derror(1))./c(1)) *max([CFitRange(3:4,:);cFit1(1,:)]))-cFit1(1,:))+((((1-c(1)+Derror(1))./(1-c(1)))*max([CFitRange(1:2,:);cFit2(1,:)]))-cFit2(1,:))) ; ((((c(1)-Derror(1))./c(1))*min([CFitRange(3:4,:);cFit1(1,:)]))-cFit1(1,:))+((((1-c(1)-Derror(1))./(1-c(1)))*min([CFitRange(1:2,:);cFit2(1,:)]))-cFit2(1,:)) ]','k--','alpha');
         [hlr, hpr]=boundedline(x2,cFit1, [ -1*(( ((c(1)+Derror(1))./c(1)) *max([CFitRange(3:4,:);cFit1(1,:)]))-cFit1(1,:)); (((c(1)-Derror(1))./c(1))*min([CFitRange(3:4,:);cFit1(1,:)]))-cFit1(1,:) ]','r-','alpha');
         [hlb, hpb]=boundedline(x2,cFit2, [ -1*((((1-c(1)+Derror(1))./(1-c(1)))*max([CFitRange(1:2,:);cFit2(1,:)]))-cFit2(1,:)) ;(((1-c(1)-Derror(1))./(1-c(1)))*min([CFitRange(1:2,:);cFit2(1,:)]))-cFit2(1,:) ]','b-','alpha');
         hlk.LineWidth=3;
         hlr.LineWidth=3;
         hlb.LineWidth=3;
        %%boundedline(x2,cFit2, [ min(CFitRange(7:8,:))-cFit2(1,:) ; -1*(max(CFitRange(7:8,:))-cFit2(1,:)) ]','b-','alpha')
        % end bootstrapping
 
        frac1 = c(1)*100;
        frac2 = 100 - frac1;
        
        disp(['  D = ', num2str(constrained_fit_vals(1),2), '  percentage = ', num2str(frac1,3),' \pm ', num2str(100*Derror(1,1),2)]);
        disp(['  D = ', num2str(c(2),2),' \pm ', num2str(100*Derror(1,2),2), '  percentage = ', num2str(frac2,3),' \pm ', num2str(100*Derror(1,1),2)]);
        
        title({'Two species 1 constrained fit;', [num2str(frac1,3),' \pm ', num2str(100*Derror(1,1),2), '% at D = ', num2str(constrained_fit_vals(1),2)],[num2str(frac2,3), ' \pm ', num2str(100*Derror(1,1),2),'% at D = ', num2str(c(2),2),' \pm ', num2str(Derror(1,2),2),] })
        set(gca, 'Layer', 'Top')%HM edit 25/08/22   
        Expected = two_gamma_1constrained_fit(c,x(1:end-1),N,constrained_fit_vals);%HM edit 25/08/22
        nparam=2;%Hm edit 25/08/22
        
    case 'two_gamma_2constrained'
        
        initGuess2 = initGuess(1);
        constrained_fit_vals = [initGuess(2)  initGuess(4)];
        lb = 0 ;
        ub = inf;
        
        c = lsqcurvefit(@(parameters,x) ...
            two_gamma_2constrained_fit(...
            parameters,x,N,constrained_fit_vals),...
            initGuess2, x, histDCount, lb, ub, curvefitoptions);
        
        cFit = two_gamma_2constrained_fit(c,x2,N, constrained_fit_vals);
        hold on;
        plot(x2, cFit, 'k--','LineWidth',3)
        
        cFit1 = one_gamma_fit([c(1),constrained_fit_vals(1)],x2,N);
        hold on;
        plot(x2, cFit1, 'r-','LineWidth',3)
        
        cFit2 = one_gamma_fit([(1-c(1)),constrained_fit_vals(2)],x2,N);
        hold on;
        plot(x2, cFit2, 'b-','LineWidth',3)
        
        frac1 = c(1)*100;
        frac2 = 100 - frac1;
        
        disp(['  D = ',  num2str(constrained_fit_vals(1),3), '  percentage = ', num2str(frac1,3)]);
        disp(['  D = ',  num2str(constrained_fit_vals(2),3), '  percentage = ', num2str(frac2,3)]);
        
        title({'Two species 2 constrained fit;', [num2str(frac1,3), '% at D = ', num2str(constrained_fit_vals(1),3)],[num2str(frac2,3), '% at D = ', num2str(constrained_fit_vals(2),3)] });
        
        
        
    case 'three_free'
        
        initGuess = initGuess([1 2 3 4 6]);
        lb = [0 0 0 0 0 ];
        ub = [inf inf inf inf inf];
        
        c = lsqcurvefit(@(parameters,x) ...
            three_gamma_fit(...
            parameters,x,N),...
            initGuess, x(1:end-1)', histDCount, lb, ub, curvefitoptions);
        
        cFit = three_gamma_fit(c,x2,N);
        
        hold on;
        plot(x2, cFit, 'k--','LineWidth',3)
               
        cFit1 = one_gamma_fit([c(1),c(3)],x2,N);
        hold on;
        plot(x2, cFit1, 'r-','LineWidth',3)
        
        cFit2 = one_gamma_fit([c(2),c(4)],x2,N);
        hold on;
        plot(x2, cFit2, 'b-','LineWidth',3)
        
        cFit3 = one_gamma_fit([(1 - c(1) - c(2)),c(5)],x2,N);
        hold on;
        plot(x2, cFit3, 'g-','LineWidth',3)
        
        %start bootstrapping
        NewD=zeros(10,5);
        for ii=1:10
            %choose random 80% of data with replacement
            RandVector=randi(length(D),round(0.8*length(D)),1);
            Dboot=D(RandVector);
            %make xboot and histDcountboot
            histDCountboot =histcounts(Dboot,binsD); %HM (and these should be plotted at rangeD+binSpacing/2
            histDCountboot = histDCountboot./(sum(histDCountboot)*binSpacing);
            %fit it using mostly the same params
            cboot = lsqcurvefit(@(parameters,x) three_gamma_fit(parameters,x,N),initGuess, x(1:end-1)', histDCountboot, lb, ub, curvefitoptions);
            NewD(ii,:)=cboot(:);  
            clear Dboot RandVector jj histDCountboot
        end
        %calc the std of the NewD
        Derror=std(NewD,1); %should calc std down columns
        %return this as an error
        % end bootstrapping
        
        frac1 = c(1)*100;
        frac2 = c(2)*100;
        frac3 = 100 - (frac1 + frac2);
        frac3error=(Derror(1,1).^2+Derror(1,2).^2)^0.5;
        
        disp(['  D = ',  num2str(c(3),3), ' \pm ', num2str(Derror(1,3),3),'  percentage = ', num2str(frac1,3),' \pm ', num2str(100*Derror(1,1),3)]);
        disp(['  D = ',  num2str(c(4),3), ' \pm ', num2str(Derror(1,4),3), '  percentage = ', num2str(frac2,3),' \pm ', num2str(100*Derror(1,2),3)]);
        disp(['  D = ',  num2str(c(5),3), ' \pm ', num2str(Derror(1,5),3), '  percentage = ', num2str(frac3,3),' \pm ', num2str(100*frac3error,3)]);
        
        title({'Three species fit;', [num2str(frac1,3), ' \pm ', num2str(100*Derror(1,1),3),'% at D = ', num2str(c(3),2), ' \pm ', num2str(Derror(1,3),2)],...
            [num2str(frac2,3),' \pm ', num2str(100*Derror(1,2),2), '% at D = ', num2str(c(4),3),' \pm ', num2str(Derror(1,4),2)],...
            [num2str(frac3,3),' \pm ', num2str(100*frac3error,2), '% at D = ', num2str(c(5),3),' \pm ', num2str(Derror(1,5),2)] },'FontSize',14);
        
        Expected = three_gamma_fit(c,x(1:end-1),N);
        nparam=5;
        
            
    case 'three_gamma_1constrained'
        
        initGuess2 = [initGuess(1) initGuess(3) initGuess(4) initGuess(6)];
        constrained_fit_vals = (initGuess(2));
        
        ub = [inf inf inf inf];
        lb = [0 0 0 0];
        
        c = lsqcurvefit(@(parameters,x) ...
            three_gamma_1constrained_fit(...
            parameters,x,N,constrained_fit_vals),...
            initGuess2, x, histDCount, lb, ub, curvefitoptions);
        
        cFit = three_gamma_1constrained_fit(c,x2,N, constrained_fit_vals);
        hold on;
        plot(x2, cFit, 'k--','LineWidth',3)
        
        cFit1 = one_gamma_fit([c(1),constrained_fit_vals(1)],x2,N);
        hold on;
        plot(x2, cFit1, 'r-','LineWidth',3)
        
        cFit2 = one_gamma_fit([c(2),c(3)],x2,N);
        hold on;
        plot(x2, cFit2, '-', 'Color',[0,0.6,0.8], 'LineWidth',3)      
             
        cFit3 = one_gamma_fit([(1 - c(1) - c(2)),c(4)],x2,N);
        hold on;
        plot(x2, cFit3, 'g-','LineWidth',3)
        
        frac1 = c(1)*100;
        frac2 = c(2)*100;
        frac3 = 100 - (frac1 + frac2);
        
        disp(['  D = ',  num2str(constrained_fit_vals(1)),  num2str(frac1,3)]);
        disp(['  D = ',  num2str(c(3),3), '  percentage = ', num2str(frac2,3)]);
        disp(['  D = ',  num2str(c(4),3), '  percentage = ', num2str(frac3,3)]);

         title({'Three species 1 constrained fit;', [num2str(frac1,3), '% at D = ', num2str(constrained_fit_vals(1))],...
            [num2str(frac2,3), '% at D = ', num2str(c(3),3)],...
            [num2str(frac3,3), '% at D = ', num2str(c(4),3)] });
       
        case 'three_gamma_2constrained'
        
        ub = [inf inf inf ];
        lb = [0 0 0 ];
        
        initGuess2 = [initGuess(1) initGuess(3) initGuess(5)];
        constrained_fit_vals = [initGuess(2) initGuess(4)];
        
        c = lsqcurvefit(@(parameters,x) ...
            three_gamma_2constrained_fit(...
            parameters,x,N,constrained_fit_vals),...
            initGuess2, x, histDCount, lb, ub, curvefitoptions);
        
        
        cFit = three_gamma_2constrained_fit(c,x2,N, constrained_fit_vals);
        hold on;
        plot(x2, cFit, 'k--','LineWidth',3)
        
        cFit1 = one_gamma_fit([c(1),constrained_fit_vals(1)],x2,N);
        hold on;
        plot(x2, cFit1, 'r-','LineWidth',3)
        
        cFit2 = one_gamma_fit([c(2),constrained_fit_vals(2)],x2,N);
        hold on;
        plot(x2, cFit2, 'b-','LineWidth',3)
        
        cFit3 = one_gamma_fit([(1 - c(1) - c(2)),c(3)],x2,N);
        hold on;
        plot(x2, cFit3, 'g-','LineWidth',3)
        
              
        disp(['  D = ',  num2str(constrained_fit_vals(1)), '  percentage = ', num2str(c(1)*100),3]);
        disp(['  D = ',  num2str(constrained_fit_vals(2)), '  percentage = ', num2str(c(2)*100,3)]);
        disp(['  D = ',  num2str(c(3),3), '  percentage = ', num2str((1- c(1) - c(2)),3)]);
        
        title({'Three species 2 constrained fit;',...
            [num2str(c(1)*100,3), '% at D = ', num2str(constrained_fit_vals(1))],...
            [num2str(c(2)*100,3), '% at D = ', num2str(constrained_fit_vals(2))],...
            [num2str((1- c(1) - c(2))*100,3), '% at D = ', num2str(c(3),3)] })
       
        
    case 'three_gamma_3constrained'
        %three constrained fit
        
        initGuess2 = [initGuess(1) initGuess(3)];
        constrained_fit_vals = [initGuess(2) initGuess(4) initGuess(6)];
        
        ub = [inf inf ];
        lb = [0 0 ];
        
        c = lsqcurvefit(@(parameters,x) ...
            three_gamma_3constrained_fit(...
            parameters,x,N,constrained_fit_vals),...
            initGuess2, x, histDCount, lb, ub, curvefitoptions);
        
        cFit = three_gamma_3constrained_fit(c,x2,N, constrained_fit_vals);
        hold on;
        plot(x2, cFit, 'k--','LineWidth',3)
        
        cFit1 = one_gamma_fit([c(1),constrained_fit_vals(1)],x2,N);
        hold on;
        plot(x2, cFit1, 'r-','LineWidth',3)
        
        cFit2 = one_gamma_fit([c(2),constrained_fit_vals(2)],x2,N);
        hold on;
        %plot(x2, cFit2, 'b-','LineWidth',3)
        plot(x2, cFit2, '-', 'Color',[0,0.6,0.8], 'LineWidth',3)   
        
        cFit3 = one_gamma_fit([(1 - c(1) - c(2)),constrained_fit_vals(3)],x2,N);
        hold on;
        plot(x2, cFit3, 'g-','LineWidth',3)
        
        frac1 = c(1)*100;
        frac2 = c(2)*100;
        frac3 = 100 - (frac1 + frac2);
        
        disp(['  D = ',  num2str(constrained_fit_vals(1),3), '  percentage = ', num2str(frac1,3)]);
        disp(['  D = ',  num2str(constrained_fit_vals(2),3), '  percentage = ', num2str(frac2,3)]);
        disp(['  D = ',  num2str(constrained_fit_vals(3),3), '  percentage = ', num2str(frac3,3)]);
        
        title({'Three species 3 constrained fit;', [num2str(frac1,3), '% at D = ', num2str(constrained_fit_vals(1),3)],...
            [num2str(frac2,3), '% at D = ', num2str(constrained_fit_vals(2),3)],...
            [num2str(frac3,3), '% at D = ', num2str(constrained_fit_vals(3),3)] });
        
end
end

%different fitting equations (1,2,3 free and 2 and 3 constranined fit
%options)
function F = one_gamma_fit( parameters, x, n)

A(1) = parameters(1);
D(1) = parameters(2);

F= A(1)*(n/D(1))^n*x.^(n-1).*exp(-n*x/D(1))/factorial(n-1);

end

function F = one_gamma_fit_no_amplitude( D, x, n)
%equation for the pdf of a single D species
%no amplitude as area should = 1

F= (n/D)^n*x.^(n-1).*exp(-n*x/D)/factorial(n-1);


end

function F = two_gamma_fit( parameters, x, n)

A(1) = parameters(1);
D(1) = parameters(2);
D(2) = parameters(3);

F= A(1)*(n/D(1))^n*x.^(n-1).*exp(-n*x/D(1))/factorial(n-1) + ...
    (1 - A(1))*(n/D(2))^n*x.^(n-1).*exp(-n*x/D(2))/factorial(n-1);


end
function F = two_gamma_1constrained_fit( parameters, x, n, b)

A(1) = parameters(1);
D(1) = parameters(2);

F= A(1)*(n/b)^n*x.^(n-1).*exp(-n*x/b)/factorial(n-1) + ...
    (1 - A(1))*(n/D(1))^n*x.^(n-1).*exp(-n*x/D(1))/factorial(n-1);


end
function F = two_gamma_2constrained_fit( parameters, x, n, b)

%free amplitudes:
A(1) = parameters(1);

F = A(1)*(n/b(1))^n*x.^(n-1).*exp(-n*x/b(1))/factorial(n-1) + ...
    (1-A(1))*(n/b(2))^n*x.^(n-1).*exp(-n*x/b(2))/factorial(n-1);

end
function F = three_gamma_fit( parameters, x, n)

A(1) = parameters(1);
A(2) = parameters(2);
D1 = parameters(3);
D2 = parameters(4);
D3 = parameters(5);

F = A(1)*(n/D1)^n*x.^(n-1).*exp(-n*x/D1)/factorial(n-1) + ...
    A(2)*(n/D2)^n*x.^(n-1).*exp(-n*x/D2)/factorial(n-1) + ...
    (1-A(1)-A(2))*(n/D3)^n*x.^(n-1).*exp(-n*x/D3)/factorial(n-1);

end
function F = three_gamma_1constrained_fit( parameters, x, n,b)

b1= b(1);
A(1) = parameters(1);
A(2) = parameters(2);
D1 = parameters(3);
D2 = parameters(4);

F = A(1)*(n/b1)^n*x.^(n-1).*exp(-n*x/b1)/factorial(n-1) + ...
    A(2)*(n/D1)^n*x.^(n-1).*exp(-n*x/D1)/factorial(n-1) + ...
    (1-A(1)-A(2))*(n/D2)^n*x.^(n-1).*exp(-n*x/D2)/factorial(n-1);

end
function F = three_gamma_2constrained_fit( parameters, x, n, b)

b1= b(1);
b2= b(2);
A(1) = parameters(1);
A(2) = parameters(2);
D = parameters(3);

F = A(1)*(n/b1)^n*x.^(n-1).*exp(-n*x/b1)/factorial(n-1) + ...
    A(2)*(n/b2)^n*x.^(n-1).*exp(-n*x/b2)/factorial(n-1) + ...
    (1 - A(1) - A(2))*(n/D)^n*x.^(n-1).*exp(-n*x/D)/factorial(n-1);

end
function F = three_gamma_3constrained_fit( parameters, x, n, b)

b1= b(1);
b2= b(2);
b3= b(3);
A(1) = parameters(1);
A(2) = parameters(2);

%  F = A(1)*(n/b1)^n*x.^(n-1).*exp(-n*x/b1)/factorial(n-1) + ...
%     A(2)*(n/b3)^n*x.^(n-1).*exp(-n*x/b3)/factorial(n-1) + ...
%     (1 - A(1) - A(2))*(n/b2)^n*x.^(n-1).*exp(-n*x/b2)/factorial(n-1);


F = A(1)*(n/b1)^n*x.^(n-1).*exp(-n*x/b1)/factorial(n-1) + ...
    A(2)*(n/b2)^n*x.^(n-1).*exp(-n*x/b2)/factorial(n-1) + ...
    (1 - A(1) - A(2))*(n/b3)^n*x.^(n-1).*exp(-n*x/b3)/factorial(n-1);

end