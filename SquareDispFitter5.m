function [FittedParams, rchi2,Percents,DiffCoeff] = SquareDispFitter5(MSD,NUse,Range)
% SquareDispFitter5
% This code fits the distribution of squared displacements for a the first 4 time lag
% with NUse distributions.  Make sure you set the 
% parameters at the end of this section appropriately for your data.
% The fit and residuals plot is best displayed on a log scale on x but this is
% difficult to do programatically. Ideally open the editor and do this. In
% this graph the time lags are plotted sequentially, so the first line is 1
% step, second is 2 step etc.
% In these fits there must be at least 5% of the data in each fitted population.
% chooses which displacements go together based on percentages - if you
% notice that grouping by species in the MSD vs T plots seems wrong please
% screenshot and send to me. The part of the code that decides which
% MSDs at a time point go with which at the next time point is not
% sophisticated but works in all samples tried so far.
%
% INPUTS
% MSD         the msd values from the linear fitting (MSD_vs_Lag4)in um^2
% NUse        The number of distribuions to use in this fitting. Possible
%             values 1-4
% Range       to fit eg [10e-5 2e-2]. You will be told how much of your
%             data is excluded by the range you have chosen. You probably
%             want to include all of it, although missing a few percent
%             might be ok if there are outliers.
%
% OUTPUTS
% fitresult   each row is a timepoint of frac1,msd1,frac2,msd2 etc for the
%             graph with the diffusion coefficients are calculated from. Probably not needed for most uses. 
% rchi2       Sum of the reduced chi^2 for all four fits
% percents    the percentage of molecules at each diffusion coefficient.
%             Columns correspond to columns of DiffCoeff
% DiffCoeff   the fitted diffusion coefficients found in um^2/s.
%
% e.g. get Percents =[91.4030  8.5970] DiffCoeff =[0.0003    0.0086] means
% 91.4% at 0.0003um^2/s, 8.6% at 0.0086um^2/s
%
% EXAMPLECODE
% [~, rchi2,Percents,DiffCoeff] = SquareDispFitter4(MSD,2,[10^(-4) 10^(-2)])
%
% Helen Miller September 2021.

%%PARAMETERS TO SET FOR A DATA SET
pixel =0.117; % length per pixel in um (0.096 = 96nm)
dT = 0.35; % time per frame in seconds (1/frames per sec, should be the 
%cycle time, not the exposure time)
Npoints=100; %number of points in the fit; too high and number of params 
%will make no difference, too low and fit will be bad. 100 is recommended
%%

%initialise variables
rchi2=zeros(1,4);
FittedParams=zeros(4,2*NUse); %variable to put percentages and D into
%loop over the different columns of MSD
for jj=1:4
    %Make a plot of the cumulative distribution function of the MSD for values
    %within the fit range
    rows=find(MSD(:,jj)<Range(2));
    temp=MSD(rows,jj);
    rows2=find(temp>Range(1));
    temp2=temp(rows2);
    %tell the user how many data values were excluded
    percentremoved=100*(length(MSD(:,jj))-length(rows2))./length(MSD(:,jj));
    disp(strcat(num2str(percentremoved),'% of data removed by your range limits'))
    %take the correct column of the distribution, order it according to sizeU
    MSDCol=temp2;
    MSDOrder=sort(MSDCol); %sort it into order
    NRows=length(MSDCol(:,1)); %ready for normalisation
    %split the log curve into 100 points equally spaced between highest and lowest
    CumulativeN=zeros(1+Npoints,1); %initialise variable for cumulative sum
    xvals=zeros(1+Npoints,1);
    interval=(log10(Range(2))-log10(Range(1)))./Npoints;
    for ii=1:Npoints
       xvals(ii+1)=10^(log10(Range(1))+ii*interval);
       rowsinc=find(MSDOrder<=xvals(ii+1));
       CumulativeN(ii+1)=length(rowsinc(:,1))./NRows; %normalised
    end
    %Now we need to fit it
    [xData, yData] = prepareCurveData( xvals, CumulativeN );
    [fitresult,Expected,nparam]=FitIt(xData,yData,NUse);
    %sort the outputs of the fit result into a variable based on the number
    %of distributions fitted
    if NUse==1
       FittedParams(jj,:)=[100 fitresult.a];
    elseif NUse==2
        %determine if b or c was larger
        if fitresult.a>=0.5
            FittedParams(jj,:)=[100*fitresult.a fitresult.b 100*(1-fitresult.a) fitresult.c]; 
        elseif fitresult.a<0.5
            FittedParams(jj,:)=[100*(1-fitresult.a) fitresult.c 100*fitresult.a fitresult.b]; 
        end
    elseif NUse==3
        %create a vector of the MSDs found
        testvector=[fitresult.b fitresult.d fitresult.g];
        %and of the percentages
        percentvector=100*[fitresult.a fitresult.c 1-fitresult.a-fitresult.c];
        %sort the MSDS into the right order
        [sorted,indices]=sort(percentvector);
        %sort the D using the same indices
        Dsorted=testvector(indices);
        %create the output vector
        FittedParams(jj,:)=[sorted(1,1) Dsorted(1,1) sorted(1,2) Dsorted(1,2) sorted(1,3) Dsorted(1,3)]; 
    elseif NUse==4
        %create a vector of the MSDs found
        testvector=[fitresult.b fitresult.d fitresult.g fitresult.k];
        %and of the percentages
        percentvector=100*[fitresult.a fitresult.c fitresult.f 1-fitresult.a-fitresult.c-fitresult.f];
        %sort the MSDS into the right order
        [sorted,indices]=sort(percentvector);
        %sort the percentages using the same indices
        Dsorted=testvector(indices);
        %create the output vector
        FittedParams(jj,:)=[sorted(1,1) Dsorted(1,1) sorted(1,2) Dsorted(1,2) sorted(1,3) Dsorted(1,3) sorted(1,4) Dsorted(1,4)]; 
    end
    rchi2(1,jj)=sum(((yData-Expected).^2))./(Npoints+1-nparam);
    if jj==1
        %create the figure and put hold on
        % Plot fit with data.
        figure( 'Name', 'Fit and residuals' );
        h = plot( fitresult, xData, yData ); hold on
        subplot( 3, 1, 1:2 );
        h = plot( fitresult, xData, yData ); hold on
        grid on
        % Plot residuals.
        subplot( 3, 1, 3 );
        h = plot( fitresult, xData, yData, 'residuals' ); hold on
        grid on
    else
        %just add to the figure
        subplot( 3, 1, 1:2 );
        h = plot( fitresult, xData, yData ); hold on         
        xlabel('Squared Displacement (\mum^2)');
        ylabel('Cumulative probability (P(r^2,t_{lag})');
        subplot( 3, 1, 3 );
        h = plot( fitresult, xData, yData, 'residuals' );
        xlabel('Squared Displacement (\mum^2)'); ylabel('Residuals');
    end
    %legend( h, 't=1', 'Fit to t=1','t=2', 'Fit to t=2','t=3', 'Fit to t=3','t=4', 'Fit to t=4', 'Location', 'NorthEast', 'Interpreter', 'none' );
end
%sum of r^2 (across all the fits)
disp(strcat('Total residual is ',num2str(sum(rchi2))));

%Now calculate the diffusion coefficients by fitting linearly
% Set up fittype and options.
fit1 = fittype( 'poly1' );
% Fit model to data fro each one
DiffCoeff=zeros(1,NUse);
%summary plot of the results
colour=[0 0 1;1 0 0; 0 0 0;0.929411768913269 0.694117665290833 0.125490203499794];
figure;
for zz=1:NUse
    % Display bubble chart with axis labels and legend
    [xData, yData] = prepareCurveData( dT*(1:4), FittedParams(:,2*zz)' );
    [fitresult1, ~] = fit( xData, yData, fit1 );
    DiffCoeff(1,zz)=(fitresult1.p1)./4; %I think the gradient should be =4D
    scatter(dT*(1:4),FittedParams(:,zz*2),10*FittedParams(:,(zz*2)-1),colour(zz,:),'filled'); hold on
    plot(fitresult1);
    %plot(1:4,FittedParams(:,zz*2),'Color',colour(zz,:));hold on
end
xlabel('Time lag(seconds)')
ylabel('MSD (\mum^2s^-^1)')

%mean %
Percents=mean(FittedParams(:,1:2:end-1)); %mean percentages
end

function [fitresult,Expected,nparam]=FitIt(xData,yData,NoDist)
% Set up fittype with NoDist number of populations and options.
if NoDist==4
    ft = fittype( '1-(a*exp(-x/b)+c*exp(-x/(d))+f*exp(-x/(g))+(1-a-c-f)*exp(-x/(k)))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [0.05 0 0.05 0 0.05 0  0];
    opts.Upper = [0.85 inf 0.85 inf 0.85 inf  inf];
    opts.StartPoint = [0.2 0.0001 0.2 0.001 0.2 0.01 0.1];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    Expected=1-(fitresult.a*exp(-xData/fitresult.b)+fitresult.c*exp(-xData/fitresult.d)+fitresult.f*exp(-xData/fitresult.g)+(1-fitresult.f-fitresult.a-fitresult.c)*exp(-xData/fitresult.k));
    nparam=7;
elseif NoDist==3
    ft = fittype( '1-(a*exp(-x/b)+c*exp(-x/(d))+(1-a-c)*exp(-x/(g)))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [0.05 0 0.05 0  0];
    opts.Upper = [0.9 inf 0.9 inf  inf];
    opts.StartPoint = [0.3 0.0001 0.3 0.001  0.01];
   
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    Expected=1-(fitresult.a*exp(-xData/fitresult.b)+fitresult.c*exp(-xData/fitresult.d)+(1-fitresult.c-fitresult.a)*exp(-xData/fitresult.g));
    nparam=5;
elseif NoDist==2
    ft = fittype( '1-(a*exp(-x/b)+(1-a)*exp(-x/(c)))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [0.05 0 0];
    opts.Upper = [0.95 inf inf];
    opts.StartPoint = [0.1 0.0001 0.001];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    Expected=1-(fitresult.a*exp(-xData/fitresult.b)+(1-fitresult.a)*exp(-xData/fitresult.c));
    nparam=3;
elseif NoDist==1
    ft = fittype( '1-exp(-x/a)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [0];
    opts.Upper = [inf];
    opts.StartPoint = [0.0001];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    Expected=1-exp(-xData/fitresult.a);
    nparam=1;
end
end
