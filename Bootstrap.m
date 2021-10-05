function [fitresult, rchi2,Percents,DiffCoeff,Errors]=Bootstrap(Nsample, Data, range,DisttoFit)
% Bootstrap
% This code runs a bootstrap (take a sample with replacement of the data
% 100 times; do the fitting,standard deviation of the fitted parameter
% from the samples is the standard error on the real fitted parameter
%
% INPUTS
% Nsample     number of samples to take; 25-200 recommended by Efron and
%             Tibshirani. Some references suggest little added accuracy for
%             >100 samples. HM uses 100.
% Data        the msd values from the linear fitting (MSD_vs_Lag4)in um^2
%             (MSDs in SquareDispFitter4)
% DisttoFit   The number of distribuions to use in this fitting. Possible
%             values 1-4 (NUse in SquareDispFitter4)
% Range       to fit eg [10e-5 2e-2]. You will be told how much of your
%             data is excluded by the range you have chosen. You probably
%             want to include all of it, although missing a few percent
%             might be ok if there are outliers. (same as
%             SquareDispFitter4)
%
% OUTPUTS
% fitresult   each row is a timepoint of frac1,msd1,frac2,msd2 etc for the
%             graph with the diffusion coefficients are calculated from. Probably not needed for most uses. 
% rchi2       Sum of the reduced chi^2 for all four fits
% percents    the percentage of molecules at each diffusion coefficient.
%             Columns correspond to columns of DiffCoeff
% DiffCoeff   the fitted diffusion coefficients found in um^2/s.
% Errors      on the percentages and then diffusion coefficients found. In
%             the same units as the returned percentages and diffusion
%             coefficients
%
% e.g. get Percents =[91.4030  8.5970] DiffCoeff =[0.0003    0.0086] Errors=[0.4 0.4 0.00005 0.0007]means
% 91.4 +/- 0.4 % at 0.0003 +/- 0.00005 um^2/s, 8.6% at 0.0086 +/- 0.0007um^2/s
%
% EXAMPLECODE
% [~, rchi2,Percents,DiffCoeff,Errors]=Bootstrap(100, MSD_Fixed_25_filt, [10^(-5) 10^(-1)],1)
%
% Helen Miller October 2021.

%initialise variable
Outputs=zeros(Nsample,2*DisttoFit);
NumSamples=length(Data(:,1));
for n=1:Nsample
    %make a subset via bootstrapping (same length as input data, random
    %sample with replacement
    X = ceil(NumSamples*rand(NumSamples,1));
    [~, ~,Percentstemp,DiffCoefftemp] = SquareDispFitter5(Data(X,:),DisttoFit,range);
    Outputs(n,:)=[Percentstemp,DiffCoefftemp];
    clear Percentstemp DiffCoefftemp
end

%Now calculate the standard deviations of these resamples (do it manually
%because need SAMPLE STD here and Matlab only has population built in )
Means=mean(Outputs,1);
Diffs=Outputs-Means;
Sums=sum(Diffs.^2,1);
Stds=sqrt(Sums./(Nsample-1));
Errors=Stds;
%Do full data fit
[fitresult, rchi2,Percents,DiffCoeff] = SquareDispFitter5(Data,DisttoFit,range);

end
