
%% put data into matlab format
savename = TDT_MATformatBDH_KTB(2,'hi'); 
    %channel 1 or 2, 'hi'=6khz FP sampling, 'lo'=101hz. Use 'lo' if possible
    %KTB data is 'hi'
    %savename = FormatForFPpipe(sig,ref);  
    %if data already converted, feed to FP pipe by preloading in sig and ref .mat

savename=strcat(savename,'.mat');
load(savename)

%% debleach, calc df/f, subtract 405 from 472
debleach_flag=1; 
numExp =0;  %0 = cubic poly fit, 1=single exp fit,  2=dbl exp fit
startfit = 1*6000; % start exp fit at this min (*100hz & 60sec)
endfit =min(length(sig_472_RS),length(sig_405_RS)); %50*6000; %length(sig_472_RS);  % end fit
sig_472_RS = sig_472_RS(1:endfit);
sig_405_RS = sig_405_RS(1:endfit);
timeFP_RS = timeFP_RS(1:endfit);
timeFP_frames = timeFP_frames(1:endfit);
baseline_method=2   % 1=main peak of histogram for dataset 
                    % 2=median val
                    % 3=manually assigned baseline (set manualbase)
manualbase=102;%53
dfof=debleachBDH(timeFP_RS,sig_472_RS,'472',savename,debleach_flag,numExp,startfit,endfit,baseline_method,manualbase);
dfof_control=debleachBDH(timeFP_RS,sig_405_RS,'405',savename,debleach_flag,numExp,startfit,endfit,baseline_method,manualbase);

% subtract reference signal
subtract_method='Subtract';
dfofCorr=subtract_refBDH(timeFP_RS,dfof,dfof_control,savename,subtract_method);
    % subtract_method='MinResid'  use bestfit model
    % subtract_method='Subtract' just subtract them
    % subtract_method='None' no correction with 405
    % subtract_method= number, manually add this value
    
%% smooth, z-score
dfofCorr=filtfilt(ones(1,100)/100,1,dfofCorr);
%z score
ZdfofCorr = (dfofCorr - mean(dfofCorr(startfit:length(dfofCorr))))  / std(dfofCorr(startfit:length(dfofCorr)));

savename2=strcat(savename,'_dFFcorr');
save(savename2,'dfofCorr','ZdfofCorr');%,'camtime','timeFP_RS');

% ALT OPTION for z-scoring, e.g. track changes in global signal by 
% z-scoring a reference part of an FP trace (startZ, endZ), then find fraction of points
% >1SD above mean in a readout part of the trace (startRO, endRO)
% 100 Hz x 60 s = 6000 pts/min
%    startZ=15*6000; endZ=20*6000; stdevThresh = 1;
%    startRO=50*6000; endRO=60*6000;
    
%    ZdfofCorrPart = (dfofCorr - mean(dfofCorr(startZ:endZ)))  / std(dfofCorr(startZ:endZ));
    
%    ProcZpart1 = ZdfofCorrPart-stdevThresh;
 %   numptsOvThreshBase = find(ProcZpart1(startZ:endZ) >0);
%    numptsOvThreshReadout = find(ProcZpart1(startRO:endRO)>0);
    
 %   FractBase=length(numptsOvThreshBase)/(endZ-startZ);
%    FractRO = length(numptsOvThreshReadout)/(endRO-startRO);

%% scroll to bottom if using video without Synapse timestamps

%% new PSTH
%evTime is vector holding frames of interest
%Assume 100Hz resampling of FP data
%use dfofCorr or ZdfofCorr in FPpsth arguments
%timeFP_frames=timeFP_RS; %if using video not encoded by Synapse
    Pre=300;
    Post=300;
    BaseStart=290;%Pre*.9;
    BaseEnd=100;%Pre*.1;
    Skip=400;
    Norm =1; %=0 no normalization to baseline,   =1 normalize to baseline
[PSTHarrNew,markFP,markFPfr]=FPpsth(camtime,timeFP_RS,timeFP_frames,ZdfofCorr,evTime,Pre,Post,BaseStart,BaseEnd,Skip,Norm); %camtime,FP time, proc FP, events, pre, post, start baseline @ -, end baseline @ -

ylim([-0.14 .55]);

%plot FP vs time with events marked
%figure;plot(timeFP_RS,dfofCorr);hold on; plot(timeFP_RS,markFP);

%plot FP vs frame number with events marked
figure;plot(timeFP_frames,dfofCorr);hold on; plot(timeFP_frames,markFPfr);


