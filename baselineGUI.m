%filepath=uigetdir;
%% put data into matlab format
savename = TDT_MATformatBDH_KTB (1,'hi'); 
%I just added ,'' as an extra argument. May have to delete this.
    %channel 1 or 2 (1 is A group, 2 is B group), 'hi'=6khz FP sampling, 'lo'=101hz. Use 'lo' if possible
    
    %savename = FormatForFPpipe(sig,ref);  
    %if data already cstrsplitonverted, feed to FP pipe by preloading in sig and ref .mat

prompt = {'Starting frame (X 6000): ', 'Ending frame (X 6000): ', 'Manual baseline (automatic if empty): '};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1', '30', ''};
answer = inputdlg(prompt,dlgtitle,dims,definput)

startFrameMs = (str2double(answer(1)) * 6000) 
endFrameMs = (str2double(answer(2)) * 6000)
baseLineMethodInput = 2
manualBaselineInput = 0
if (answer(3) ~= "")
  manualBaselineInput = str2double(answer(3))
  baseLineMethodInput = 3
end  

savename=strcat(savename,'.mat');
load(savename)

%% debleach, calc df/f, subtract 405 from 472
debleach_flag=1; 
numExp =1;  %0 = cubic poly fit, 1=single exp fit,  2=dbl exp fit 3=cam method
startfit = startFrameMs; % start exp fit at this min (*100hz & 60sec)
endfit = endFrameMs; %length(sig_472_RS); %60*6000%length(sig_472_RS); %50*6000; %length(sig_472_RS);  % end fit
baseline_method=baseLineMethodInput   % 1=main peak of histogram for dataset 
                    % 2=median val
                    % 3=manually assigned baseline (set manualbase)
manualbase=manualBaselineInput;%53
dfof=debleachBDH(timeFP_RS,sig_472_RS,'472',savename,debleach_flag,numExp,startfit,endfit,baseline_method,manualbase);
dfof_control=debleachBDH(timeFP_RS,sig_405_RS,'405',savename,debleach_flag,numExp,startfit,endfit,baseline_method,manualbase);

% subtract reference signal
subtract_method='Subtract';

% make sure dfof and dfof_control have the same length

dfof_size=size(dfof);
dfof_size=max(dfof_size);
dfof_control_size=size(dfof_control);
dfof_control_size=max(dfof_control_size);

if(dfof_size < dfof_control_size)
    dfof_control=dfof_control(1:length(dfof));
else
    dfof=dfof(1:length(dfof_control));
end
timeFP_RS=timeFP_RS(1:length(dfof));

dfofCorr=subtract_refBDH(timeFP_RS,dfof,dfof_control,savename,subtract_method);
    % subtract_method='MinResid'  use bestfit model
    % subtract_method='Subtract' just subtract them
    % subtract_method='None' no correction with 405
    % subtract_method= number, manually add this value

    %% alternate 405 correction, scale 405 to match 472
     %   raw472=sig_472_RS;
     %   scaled_405(1:length(raw472))=0;
     %   scaled_405=controlFit(sig_472_RS,sig_405_RS);
       % scaled_405=scaled_405';
     %   dfofCorrScale=(sig_472_RS-scaled_405)/mean(scaled_405);
     %   subtract_method=0.008;
     %   dfofCorrScale=subtract_refBDH(timeFP_RS,dfofCorrScale',scaled_405,savename,subtract_method);
     %   dfofCorr=dfofCorrScale; % depending correction looks better
        %figure;hold on; plot(dfofCorrScale); plot(dfofCorr); plot(dfof)

%% close all the figures
% close all;
%% smooth, z-score
dfofCorr=filtfilt(ones(1,100)/100,1,dfofCorr);
%z score
ZdfofCorr = (dfofCorr - mean(dfofCorr(startfit:length(dfofCorr))))  / std(dfofCorr(startfit:length(dfofCorr)));
%z score using part of trace; 100 Hz x 60 s = 6000 pts/min
    %startZ=15*6000; endZ=20*6000; stdevThresh = 1;
%     start_time = 1; end_time = 30;
%     startZ = startFrameMs; endZ = endFrameMs; stdevThresh = 1;
%     startRO=50*6000; endRO=60*6000;
%     
%     
%     ZdfofCorrPart = (dfofCorr(startZ:endZ) - mean(dfofCorr(startZ:endZ)))  / std(dfofCorr(startZ:endZ));
%     y = mean (ZdfofCorrPart);
%     disp(y); 
%     
      start_time = str2double(answer(1)); end_time = str2double(answer(2));
    startZ=start_time*6000; endZ=end_time*6000; stdevThresh = 1;
    startRO=50*6000; endRO=60*6000;
    
    
    ZdfofCorrPart = (dfofCorr(startZ:endZ) - mean(dfofCorr(startZ:endZ)))  / std(dfofCorr(startZ:endZ));
    y = mean (ZdfofCorrPart);
    disp(y); 
%ZdfofCorrPart = (dfofCorr - mean(dfofCorr(startZ:endZ)))  /std(dfofCorr(startZ:endZ));
% the line above only changed the section for normalization for the full length of the data, but didn't cut the data   


% time_test = linspace(start_time,end_time,length(ZdfofCorrPart));
% figure;plot(time_test,ZdfofCorrPart,'LineWidth',1.5);
% set(gca,'FontSize',18);set(gcf,'Color','White'); 
% title('z-scored data');

    %ProcZpart1 = ZdfofCorrPart-stdevThresh;camcomment
    %numptsOvThreshBase = find(ProcZpart1(startfit:endfit) >0);camcomment
    %numptsOvThreshReadout = find(ProcZpart1(startfit:endfit)>0);camcomment
    
    %FractBase=length(numptsOvThreshBase)/(endfit-startfit);camcomment
    %FractRO = length(numptsOvThreshReadout)/(endfit-startfit);camcomment
    %numptsOvThreshBase = find(ProcZpart1(startZ:endZ) >0);
    %numptsOvThreshReadout = find(ProcZpart1(startRO:endRO)>0);
    
    %FractBase=length(numptsOvThreshBase)/(endZ-startZ);
    %FractRO = length(numptsOvThreshReadout)/(endRO-startRO);
    

savename2=strcat(savename,'_dFFcorr');
save(savename2,'dfofCorr','ZdfofCorr');%,'camtime','timeFP_RS');


time_test = linspace(start_time,end_time,length(ZdfofCorrPart));
y = mad(ZdfofCorrPart,1,'all');
thresholds = [y*1,y*1.25,y*1.5,y*1.75,y*2,y*2.25,y*2.5,y*2.75,y*3,y*3.25,y*3.5,y*3.75,y*4,y*4.25,y*4.5,y*4.75,y*5,y*5.25,y*5.5,y*5.75,y*6,y*6.25,y*6.5,y*6.75,y*7,y*7.25,y*7.5,y*7.75,y*8,y*8.25,y*8.5,y*8.75,y*9,y*9.25,y*9.5,y*9.75,y*10, y+1,y+1.25,y+1.5,y+1.75,y+2,y+2.25,y+2.5,y+2.75,y+3,y+3.25,y+3.5,y+3.75,y+4,y+4.25,y+4.5,y+4.75,y+5,y*2.91];

[pks,locs,width] = findpeaks(ZdfofCorrPart);
for i = 1:length(thresholds);
    [row,col] = find(pks>thresholds(i));
    target_pks = pks(col);
    target_locs = locs(col);
   
     if i==1;
        str_number = strcat('# of peaks above the threshold  is', {' '},num2str(length(target_pks)));%{' '} is just for space
        time_test_pk = time_test(target_locs);
        
        figure;plot(time_test,ZdfofCorrPart,'LineWidth',1.5); 
        hold on; yline(thresholds(3),'Color','r','LineWidth',2.5);
        hold on; plot(time_test_pk,target_pks,'*');text(mean(time_test)-0.5,max(pks)+1,str_number);
        title('z-scored data');
        set(gca,'FontSize',18);set(gcf,'Color','White');
        
        cd '/Users/kevin/Desktop/FP'; %folder pathway
        saveas(gcf,strcat(savename,'threshold'))
    end
    peak_numbers(i,:) = length(target_pks);
end

% for i = 1:length(thresholds);
%     for c = 2:length(ZdfofCorrPart)
%       [row_all,col_all] = find(ZdfofCorrPart(c) < thresholds(i) && ZdfofCorrPart(c-1) > thresholds(i));
%       tot_events(i,:) = length(col_all);
%     end
% end
for i = 1:length(thresholds);
    event_number(i,:) = 0;
   for c = 2:length(ZdfofCorrPart)
      if ZdfofCorrPart(c) > thresholds(i) && ZdfofCorrPart(c-1) <thresholds(i)
        event_number(i,:)= event_number(i,:) + 1;
      end
   end

end
% find points above the threshold
for i = 1:length(thresholds);
    [row_all,col_all] = find(ZdfofCorrPart>thresholds(i));
    ZdfofCorrPart_value = ZdfofCorrPart(col_all);
    tot_num_thresh(i,:) = length(col_all);
    mean_ZdfofCorrPart_value(i,:) = mean(ZdfofCorrPart_value);
end

for i = length(thresholds);
    total_AUC(i,:) = trapz(ZdfofCorrPart);
end

%output to an excel sheet
cd '/Users/kevin/Desktop/FP';
combined_info = [peak_numbers,tot_num_thresh,mean_ZdfofCorrPart_value];
T = table(tot_num_thresh,mean_ZdfofCorrPart_value,peak_numbers,event_number,'VariableNames',{'tot_num_thresh','mean_ZdfofCorrPart_value','peak_numbers','event_number'});
writetable(T,'combined_info_threshold.xlsx');


%% scroll to bottom if using video without Synapse timestamps

%% heatmap, esp for 3CT
% create variable trace(x,y), use nose (or center?)
% animal tracking can start at arbitrary pt, must continue to end of expt
% cam timestamps are 'missedTSnum' less than vid frames - vid keeps
% acquiring after Synapse stops timestamping
%lowRes = 1; %factor to reduce map density. native trace should be in mm
%startFr=1;%100;%20618;%23500;  %22865; 31865; 40864
% 'dur' = minutes after startFr to track  % b=frame for dur
%dur=15; [a b]=min(abs(camtime/60 - dur - camtime(startFr)/60)); b
%endFr=b;%length(camtime);%b;
%[mapDec,thismanyhitsDec,mapMaxDec] = heatmap3ct(trace,ZdfofCorr,camtime,timeFP_frames,startFr,endFr,lowRes,missedTSnum);

% ZdfofCorr=DFForZ; trace=traceOrig; missedTSnum=12;


%% dist from soc cup vs. max DFF
%clear DistFrMouse mapMaxVal;
%o=1;
%Mcup=[11,14];
 %   for m=1:size(mapMaxDec,1)
  %      for n=1:size(mapMaxDec,2)
%            DistFrMouse(o)=sqrt((Mcup(1) - m )^2 + (Mcup(2)-n)^2);
 %%           mapMaxVal(o)=mapMaxDec(m,n);
  %          o=o+1;
  %      end
  %  end
  %  
  %  DistFrMouseBin(1:ceil(max(DistFrMouse)))=0;
  %  for q=1:length(DistFrMouseBin); %each element i = i cms from mouse
  %      DistFrMouseBin(q) = max(  mapMaxVal(  find(~floor(abs(DistFrMouse - (q-1))))   ) );
  %      q=q+1;
  %  end
  %  DistFrMouseBin=DistFrMouseBin';
  %  figure;plot(DistFrMouse,mapMaxVal);hold on; plot(DistFrMouseBin,'black')
  %  save(strcat(savename,'DistVdff'),'DistFrMouseBin','Mcup');
%%
%figure;plot(timeFP_RS(1:1:length(timeFP_RSretime))/60,dfofCorr(1:1:length(timeFP_RS)),'black')    
%xlim([1 45]); ylim([-.12 .9])
%%
    % quantify total Max dF/F in an 18 pix square @ 1 pixel/cm
 %   clear
 %   figure; imagesc(mapMaxDec); colormap jet
 %   avgheat=mapDec./thismanyhitsDec;
 %   MaxLRavgLR(1:8) = 0;
    
    % define center of square for summing dF/F
 %   LcupR=14; LcupC=13; %left or up, Row and Column
 %   RcupR=60; RcupC=12; %right or down, Row and Column
    
 %   MaxLRavgLR(1)=sum(sum(mapMaxDec(LcupR-9:LcupR+9,LcupC-9:LcupC+9))); %left cup
 %   MaxLRavgLR(2)=sum(sum(mapMaxDec(RcupR-9:RcupR+9,RcupC-9:RcupC+9))); %right cup
 %   MaxLRavgLR(3)=sum(sum(avgheat(LcupR-9:LcupR+9,LcupC-9:LcupC+9))); %left cup
 %   MaxLRavgLR(4)=sum(sum(avgheat(RcupR-9:RcupR+9,RcupC-9:RcupC+9))); %right cup
 %   MaxLRavgLR(5)=sum(sum(mapDec(LcupR-9:LcupR+9,LcupC-9:LcupC+9))); %left cup
 %   MaxLRavgLR(6)=sum(sum(mapDec(RcupR-9:RcupR+9,RcupC-9:RcupC+9))); %right cup
 %   MaxLRavgLR(7)=sum(sum(thismanyhitsDec(LcupR-9:LcupR+9,LcupC-9:LcupC+9))); %left cup
 %   MaxLRavgLR(8)=sum(sum(thismanyhitsDec(RcupR-9:RcupR+9,RcupC-9:RcupC+9))); %right cup
    
%% new PSTH
%evTime is vector holding frames of interest
%Assume 100Hz resampling of FP data
%use dfofCorr or ZdfofCorr in FPpsth arguments
%timeFP_frames=timeFP_RS; %if using video not encoded by Synapse
%    Pre=300;
%    Post=400;
%    BaseStart=290;%Pre*.9;
%    BaseEnd=100;%Pre*.1;
%    Skip=400;
%    Norm =1; %=0 no normalization to baseline,   =1 normalize to baseline
%[PSTHarrNew,markFP,markFPfr]=FPpsth(camtime,timeFP_RS,timeFP_frames,ZdfofCorr,evTime,Pre,Post,BaseStart,BaseEnd,Skip,Norm); %camtime,FP time, proc FP, events, pre, post, start baseline @ -, end baseline @ -

%ylim([-0.14 .55]);

%plot FP vs time with events marked
%figure;plot(timeFP_RS,dfofCorr);hold on; plot(timeFP_RS,markFP);

%plot FP vs frame number with events marked
%figure;plot(timeFP_frames,dfofCorr);hold on; plot(timeFP_frames,markFPfr);


%% needs a [R,2] matrix of beginning and end frames for each bout
%[boutArray,boutInteg,CumBoutInteg] = boutDFF(bouts,dfofCorr);
%%


% plot mean of a subset of PSTHarrNew
%bb=1
%hold on; plot(mean(PSTHarrNew(:,bb:bb+9),2));



%% spectrogram
%figure;
%spectrogram(dfofCorr,kaiser(30000,10),24000,100000,100,'yaxis');ylim([0 1]);caxis([-40 0])
% 5 min window, overlaps by 4min (should be 1div/min), 100K pt freq
% resolution to get lower registers, 100Hz sampling


%% manual demodulate: check contoperator


        %% for video acquired in a separate program from FP
        % trim FP to match video length
 %       delayVid=10.2; %in seconds
  %      timeFP_RS=timeFP_RS-delayVid;
   %     camtime=camtime-delayVid;
   %     TrimFPtimeInd = timeFP_RS>=0; TrimCamtimeInd=camtime>=0;
   %     camtime = camtime(TrimCamtimeInd);
   %     timeFP_RS = timeFP_RS(TrimFPtimeInd); 
   %     dfofCorr=dfofCorr(TrimFPtimeInd);
   %     ZdfofCorr=ZdfofCorr(TrimFPtimeInd);
    
   %     savename3=strcat(savename,'_dFFcorrTrim');
   %     save(savename3,'dfofCorr','ZdfofCorr','camtime','timeFP_RS');
