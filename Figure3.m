clearvars; clear all

dataFolder ='/Volumes/group_cogneuro/Matt/Data/Repeition suppresion IEM';
addpath([dataFolder '/Functions']);
partcount=1;
for participant =[1:8 10:11 13:17]
    %% loads the EEG data
    cd([dataFolder '/EEG'])
    
    load('-mat',strcat('/P-',num2str(participant),'(Combinedcleaned).set'))
    
    %% laods the behavioural data
    Conditions1=csvread(strcat(dataFolder,'/Session 1/SubS1-',sprintf('%05d',participant),'.csv'));
    if participant  ==1
        Conditions1=Conditions1(1:2700,:);
    end
    Conditions2=csvread(strcat(dataFolder,'/Session 2/SubS2-',sprintf('%05d',participant),'.csv'));
    Conditions =[Conditions1 ;Conditions2];
    
    for trial =1:5400
        temp =EEG.data(1:64,:,trial);
        TrialsOfInterest(trial,:)=max(max(abs(temp(:,:))'))<100;
        
    end
    
    ConditionsNew=Conditions(TrialsOfInterest,:);
    
    %%  baseline corrects the EEG data and changes its format
    baseind=EEG.times > 500 & EEG.times < 600;
    Data= bsxfun(@minus,EEG.data,mean(EEG.data(:,baseind,:),2));
    Data=permute(Data,[3 1 2]);
    Data=Data(TrialsOfInterest,1:64,:);
    Data= bsxfun(@minus,Data,mean(Data,2));
    
    %% makes the regression matrix
    clear chan_resp
    
    n_ori_chans = 9;
    
    make_basis_function = @(xx,mu) (cosd(xx-mu)).^(n_ori_chans-mod(n_ori_chans,2));
    
    
    xx = linspace(1,180,180);
    basis_set = nan(180,n_ori_chans);
    chan_center = linspace(180/n_ori_chans,180,n_ori_chans);
    
    for cc = 1:n_ori_chans
        basis_set(:,cc) = make_basis_function(xx,chan_center(cc));
    end
    % Build stimulus mask
    
    Angs =Conditions(TrialsOfInterest,13);
    Angs = mod(Angs,180);
    Angs(Angs==0) = 180; % so that we can use positive integer to index into stimulus mask vector
    
    stim_mask = zeros(length(Angs),length(xx));
    
    for tt = 1:size(stim_mask,1)  % loop over trials
        stim_mask(tt,Angs(tt))=1;
        
    end
    
    % Generate design matrix
    
    trnX = stim_mask*basis_set;
    
    %% makes sure there's an even number of trials per orientation
    
    trn_ou = unique(Angs);
    
    
    
    trn_repnum = nan(size(Data));
    n_trials_per_orientation = nan(length(trn_ou),1);
    for ii = 1:length(trn_ou)  % for each orientation in training set
        thisidx = Angs==trn_ou(ii);
        trn_repnum(thisidx) = 1:(sum(thisidx));
        n_trials_per_orientation(ii) = sum(thisidx);
        clear thisidx;
    end
    trn_repnum(trn_repnum>min(n_trials_per_orientation))=NaN;
    trng_cv = Angs(~isnan(trn_repnum));
    ConditionsNewNew =ConditionsNew(~isnan(trn_repnum),:);
    
    % cull all those trials from trn data (renamed trn_cv here for convenience)
    trn_cv_coeffs=Data(~isnan(trn_repnum),:,:);
    trnX_cv = trnX(~isnan(trn_repnum),:);
    trn_repnum = trn_repnum(~isnan(trn_repnum));
    
    %% set analysis time windows
    winsize=16;
    trainwins      = [EEG.times(1):1000/256:(EEG.times(end)-winsize)]';
    trainwins(:,2) = trainwins(:,1) + winsize;
    ntrains        = size(trainwins,1);
    %%
    chan_resp_cv_coeffs = nan(size(trn_cv_coeffs,1),length(chan_center),ntrains);
    
    n_reps = max(trn_repnum(:));
    for tt = 1:length(trainwins)
        itrain   = EEG.times >= trainwins(tt,1) & EEG.times <= trainwins(tt,2);
        
        TimeMean = mean(trn_cv_coeffs(:,:,itrain),3)';
        
        for ii = 1:n_reps
            trnidx = trn_repnum~=ii;
            tstidx = trn_repnum==ii;
            
            thistrn = TimeMean(:,trnidx);
            thistst = TimeMean(:,tstidx);
            C1 =trnX_cv(trnidx,:)';
            
            W   = squeeze(thistrn*C1'*pinv(C1*C1')); %OLS
            chan_resp_cv_coeffs(tstidx,:,tt)  = (pinv(W'*W)*W'*thistst)'; %multiply test data by weights from training data
            
            
        end
        sprintf('Time : %2.0f/%2.0f',tt,length(trainwins))
    end
    
    
    
    
    %% now let's plot coregistered reconstructions
    
    targ_ori = chan_center(round(length(chan_center)/2));
    targ_ori_idx = find(chan_center==targ_ori);
    
    chan_resp_cv_coeffs_shift = nan(size(chan_resp_cv_coeffs));
    for ii = 1:length(trn_ou)
        thisidx = trng_cv==trn_ou(ii);
        
        chan_resp_cv_coeffs_shift(thisidx,:,:) = circshift(chan_resp_cv_coeffs(thisidx,:,:), targ_ori_idx-find(trn_ou(ii)==chan_center) , 2 );
    end
    %%
    GeneralizationStore=nan(4,size(chan_resp_cv_coeffs_shift,2),size(chan_resp_cv_coeffs_shift,3));
    for cond =1:4
        if cond ==1
            trials=find(ConditionsNewNew(:,5)==1& ConditionsNewNew(:,3)==1& ConditionsNewNew(:,6)==1);
        elseif cond ==2
            trials=find(ConditionsNewNew(:,5)==2 & ConditionsNewNew(:,3)==1& ConditionsNewNew(:,6)==1);
        elseif cond ==3
            
            trials=find(ConditionsNewNew(:,5)==2 & ConditionsNewNew(:,3)==2& ConditionsNewNew(:,6)==1);
        elseif cond ==4
            
            trials=find(ConditionsNewNew(:,5)==1 & ConditionsNewNew(:,3)==2& ConditionsNewNew(:,6)==1);
        end
        GeneralizationStore(cond,:,:,:)=squeeze(mean(chan_resp_cv_coeffs_shift(trials,:,:,:),1));
    end
OverallStore(partcount,:,:,:,:)=GeneralizationStore;
partcount=partcount+1;
end


%% plot the results
chan_center= 20:20:180;
lb = [-10 60 5 -3];
ub = [10 120 90 3];
f=@(c,xdata)c(1).*exp(-0.5.*((xdata-c(2))./c(3)).^2) +c(4);
x0=[.5 100 30 0.2];
options = optimoptions('lsqcurvefit','Display','off');
Fits= zeros(15,4,EEG.pnts,4);
colt=1;


smoothtime=16;
ftype = 'gaussian'; % smooth the data with a gaussian filter (for regression analysis, not ERP/TF)

fsize = round(smoothtime*(256/1000)); % 32 ms = 8 samples at 250 Hz


for part = 1:15
    for cond =1:4
        
        temp = squeeze(OverallStore(part,cond,:,:));
        temp    = filtfast(temp,2,[],ftype,fsize);
        
        
        parfor time =1:size(temp,2)
            tmean = temp(:,time);
            
            [psy,resnorm,residual,exitflag,output] = lsqcurvefit(f,x0,chan_center,...
                tmean',lb,ub,options);
            Fits(part,cond,time,:)=psy;
            
        end
        sprintf('part: %2.0f | cond: %2.0f',part,cond)
    end
end

%%
clear h
close all
clusterNPerms=5000;
var=1;
a=figure;
a.Color = [1 1 1];
Cols=linspecer(2);
     winsize=16;

    trainwins      = [EEG.times(1):1000/256:(EEG.times(end)-winsize)]';
    
    
Times2plot=trainwins>500 & trainwins<1200;
Time=trainwins(Times2plot)-600;
for cond =1:2
    if cond ==1
        
        Cond1 = squeeze((Fits(:,1,Times2plot,var)))...
            +squeeze((Fits(:,3,Times2plot,var)));
        Cond2 = squeeze((Fits(:,2,Times2plot,var)))...
            +squeeze((Fits(:,4,Times2plot,var)));
        n=1:2;
    elseif cond ==2
        
        Cond1 = squeeze((Fits(:,1,Times2plot,var)))...
            +squeeze((Fits(:,2,Times2plot,var)));
        Cond2 = squeeze((Fits(:,3,Times2plot,var)))...
            +squeeze((Fits(:,4,Times2plot,var)));
        n=4:5;
    else
        
        Cond1 = squeeze((Fits(:,1,Times2plot,var)))...
            -squeeze((Fits(:,2,Times2plot,var)));
        Cond2 = squeeze((Fits(:,3,Times2plot,var)))...
            -squeeze((Fits(:,4,Times2plot,var)));
        
    end
    subplot(2,3,n);
    
    %%
    Cond1=(Cond1/2);
    Cond2=(Cond2/2);
    
    
    %%
    hold on
    
    [datobs, datrnd] = cluster_test_helper((Cond1-Cond2)', clusterNPerms);
    [h_mem, p_mem, ~] = cluster_test(datobs,datrnd,0,0.05,0.05);
    fprintf('Var : %2.0f | Cond %2.0f = %2.6f \n',var,cond,min(p_mem))
    pclustu = unique(p_mem);
    npclust = nnz(pclustu < 0.05);
    for ipclust = 1:npclust % extract time range of each significant cluster and show in figure
        currind  = p_mem == pclustu(ipclust);
        fill([min(Time(currind)),min(Time(currind)),max(Time(currind)),max(Time(currind))],[2,2,2,2],[0 0 0],'EdgeColor','none')
        h1=fill([min(Time(currind)),min(Time(currind)),max(Time(currind)),max(Time(currind))],[0,2,2,0],[0 0 0],'EdgeColor','none');
        set(h1,'facealpha',.2);
    end
    Time(currind)
    Cond1Mean=mean(Cond1,1);
    Cond2Mean=mean(Cond2,1);
    clear Cond1SE Cond2SE
    Cond1SE=std(Cond1)/sqrt(15);
    Cond2SE=std(Cond2)/sqrt(15);
    
    
    h{1}=shadedErrorBar(Time,Cond1Mean',Cond1SE',{'Color',Cols(1,:),'LineWidth',2},.5);
    h{2}=shadedErrorBar(Time,Cond2Mean,Cond2SE,{'Color',Cols(2,:),'LineWidth',2},.5);
    
    xlim([-100 600])
    if cond ==1
        c=legend([h{1}.mainLine,h{2}.mainLine],'Repeat','Alternate');
        c.FontSize=12;
        c.EdgeColor=[ 1 1 1];
        
        title('Repetition Suppression')
    else
        c=legend([h{1}.mainLine,h{2}.mainLine],'Expected','Unexpected');
        c.FontSize=12;
        c.EdgeColor=[ 1 1 1];
        
        title('Expectation')
        
    end
    
    axes1=gca;
    set(axes1,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
        'YAxisLocation','left','LineWidth',1.5,...
        'XMinorTick','on','YMinorTick','on'...
        ,'TickLength',[.015 .015],'Layer','bottom');
    if cond ==2
        xlabel('Time (ms)','FontWeight','bold');
    end
    ylim([-.25 1.3])
    ylabel('Channel amplitude (a.u.)','FontWeight','bold');
    
    box off
end

%%
TimeOI1= ([79 185])+600;
TimeOI1ind=dsearchn(trainwins,TimeOI1');

MeanTime1=mean(OverallStore(:,:,:,TimeOI1ind),4);


cols=linspecer(2);
LineStyle='-';
CorrectOri=[20:20:180]-90;
X1=-89:90;

lb = [-10 -90 5 -3];
ub = [10 90 90 3];
x0=[.5 0 30 0.2];


for cond =1:2
    
        tempT =MeanTime1;
 
    
    if cond ==1
        
        Cond1 = squeeze((tempT(:,1,:)))...
            +squeeze((tempT(:,3,:)));
        Cond2 = squeeze((tempT(:,2,:)))...
            +squeeze((tempT(:,4,:)));
        n=3;
    elseif cond ==2
        
        
        Cond1 = squeeze((tempT(:,1,:)))...
            +squeeze((tempT(:,2,:)));
        Cond2 = squeeze((tempT(:,3,:)))...
            +squeeze((tempT(:,4,:)));
        n=6;
    end
    subplot(2,3,n);
    
    Cond1=Cond1/2;
    Cond2=Cond2/2;
    Cond1Mean=mean(Cond1,1);
    Cond2Mean=mean(Cond2,1);
    Cond1SE=std(Cond1,1)/sqrt(15);
    Cond2SE=std(Cond2,1)/sqrt(15);
    hold on
    
    errorbar(CorrectOri,Cond1Mean,Cond1SE,'o',...
        'Color',cols(1,:),'MarkerFaceColor',cols(1,:),'CapSize',0)
    
    errorbar(CorrectOri,Cond2Mean,Cond2SE,'o','Color',...
        cols(2,:),'MarkerFaceColor',cols(2,:),'CapSize',0)
    
    [psy,resnorm,residual,exitflag,output] = lsqcurvefit(f,x0,CorrectOri,...
        Cond1Mean,lb,ub,options);
    
    plot(X1,f(psy,X1),'MarkerFaceColor',cols(1,:))
    
    
    [psy,resnorm,residual,exitflag,output] = lsqcurvefit(f,x0,CorrectOri,...
        Cond2Mean,lb,ub,options);
    
    plot(X1,f(psy,X1),'MarkerFaceColor',cols(2,:))
    
    ylim([0 .7])
    xticks(-90:90:90)
    if cond ==1
        c=legend('Repeat','Alternate');
        
        title('Repetition Suppression')
    else
        c=legend('Expected','Unexpected');
        
        title('Prediction Error')
    end
    
    axes1=gca;
    set(axes1,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
        'YAxisLocation','left','LineWidth',1.5,...
        'XMinorTick','on','YMinorTick','on'...
        ,'TickLength',[.015 .015],'Layer','bottom');
    if cond ==2
        xlabel('Channel orientation','FontWeight','bold');
    end
    ylabel('Channel amplitude (a.u.)','FontWeight','bold');
    
    box off
end
set(gcf,'Position',[   416   506   591   517])

%%

TimeOI1= [.79 .165]+.6;
TimeOI2= [.286 .3]+.6;
TimeOI1ind=dsearchn(EEG.times',TimeOI1');
TimeOI2ind=dsearchn(EEG.times',TimeOI2');



temp = squeeze(OverallStore(part,cond,:,:));

MeanTime1=mean(OverallStore(:,:,:,TimeOI1ind),4);
MeanTime2=mean(OverallStore(:,:,:,TimeOI2ind),4);

chan_center= 20:20:180;
lb = [-2 60 5 -1];
ub = [2 120 90 1];
X1=1:180;
f=@(c,xdata)c(1).*exp(-0.5.*((xdata-c(2))./c(3)).^2) +c(4);
x0=[.5 100 30 0.2];
options = optimoptions('lsqcurvefit','Display','off');
FitstoMean= zeros(15,2,2,2,4);
cols = linspecer(20);
for cond =1:2
    for time =1:2
        if time ==1
            tempT =MeanTime1;
        else
            tempT =MeanTime2;
        end
        
        subplot(2,1,cond);
        if cond ==1
            
            Cond1 = squeeze((tempT(:,1,:)))...
                +squeeze((tempT(:,3,:)));
            Cond2 = squeeze((tempT(:,2,:)))...
                +squeeze((tempT(:,4,:)));
        elseif cond ==2
            Cond1 = squeeze((tempT(:,1,:)))...
                +squeeze((tempT(:,2,:)));
            Cond2 = squeeze((tempT(:,3,:)))...
                +squeeze((tempT(:,4,:)));
        end
        
        for part = 1:15
            temp =Cond1(part,:);
            
            [psy,resnorm,residual,exitflag,output] = lsqcurvefit(f,x0,chan_center,...
                temp,lb,ub,options);
            
            FitstoMean(part,cond,1,time,:)=psy;
            temp =Cond2(part,:);
            
            [psy,resnorm,residual,exitflag,output] = lsqcurvefit(f,x0,chan_center,...
                temp,lb,ub,options);
            FitstoMean(part,cond,2,time,:)=psy;
            
            
            
        end
    end
    
end
%%
clc
var=1;
time =1;
cond=2;
[H,P,CI,STATS]=ttest(squeeze(FitstoMean(:,cond,1,time,var)),squeeze(FitstoMean(:,cond,2,time,var)))

squeeze(FitstoMean(:,cond,:,time,var))

%


