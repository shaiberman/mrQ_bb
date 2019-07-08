function mrQ=mrQ_Call_AntsAlign_forSEIRMAP_SPGR_bb(mrQ)
% mrQ=mrQ_Call_AntsAlign_forSEIR_SPGR(mrQ)
%
% This function will try to register the SPGR T1 to a fitted T1 SEIR using
% ANTS. We will make sure that the registration is good by the gradient
% between the two images. in cases the registration is off up to a
% threshold (mrQ.QuantAntsThresh). Note that this threshold might be
% different between scanners and depends on image quality and epi
% artifacts.  We find 0.65 to be good on GE Discovery MR750 and 0.8 for
% Siemens Skyra. We will try different starting point. To do so we will
% register the SEIR image each time to a different inversion time image.
% Additionally, we will try to register the SPGR to an ACPC space or to
% native space. Last we will crop the SPGR image and remove dark area that
% can effect the registration.
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2016

if ~isfield(mrQ,'QuantAntsThresh')
    mrQ.QuantAntsThresh=0.65;
end
if ~isfield(mrQ,'ants_bm')
    mrQ.ants_bm=0;
end

if ~isfield(mrQ,'testR1_BM')
    mrQ.testR1_BM=0;
end

% In some cases we get better BrainMask for the SEIR map with R1.
%if this flag is on we will try registration with the BrainMask from R1.
testR1_BM=mrQ.testR1_BM;

tryR1mask=0; % for other iterations don't use BrainMask from R1.
use_TtSp=0;
if testR1_BM==1  % if testR1_BM is on there will be 3 trials
   error('the R1 BM mask is not relevant when a T1 map is given rather that SEIR ITs')
else
    iter_num=2;  % normally there will be be 2 trials
end

iter_start =1;

if isfield(mrQ,'Ants_Info') % If you rerun with flag testR1_BM and the previous iterations allready exists
    if isfield(mrQ.Ants_Info,'QuantAnts2High') && testR1_BM==1
        if any(mrQ.Ants_Info.QuantAnts2High(:)>0 & mrQ.Ants_Info.QuantAnts2High(:)<inf)
            iter_start=3;
            Ants_Info.QuantAnts2High=mrQ.Ants_Info.QuantAnts2High;
        end
    end
end
Ants_Info.QuantAnts2High(1,iter_start:iter_num)=inf;

for jj=iter_start:iter_num
    if jj==1
       
        spgr_initDir=mrQ.InitSPGR.spgr_initDir;
        SPGRFieldName='InitSPGR';
        useBM=mrQ.ants_bm;
        msg='[]';
        
    elseif jj==2
        %  In previous iteration we try to register with/without a bm - the opposite of what was done until now.
        
        spgr_initDir=mrQ.InitSPGR.spgr_initDir;
        SPGRFieldName='InitSPGR';
        useBM=~(mrQ.ants_bm);
        msg='change the BM input ';
        
    end
    
    if jj>1
        display(sprintf(['The SPGR - EPI registration is not good enough \nWe will ',msg,' and then try again.']));
    end
    
    Dir_SEIR = fullfile(mrQ.outDir, 'SEIR');
    if ~exist(Dir_SEIR,'dir'),mkdir(Dir_SEIR);end
    % make link to 1 map
    
    mapLink = fullfile(Dir_SEIR,'T1_SEIR_T1_BM.nii.gz');
    mrQ.SEIRfits{1}.Dir_SEIR_Align2Im=Dir_SEIR;
    mrQ.SEIRfits{1}.SEIR_fit_dir=Dir_SEIR;
    mrQ.SEIRfits{1}.SEIR_epi_T1bmfile = mapLink;
    mrQ.SEIRfits{1}.SEIR_epi_T1file = mapLink;
    mrQ.SEIRfits{1}.SEIR_epi_fitFile = mrQ.T1map_seir;
    
    mrQ.SEIRfitDone(1)=1;
    mrQ.Ants_Info=Ants_Info;
    
    bm = fullfile(Dir_SEIR,'T1_SEIR_BM.nii.gz');
    mrQ.SEIR_epi_Maskfile = bm;

    save(mrQ.name,'mrQ');
       
%%   Register high-resolution EPI image to low-resolution aligned T1 image

% create a temp folder for the ...

AntsPath = fullfile(Dir_SEIR,['ANTS_' num2str(jj)]);
if ~exist(AntsPath,'dir'), mkdir(AntsPath); end

if useBM==1    % Register based on brain only- skull stripped
    SPGR_T1_file=mrQ.LinFit.T1_LFit;
    SEIR_T1file=mrQ.SEIRfits{1}.SEIR_epi_T1bmfile;
elseif useBM==0   % Register without EPI mask - with skull
    SPGR_T1_file=mrQ.LinFit.T1_LFit_HM;
    SEIR_T1file = mrQ.SEIRfits{1}.SEIR_epi_T1file;
end    

if use_TtSp==1 && useBM==0
    SPGR_T1_file=mrQ.LinFit_NtSp.T1_LFit_HM;
elseif use_TtSp==1 && useBM==1
    SPGR_T1_file=mrQ.LinFit_NtSp.T1_LFit;
end

[WARP_SPGR_EPI,  T1_spgr_epi]= mrQ_NLANTS_SPGR2EPI(SEIR_T1file,SPGR_T1_file,[],AntsPath,{SPGR_T1_file});

MovingScaleConstat=1000;% to compare SPGR T1 in sec to SEIR T1 in ms.
% check the registration:
T1_spgr_epi = T1_spgr_epi{1};
[Ants_Info.QuantAntsScore]=mrQ_QuantAnts(SEIR_T1file,T1_spgr_epi,MovingScaleConstat);

% We will keep the files only if the are the current best Ants registration
if min(Ants_Info.QuantAnts2High(:)) > Ants_Info.QuantAntsScore
    
    
    Ants_Info.WARP_SPGR_EPI = WARP_SPGR_EPI;
    Ants_Info.T1_spgr_epi=T1_spgr_epi;
    
    % Saving the epi alignment parameters
    Ants_Info.SEIR_SPGR_Curent_AlignNums=[jj];
    Ants_Info.SEIR_SPGR_Curent_AlignDirs={Dir_SEIR, spgr_initDir};
    Ants_Info.SPGRFieldName=SPGRFieldName;
    
    %saving the SEIR files
    mrQ.SEIR_epi_T1file = mrQ.SEIRfits{1}.SEIR_epi_T1bmfile;
    %     mrQ.SEIR_epi_resnormfile = mrQ.SEIRfits{1}.SEIR_epi_resnormfile;
    %     mrQ.SEIR_epi_fitFile = mrQ.SEIRfits{1}.SEIR_epi_fitFile;
    %     mrQ.SEIR_epi_M0file = mrQ.SEIRfits{1}.SEIR_epi_M0file;
    %     mrQ.SEIR_epi_Maskfile =mrQ.SEIRfits{1}.SEIR_epi_Maskfile;
    %     mrQ.SEIR_epi_T1bmfile_fromR1=mrQ.SEIRfits{ii}.SEIR_epi_T1bmfile_fromR1;
    %     mrQ.SEIR_epi_T1bmfile=mrQ.SEIRfits{1}.SEIR_epi_T1bmfile;
    %     Ants_Info.Use_BrainMask = useBM;
end

%   check if the current (best) registration is good enough
if Ants_Info.QuantAntsScore< mrQ.QuantAntsThresh
    Ants_Info.QuantAnts2High(jj)=0;
    mrQ.Ants_Info=Ants_Info;
    save(mrQ.name,'mrQ');
    break
else
    % if the registration quality did not pass the threshold, we will
    % keep the values of the registration quality and try again
    
    Ants_Info.QuantAnts2High(jj)=Ants_Info.QuantAntsScore;
    mrQ.Ants_Info=Ants_Info;
    save(mrQ.name,'mrQ');
    
end
 
end