function mrQ_run_onlyspgrT1(dir,outDir,inputData_spgr,inputData_seir,B1file, varArgIn)
% mrQ_run(dir,outDir,useSUNGRID,refFile,inputData_spgr,inputData_seir,B1file)

%% I. Create the initial structure

if notDefined('outDir')
    outDir = fullfile(dir,'mrQ');
end

% Creates the name of the output directory
if ~exist(outDir,'dir'); mkdir(outDir); end

% Creates the mrQ structure
mrQ = mrQ_Create(dir,[],outDir);

%     Create a file containing mrQpath, named after its 'ID' (taken from
%     its tempname). This allows for an easy use of SunGrid.
mrQ_createIDfile(mrQ);

% Set other parameters, such as SUNGRID and fieldstrength

if ~notDefined('varArgIn')
    if ~isempty(varArgIn)
        for ii = 1:2:numel(varArgIn)-1
            % Check to make sure that the argument is formatted properly
            mrQ = mrQ_Set(mrQ, varArgIn{ii}, varArgIn{ii+1});
        end
    end
end

%% II. Arrange the SPGR
% A specific arrange function for nimsfs, nifti, or using input for user

if ~isfield(mrQ,'ArrangeSPGR_Date');
    
    if ~notDefined('inputData_spgr')
        mrQ = mrQ_arrangeSPGR_nimsfs(mrQ,inputData_spgr);
    else
        mrQ = mrQ_arrangeSPGR_nimsfs(mrQ);
        
    end
else
    fprintf('Data was already arranged on %s \n',mrQ.ArrangeSPGR_Date)
end


%% IV. Initiate and align SPGR
%  parameters for aligning SPGR

%load(name);

if isfield(mrQ,'SPGR_init_done');
else
    mrQ.SPGR_init_done=0;
end

if     mrQ.SPGR_init_done==0
    
    % Keeps track of the variables we use.
    % For details, look inside the function.
    [mrQ.InitSPGR]=mrQ_initSPGR(mrQ.SPGR,mrQ.refIm,mrQ.mmPerVox,mrQ.interp,mrQ.skip,[],mrQ);
    mrQ.SPGR_init_done=1;
    
    save(mrQ.name,'mrQ');
    fprintf('\n  init SPGR - done!           \n');
else
    fprintf(' \n Loading init SPGR data            \n');
    
end

%%  V. Fit SPGR PD

if ~isfield(mrQ,'SPGR_LinearT1fit_done');
    
    mrQ.SPGR_LinearT1fit_done=0;
end

% clobber is implemented inside (we can add this to the inputs)
if (mrQ.SPGR_LinearT1fit_done==0);
    
    [mrQ.LinFit]=mrQfit_T1M0_Lin(mrQ,mrQ.InitSPGR.spgr_initDir);
    
    mrQ.SPGR_LinearT1fit_done=1;
    
    save(mrQ.name,'mrQ');
    
    fprintf('\n Fit linear T1 SPGR  - done!              \n');
else
    fprintf('\n Loading linearly fitted SPGR T1                \n');
    
end
