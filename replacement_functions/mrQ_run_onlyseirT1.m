function mrQ_run_onlyseirT1(dir,outDir,inputData_spgr,inputData_seir,B1file, varArgIn)
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


%% III. Perform SEIR fit
if ~notDefined('B1file')
    if exist(B1file,'file')
        mrQ.B1FileName=B1file;
        fprintf('Using the B1 map:  %s \n',B1file);
    else
        error('Can not find the B1 map:  %s \n',B1file);
    end
end

if  ~isfield(mrQ,'B1FileName')
    % Checks if B1 was defined by the user.
    % If not, we will use the SEIR data to map it.
    
    if isfield(mrQ,'SEIR_done');
    else
        mrQ.SEIR_done=0;
    end
    %%
    %
    %%
    if (mrQ.SEIR_done==0);
        
        % ARRANGE SEIR data
        if ~isfield(mrQ,'ArrangeSEIR_Date')
            
            if ~notDefined('inputData_seir')
                mrQ = mrQ_arrangeSEIR_nimsfs(mrQ,inputData_seir);
            else
                mrQ = mrQ_arrangeSEIR_nimsfs(mrQ);
            end
        else
            fprintf('Data was already arranged on %s \n',mrQ.ArrangeSEIR_Date)
        end
        %% align SEIR
        Dir_SEIR_Align = fullfile(mrQ.outDir, 'SEIR');
        SEIR_fit_dir = fullfile(Dir_SEIR_Align, 'fitT1_GS');
        if ~exist( SEIR_fit_dir,'dir'),     mkdir(SEIR_fit_dir);end
        
        [~, ~, ~, SEIRsaveData]=mrQ_initSEIR(mrQ,Dir_SEIR_Align,mrQ.alignFlag);    
        
        mrQ.SEIRfits = mrQ_fitSEIR_T1(Dir_SEIR_Align,[],0);
    end      
end