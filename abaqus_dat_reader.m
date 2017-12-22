function [Results, Step, Err] = abaqus_dat_reader(File)
    % Michael Jandron
    % Brown University
    % June 21, 2017
    %
    % ****************** ABAQUS .DAT FILE READER **********************
    %
    % Parses .dat file and saves .mat file of all results generated using
    % the following keywords:
    %  *NODE PRINT  (Works with *STATIC, *FREQUENCY, *STEADY STATE DYNAMIC, *COUPLED TEMPERATURE-DISPLACEMENT)
    %  *EL PRINT (Works with *STATIC, *FREQUENCY)
    %  *ENERGY PRINT (Works with *STATIC)
    % Additionally it parses *FREQUENCY data generated using the LANCZOS solver.
    %
    % Rev: 06/22/17 - Added ability to read coupled-temperature displacement steps
    %               - Now can read 8 variables per line, before it was 6
    % Rev: 07/01/17 - Added ability to read *ENERGY PRINT for static analyses
    % Rev: 07/09/17 - Added ability to read *EL PRINT for static analyses
    %               - Changed some of the variable names for node print and
    %                 energy to differentiate it from element print data
    % Rev: 08/11/17 - Added ability to read *NODE PRINT and *EL PRINT
    %                 inside *FREQUENCY steps (handy for visualizing mode shapes)
    %                 stored to Step.Mode structure
    % Rev: 12/22/17 - General update for git
	%

    if nargin==0
        [fname,path] = uigetfile('*.dat','Please select dat file');
        File = fullfile(path,fname);
    end
    
    fprintf('Reading file...')
    fid = fopen(File,'r');
    C = textscan(fid,'%s','Delimiter','\r\n','MultipleDelimsAsOne',1);
    fclose(fid);
    fprintf('complete\n')
    
    C = upper(C{1});
    N = length(C);

    % first find all steps
    fprintf('Identifying all steps...')
    S = 0;
    for k = 1:N
        isStep = contains(C{k},'S T E P');
        if isStep
            S = S + 1;
            Step(S).num = sscanf(C{k},'S T E P %d');
            Step(S).header = C{k};
            Step(S).position = k;
            Step(S).isStatic = contains(C{k},'S T A T I C');
            Step(S).isCoupledTempDisp = contains(C{k},'C O U P L E D - T E M P E R A T U R E - D I S P L A C E M E N T');
            Step(S).isEigen = contains(C{k},'E I G E N V A L U E S');
            Step(S).isSSD = contains(C{k},'S T E A D Y   S T A T E   A N A L Y S I S');
            %fprintf('%s\n',C{k})
        end
    end
    fprintf('complete\n')
    
    Results = []; Err = [];

    StaticSteps = find([Step(:).isStatic]);
    try
    if ~isempty(StaticSteps)
        fprintf('Processing static steps...')
        Step = ReadStaticSteps(C,Step,StaticSteps);
        fprintf('complete\n')
        fprintf('Condensing static steps...')
        Results = CondenseStaticSteps(Step,StaticSteps,Results);
        fprintf('complete\n')
    end
    catch Err
        fprintf('error in static routines, data may be unreliable\n')
    end
    
    CoupledTempDispSteps = find([Step(:).isCoupledTempDisp]);
    try
    if ~isempty(CoupledTempDispSteps)
        fprintf('Processing coupled temp-disp steps...')
        Step = ReadStaticSteps(C,Step,CoupledTempDispSteps); % same format as static
        fprintf('complete\n')
        fprintf('Condensing coupled temp-disp steps...')
        Results = CondenseStaticSteps(Step,CoupledTempDispSteps,Results); % same format as static
        fprintf('complete\n')
    end
    catch Err
        fprintf('error in static routines, data may be unreliable\n')
    end

    EigenSteps = find([Step(:).isEigen]);
    try
    if ~isempty(EigenSteps)
        fprintf('Processing eigenvalue steps...')
        [Step, EigenSteps] = ReadEigenSteps(C,Step,EigenSteps);
        fprintf('complete\n')
        fprintf('Consdensing eigenvalue steps...')
        Results = CondenseEigenSteps(Step,EigenSteps,Results);
        fprintf('complete\n')
    end
    catch Err
        fprintf('error in eigenvalue routines, data may be unreliable\n')
    end
    
    SSDSteps = find([Step(:).isSSD]);
    try
    if ~isempty(SSDSteps)
        fprintf('Processing SSD steps...')
        Step = ReadSSDSteps(C,Step,SSDSteps);
        fprintf('complete\n')
        fprintf('Condensing SSD steps...')
        Results = CondenseSSDSteps(Step,SSDSteps,Results);
        fprintf('complete\n')
    end
    catch Err
        fprintf('error in SSD routines, data may be unreliable\n')
    end
    
    fprintf('Saving matlab data file...')
    save(sprintf('%s.mat',File),'Step','Results');
    fprintf('complete\n')
    
return
function Step = ReadStaticSteps(C,Step,StaticSteps)
    try
        for S = StaticSteps

            if length(Step)>=S+1
                EndPos = Step(S+1).position;
            else
                EndPos = length(C);
            end
            % find all increments and table information
            I = 0; NPT = 0; EPT = 0;
            Step(S).hasEnergy = 0; Step(S).hasNodePrint = 0; Step(S).hasElPrint = 0;
            
            for LinePos = Step(S).position : EndPos
                isInc = contains(C{LinePos},'INCREMENT') && contains(C{LinePos},'SUMMARY');
                if isInc
                    I = I + 1; NPT = 0; EPT = 0;
                    Step(S).Inc(I).Num = sscanf(C{LinePos},'INCREMENT %d SUMMARY');
                    Step(S).Inc(I).header = C{LinePos};
                    Step(S).Inc(I).position = LinePos;
                    tmp1 = sscanf(C{LinePos+1},'TIME INCREMENT COMPLETED %f , FRACTION OF STEP COMPLETED %f');
                    Step(S).Inc(I).TimeIncCompleted = tmp1(1);
                    Step(S).Inc(I).FracStepCompleted = tmp1(2);
                    tmp1 = sscanf(C{LinePos+2},'STEP TIME COMPLETED %f , TOTAL TIME COMPLETED %f');
                    Step(S).Inc(I).StepTimeCompleted = tmp1(1);
                    Step(S).Inc(I).TotalTimeCompleted = tmp1(2);
                end
                isEnergy = contains(C{LinePos},'APPROXIMATE ENERGY TOTALS');
                if isEnergy
                    Step(S).hasEnergy = 1;
                    Step = ReadEnergy(Step,LinePos,C,S,I);
                end
                isNodePrintTable = contains(C{LinePos},'THE FOLLOWING TABLE IS PRINTED FOR NODES BELONGING TO NODE SET');
                if isNodePrintTable
                    Step(S).hasNodePrint = 1;
                    NPT = NPT + 1;
                    Step(S).Inc(I).NodePrintTable(NPT).SetName = sscanf(C{LinePos},'THE FOLLOWING TABLE IS PRINTED FOR NODES BELONGING TO NODE SET %s');
                    Step(S).Inc(I).NodePrintTable(NPT).header = C{LinePos};
                    Step(S).Inc(I).NodePrintTable(NPT).position = LinePos;
                    Variables = textscan(C{LinePos+1},'NODE FOOT- %s %s %s %s %s %s %s %s %s %s');
                    Step(S).Inc(I).NodePrintTable(NPT).Variables = [Variables{:}]; %flatten and remove empty cells
                    innerLinePos = LinePos + 3; indx = 1;

                    while ~(contains(C{innerLinePos},'MAXIMUM') || contains(C{innerLinePos},'INCREMENT') || contains(C{innerLinePos},'ALL VALUES IN THIS TABLE ARE ZERO'))
                        DataLine1 = sscanf(C{innerLinePos},'%i %f %f %f %f %f %f %f %f');
                        Step(S).Inc(I).NodePrintTable(NPT).Nodes(indx) = DataLine1(1);
                        Step(S).Inc(I).NodePrintTable(NPT).Data(:,indx) = DataLine1(2:end);
                        innerLinePos = innerLinePos + 1; indx = indx + 1;
                    end
                end
                isElPrintTable = contains(C{LinePos},'THE FOLLOWING TABLE IS PRINTED AT THE INTEGRATION POINTS FOR ELEMENT TYPE');
                if isElPrintTable
                    Step(S).hasElPrint = 1;
                    EPT = EPT + 1;
                    Step(S).Inc(I).ElPrintTable(EPT).Type = sscanf(C{LinePos},'THE FOLLOWING TABLE IS PRINTED AT THE INTEGRATION POINTS FOR ELEMENT TYPE %s');
                    Step(S).Inc(I).ElPrintTable(EPT).SetName = sscanf(C{LinePos+1},'%s');
                    Step(S).Inc(I).ElPrintTable(EPT).header = C{LinePos};
                    Step(S).Inc(I).ElPrintTable(EPT).position = LinePos;
                    Variables = textscan(C{LinePos+2},'ELEMENT  PT FOOT- %s %s %s %s %s %s %s %s %s %s');
                    Step(S).Inc(I).ElPrintTable(EPT).Variables = [Variables{:}]; %flatten and remove empty cells
                    innerLinePos = LinePos + 4; indx = 1;

                    while ~(contains(C{innerLinePos},'MAXIMUM') || contains(C{innerLinePos},'INCREMENT') ...
                            || contains(C{innerLinePos},'ALL VALUES IN THIS TABLE ARE ZERO') || contains(C{innerLinePos},'N O D E   O U T P U T') ...
                            || contains(C{innerLinePos},'THE FOLLOWING TABLE IS PRINTED AT THE INTEGRATION POINTS FOR ELEMENT TYPE'))
                        DataLine1 = sscanf(C{innerLinePos},'%i %f %f %f %f %f %f %f %f');
                        Step(S).Inc(I).ElPrintTable(EPT).Element(indx) = DataLine1(1);
                        Step(S).Inc(I).ElPrintTable(EPT).Point(indx) = DataLine1(2);
                        Step(S).Inc(I).ElPrintTable(EPT).Data(:,indx) = DataLine1(3:end);
                        innerLinePos = innerLinePos + 1; indx = indx + 1;
                    end
                end
            end
        end
	catch ME
        if (strcmp(ME.identifier,'MATLAB:badsubscript')) && (LinePos > length(C))
            fprintf('warning: end of file encountered...')
        else
            rethrow(ME);
        end
    end
return
function Step = ReadEnergy(Step,LinePos,C,S,I)
    Step(S).Inc(I).Energy_RSE = str2double(C{LinePos+1}(44:end)); % RECOVERABLE STRAIN ENERGY
    Step(S).Inc(I).Energy_CCEE = str2double(C{LinePos+2}(44:end)); % CONTACT CONSTRAINT ELASTIC ENERGY
    Step(S).Inc(I).Energy_CCENE = str2double(C{LinePos+3}(44:end)); % CONTACT CONSTRAINT ELASTIC NORMAL ENERG
    Step(S).Inc(I).Energy_CCETE = str2double(C{LinePos+4}(44:end)); % CONTACT CONSTRAINT ELASTIC TANG.  ENERGY 
    Step(S).Inc(I).Energy_KE = str2double(C{LinePos+5}(44:end)); % KINETIC ENERGY
    Step(S).Inc(I).Energy_EW = str2double(C{LinePos+6}(44:end)); % EXTERNAL WORK     
    Step(S).Inc(I).Energy_CCDW = str2double(C{LinePos+7}(44:end)); % CONTACT CONSTRAINT DISCONTINUITY WORK 
    Step(S).Inc(I).Energy_EWPCDW = str2double(C{LinePos+8}(44:end)); % EXTERNAL WORK plus CONTACT DISCONT. WORK 
    Step(S).Inc(I).Energy_PD = str2double(C{LinePos+9}(44:end)); % PLASTIC DISSIPATION  
    Step(S).Inc(I).Energy_CD = str2double(C{LinePos+10}(44:end)); % CREEP DISSIPATION 
    Step(S).Inc(I).Energy_VD = str2double(C{LinePos+11}(44:end)); % VISCOUS DISSIPATION (IN DAMPERS ETC) 
    Step(S).Inc(I).Energy_SD = str2double(C{LinePos+12}(44:end)); % STATIC DISSIPATION (STABILIZATION)     
    Step(S).Inc(I).Energy_SDICSD = str2double(C{LinePos+13}(44:end)); % INCLUDING CONTACT STAB. DISSIPATION 
    Step(S).Inc(I).Energy_CSND = str2double(C{LinePos+14}(44:end)); % CONTACT STABILIZATION NORMAL DISSIPATION
    Step(S).Inc(I).Energy_CSTF = str2double(C{LinePos+15}(44:end)); % CONTACT STABILIZATION TANG.  DISSIPATION
    Step(S).Inc(I).Energy_ELAI = str2double(C{LinePos+16}(44:end)); % ENERGY LOST AT IMPACTS   
    Step(S).Inc(I).Energy_ECSM = str2double(C{LinePos+17}(44:end)); % ENERGY TO CONTROL SPURIOUS MODES 
    Step(S).Inc(I).Energy_ELQB = str2double(C{LinePos+18}(44:end)); % ENERGY LOST THROUGH QUIET BOUNDARIES 
    Step(S).Inc(I).Energy_EE = str2double(C{LinePos+19}(44:end)); % ELECTROSTATIC ENERGY  
    Step(S).Inc(I).Energy_EDEC = str2double(C{LinePos+20}(44:end)); % ENERGY DUE TO ELECTRICAL CURRENT   
    Step(S).Inc(I).Energy_ELFD = str2double(C{LinePos+21}(44:end)); % ENERGY LOST TO FRICTIONAL DISSIPATION   
    Step(S).Inc(I).Energy_BD = str2double(C{LinePos+22}(44:end)); % BUCKLING DISSIPATION (FOR FRAME ELEMT.) 
    Step(S).Inc(I).Energy_DD = str2double(C{LinePos+23}(44:end)); % DAMAGE DISSIPATION  
    Step(S).Inc(I).Energy_TSESP = str2double(C{LinePos+24}(44:end)); % TOTAL STRAIN ENERGY (STRESS POWER) 
    Step(S).Inc(I).Energy_EB = str2double(C{LinePos+25}(44:end));% ENERGY BALANCE   
return
function Results = CondenseStaticSteps(Step,StaticSteps,Results)
    
    % condense data and prepare to initialize 3-D array for storage
    
    for S = StaticSteps
        
        if Step(S).hasNodePrint
            % list of times (dimension 1)
            time_array = [];
            for I = 1:length(Step(S).Inc)
                if isempty(time_array)
                    time_array = Step(S).Inc(1).TimeIncCompleted;
                else
                    time_array = [time_array, time_array(end) + Step(S).Inc(I).TimeIncCompleted];
                end
            end
            % unique list of nodes (dimension 2)
            common_nodes = [];
            for I = 1:length(Step(S).Inc)
                common_nodes = unique([common_nodes,[Step(S).Inc(I).NodePrintTable(:).Nodes]]);
            end
            % unique list of variables (dimension 3)
            common_vars = [];
            for I = 1:length(Step(S).Inc)
                common_vars = unique([common_vars,[Step(S).Inc(I).NodePrintTable(:).Variables]]);
            end

            Results(S).NodePrint.Time = time_array;
            Results(S).NodePrint.Nodes = common_nodes;
            Results(S).NodePrint.Vars = common_vars;
            Results(S).NodePrint.Ntime = length(time_array);
            Results(S).NodePrint.Nnodes = length(common_nodes);
            Results(S).NodePrint.Nvars = length(common_vars);

            % initialize data array
            Results(S).NodePrint.Data = zeros(Results(S).NodePrint.Ntime,Results(S).NodePrint.Nnodes,Results(S).NodePrint.Nvars);

            % now go fetch this data
            for I = 1:length(Step(S).Inc)
                for T = 1:length(Step(S).Inc(I).NodePrintTable)
                    [m,n] = size(Step(S).Inc(I).NodePrintTable(T).Data);
                    for V = 1:m % loop over variables
                        cvar = Step(S).Inc(I).NodePrintTable(T).Variables{V};
                        Vi = find(strcmp(cvar,Results(S).NodePrint.Vars));
                        for N = 1:n % loop over nodes
                            cnode = Step(S).Inc(I).NodePrintTable(T).Nodes(N);
                            Ni = find(cnode == Results(S).NodePrint.Nodes);
                            Results(S).NodePrint.Data(I,Ni,Vi) = Step(S).Inc(I).NodePrintTable(T).Data(V,N);
                        end
                    end
                end
            end
        end
        
        % Element Print
        if Step(S).hasElPrint
            % list of times (dimension 1)
            time_array = [];
            for I = 1:length(Step(S).Inc)
                if isempty(time_array)
                    time_array = Step(S).Inc(1).TimeIncCompleted;
                else
                    time_array = [time_array, time_array(end) + Step(S).Inc(I).TimeIncCompleted];
                end
            end
            % unique list of elements (dimension 2)
            common_elements = [];
            try
                for I = 1:length(Step(S).Inc)
                    common_elements = unique([common_elements,[Step(S).Inc(I).ElPrintTable(:).Element]]);
                end
            catch ME
                % this can only happen if ALL tables were zeros I think
                % abandon this
                continue
            end
            % unique list of variables (dimension 3)
            common_vars = [];
            for I = 1:length(Step(S).Inc)
                common_vars = unique([common_vars,[Step(S).Inc(I).ElPrintTable(:).Variables]]);
            end
            % unique list of integration points (dimension 4)
            common_points = [];
            for I = 1:length(Step(S).Inc)
                common_points = unique([common_points,[Step(S).Inc(I).ElPrintTable(:).Point]]);
            end

            Results(S).ElPrint.Time = time_array;
            Results(S).ElPrint.Elements = common_elements;
            Results(S).ElPrint.Vars = common_vars;
            Results(S).ElPrint.Points = common_points;
            Results(S).ElPrint.Ntime = length(time_array);
            Results(S).ElPrint.Nel = length(common_elements);
            Results(S).ElPrint.Nvars = length(common_vars);
            Results(S).ElPrint.Npoints = length(common_points);
            

            % initialize data array
            Results(S).ElPrint.Data = zeros(...
                Results(S).ElPrint.Ntime,Results(S).ElPrint.Nel,Results(S).ElPrint.Nvars,Results(S).ElPrint.Npoints);

            % now go fetch this data
            for I = 1:length(Step(S).Inc)
                for T = 1:length(Step(S).Inc(I).ElPrintTable)
                    [m,n] = size(Step(S).Inc(I).ElPrintTable(T).Data);
                    for V = 1:m % loop over variables
                        cvar = Step(S).Inc(I).ElPrintTable(T).Variables{V};
                        Vi = find(strcmp(cvar,Results(S).ElPrint.Vars));
                        for E = 1:n % loop over elements
                            cnode = Step(S).Inc(I).ElPrintTable(T).Element(E);
                            Ei = find(cnode == Results(S).ElPrint.Elements);
                            % use integration point as an index for speed
                            % instead of searching
                            Pi = Step(S).Inc(I).ElPrintTable(T).Point(E);
                            Results(S).ElPrint.Data(I,Ei,Vi,Pi) = Step(S).Inc(I).ElPrintTable(T).Data(V,E);
                        end
                    end
                end
            end

        end
            
        % if we have energy data, we can condense that as well (this is
        % done on a step by step basis)
        if Step(S).hasEnergy
            % condense data and prepare to initialize 3-D array for storage

            % list of times (dimension 1)
            time_array = [];
            for I = 1:length(Step(S).Inc)
                if isempty(time_array)
                    time_array = Step(S).Inc(1).TimeIncCompleted;
                else
                    time_array = [time_array, time_array(end) + Step(S).Inc(I).TimeIncCompleted];
                end
            end
            Results(S).Energy.Time = time_array;
            Results(S).Energy.Ntime = length(time_array);
            
            % unique list of energies (dimension 2)
            Results(S).Energy.header = {'RECOVERABLE STRAIN ENERGY', 'CONTACT CONSTRAINT ELASTIC ENERGY', ...
                'CONTACT CONSTRAINT ELASTIC NORMAL ENERGY', 'CONTACT CONSTRAINT ELASTIC TANG.  ENERGY', ...
                'KINETIC ENERGY', 'EXTERNAL WORK', 'CONTACT CONSTRAINT DISCONTINUITY WORK', ...
                'EXTERNAL WORK plus CONTACT DISCONT. WORK', 'PLASTIC DISSIPATION', 'CREEP DISSIPATION', ...
                'VISCOUS DISSIPATION (IN DAMPERS ETC)', 'STATIC DISSIPATION (STABILIZATION)', ...
                'INCLUDING CONTACT STAB. DISSIPATION', 'CONTACT STABILIZATION NORMAL DISSIPATION', ...
                'CONTACT STABILIZATION TANG.  DISSIPATION', 'ENERGY LOST AT IMPACTS', 'ENERGY TO CONTROL SPURIOUS MODES', ...
                'ENERGY LOST THROUGH QUIET BOUNDARIES', 'ELECTROSTATIC ENERGY', 'ENERGY DUE TO ELECTRICAL CURRENT', ...
                'ENERGY LOST TO FRICTIONAL DISSIPATION', 'BUCKLING DISSIPATION (FOR FRAME ELEMT.)', ...
                'DAMAGE DISSIPATION', 'TOTAL STRAIN ENERGY (STRESS POWER)', 'ENERGY BALANCE'};
            Results(S).Energy.Nval = 25;
            
            % initialize data array
            Results(S).Energy.Data = zeros(Results(S).Energy.Ntime,Results(S).Energy.Nval);

            % now go fetch this data
            for I = 1:length(Step(S).Inc)
                Results(S).Energy.Data(I,1) = Step(S).Inc(I).Energy_RSE;
                Results(S).Energy.Data(I,2) = Step(S).Inc(I).Energy_CCEE;
                Results(S).Energy.Data(I,3) = Step(S).Inc(I).Energy_CCENE;
                Results(S).Energy.Data(I,4) = Step(S).Inc(I).Energy_CCETE;
                Results(S).Energy.Data(I,5) = Step(S).Inc(I).Energy_KE;
                Results(S).Energy.Data(I,6) = Step(S).Inc(I).Energy_EW;
                Results(S).Energy.Data(I,7) = Step(S).Inc(I).Energy_CCDW;
                Results(S).Energy.Data(I,8) = Step(S).Inc(I).Energy_EWPCDW;
                Results(S).Energy.Data(I,9) = Step(S).Inc(I).Energy_PD;
                Results(S).Energy.Data(I,10) = Step(S).Inc(I).Energy_CD;
                Results(S).Energy.Data(I,11) = Step(S).Inc(I).Energy_VD;
                Results(S).Energy.Data(I,12) = Step(S).Inc(I).Energy_SD;
                Results(S).Energy.Data(I,13) = Step(S).Inc(I).Energy_SDICSD;
                Results(S).Energy.Data(I,14) = Step(S).Inc(I).Energy_CSND;
                Results(S).Energy.Data(I,15) = Step(S).Inc(I).Energy_CSTF;
                Results(S).Energy.Data(I,16) = Step(S).Inc(I).Energy_ELAI;
                Results(S).Energy.Data(I,17) = Step(S).Inc(I).Energy_ECSM;
                Results(S).Energy.Data(I,18) = Step(S).Inc(I).Energy_ELQB;
                Results(S).Energy.Data(I,19) = Step(S).Inc(I).Energy_EE;
                Results(S).Energy.Data(I,20) = Step(S).Inc(I).Energy_EDEC;
                Results(S).Energy.Data(I,21) = Step(S).Inc(I).Energy_ELFD;
                Results(S).Energy.Data(I,23) = Step(S).Inc(I).Energy_DD;
                Results(S).Energy.Data(I,24) = Step(S).Inc(I).Energy_TSESP;
                Results(S).Energy.Data(I,25) = Step(S).Inc(I).Energy_EB;
            end
        end
    end

return
function Step = ReadSSDSteps(C,Step,SSDSteps)
    try
        for S = SSDSteps

            if length(Step)>=S+1
                EndPos = Step(S+1).position;
            else
                EndPos = length(C);
            end
            % find all increments and table information
            I = 0; T = 0;
            for LinePos = Step(S).position : EndPos
                isInc = contains(C{LinePos},'INCREMENT NUMBER') && contains(C{LinePos},'AT FREQUENCY');
                if isInc
                    I = I + 1; T = 0;
                    tmp1 = sscanf(C{LinePos},'INCREMENT NUMBER %d AT FREQUENCY (CYCLES/TIME) = %f');
                    Step(S).Inc(I).Num = tmp1(1);
                    Step(S).Inc(I).Frequency = tmp1(2);
                    Step(S).Inc(I).header = C{LinePos};
                    Step(S).Inc(I).position = LinePos;
                end
                isTable = contains(C{LinePos},'THE FOLLOWING TABLE IS PRINTED FOR NODES BELONGING TO NODE SET');
                if isTable
                    T = T + 1;
                    Step(S).Inc(I).Table(T).SetName = sscanf(C{LinePos},'THE FOLLOWING TABLE IS PRINTED FOR NODES BELONGING TO NODE SET %s');
                    Step(S).Inc(I).Table(T).header = C{LinePos};
                    Step(S).Inc(I).Table(T).position = LinePos;
                    Variables = textscan(C{LinePos+1},'NODE FOOT- %s %s %s %s %s %s %s %s %s %s');
                    Step(S).Inc(I).Table(T).Variables = [Variables{:}]; %flatten and remove empty cells
                    innerLinePos = LinePos + 3; indx = 1;
                    expression = '\d-\d';
                    while ~(contains(C{innerLinePos},'MAXIMUM') || contains(C{innerLinePos},'INCREMENT'))
                        % Check that E has not vanished from string
                        % This occurs if it is to the power of 100 - so we will add it
                        % back if that is the case
                        tline = C{innerLinePos};
                        startIndex = regexp(tline,expression);
                        for kk = 1:length(startIndex)
                            tline = sprintf('%sE%s',tline(1:startIndex(kk)+(kk-1)),tline(startIndex(kk)+1+(kk-1):end));
                        end
                        DataLine1 = sscanf(tline,'%i %f %f %f %f %f %f %f %f');
                        tline = C{innerLinePos+1};
                        startIndex = regexp(tline,expression);
                        for kk = 1:length(startIndex)
                            tline = sprintf('%sE%s',tline(1:startIndex(kk)+(kk-1)),tline(startIndex(kk)+1+(kk-1):end));
                        end
                        DataLine2 = sscanf(tline,'%i SSD %f %f %f %f %f %f %f %f');
                        if DataLine1(1) ~= DataLine2(1); error('nodes don''t match?'); end
                        if length(DataLine1) ~= length(DataLine2); error('something went wrong reading these lines'); end
                        Step(S).Inc(I).Table(T).Nodes(indx) = DataLine1(1);
                        Step(S).Inc(I).Table(T).Data(:,indx) = complex(DataLine1(2:end), DataLine2(2:end));
                        
                        innerLinePos = innerLinePos + 2; indx = indx + 1;
                    end
                end
            end
        end
	catch ME
        if (strcmp(ME.identifier,'MATLAB:badsubscript')) && (LinePos > length(C))
            fprintf('warning: end of file encountered...')
        else
            rethrow(ME);
        end
    end
    
return
function Results = CondenseSSDSteps(Step,SSDSteps,Results)
    
    % condense data and prepare to initialize 3-D array for storage
    
    for S = SSDSteps
        % list of times (dimension 1)
        freq_array = zeros(length(Step(S).Inc),1);
        for I = 1:length(Step(S).Inc)
            freq_array(I) = Step(S).Inc(I).Frequency;
        end
        % unique list of nodes (dimension 2)
        common_nodes = [];
        for I = 1:length(Step(S).Inc)
            common_nodes = unique([common_nodes,[Step(S).Inc(I).Table(:).Nodes]]);
        end
        % unique list of variables (dimension 3)
        common_vars = [];
        for I = 1:length(Step(S).Inc)
            common_vars = unique([common_vars,[Step(S).Inc(I).Table(:).Variables]]);
        end
    
        Results(S).Freq = freq_array;
        Results(S).Nodes = common_nodes;
        Results(S).Vars = common_vars;
        Results(S).Nfreq = length(freq_array);
        Results(S).Nnodes = length(common_nodes);
        Results(S).Nvars = length(common_vars);

        % initialize data array
        Results(S).Data = zeros(Results(S).Nfreq,Results(S).Nnodes,Results(S).Nvars);

        % now go fetch this data
        for I = 1:length(Step(S).Inc)
            for T = 1:length(Step(S).Inc(I).Table)
                [m,n] = size(Step(S).Inc(I).Table(T).Data);
                for V = 1:m % loop over variables
                    cvar = Step(S).Inc(I).Table(T).Variables{V};
                    Vi = find(strcmp(cvar,Results(S).Vars));
                    for N = 1:n % loop over nodes
                        cnode = Step(S).Inc(I).Table(T).Nodes(N);
                        Ni = find(cnode == Results(S).Nodes);
                        Results(S).Data(I,Ni,Vi) = Step(S).Inc(I).Table(T).Data(V,N);
                    end
                end
            end
        end
    end
        
return
function [Step, EigenSteps] = ReadEigenSteps(C,Step,EigenSteps)
    try
        for S = EigenSteps
            %fprintf('Eigen: %d\n',Step(k).num)
            LinePos = Step(S).position;
            while ~contains(C{LinePos},'E I G E N V A L U E    O U T P U T')
                LinePos = LinePos + 1;
            end
            EigStart = LinePos+3;
            while ~(contains(C{LinePos},'P A R T I C I P A T I O N   F A C T O R S') || contains(C{LinePos},'NOTE'))
                LinePos = LinePos + 1;
            end
            EigEnd = LinePos-1;
            % if we find the line ***NOTE: EIGENVALUES MARKED WITH A * APPEAR TO BE RIGID BODY MODES
            % we need to increase start of Participation factors by one line
            if contains(C{LinePos},'NOTE')
                PartStart = LinePos+3;
            else
                PartStart = LinePos+2;
            end
            while ~contains(C{LinePos},'E F F E C T I V E   M A S S')
                LinePos = LinePos + 1;
            end
            PartEnd = LinePos-1;
            EffStart = LinePos+2;
            while ~contains(C{LinePos},'TOTAL')
                LinePos = LinePos + 1;
            end
            EffEnd = LinePos-1;

            % Grab eigenvalues
            numEig = EigEnd-EigStart;
            Step(S).Eig = zeros(numEig,6);
            indx = 0;
            for kk = EigStart:EigEnd
                indx = indx + 1;
                Step(S).Eig(indx,:) = sscanf(C{kk},'%i%*c %f %f %f %f %f');
            end

            % Grab participation
            numPart = PartEnd-PartStart;
            Step(S).Part = zeros(numPart,7);
            indx = 0;
            for kk = PartStart:PartEnd
                indx = indx + 1;
                Step(S).Part(indx,:) = sscanf(C{kk},'%i %f %f %f %f %f %f');
            end

            % Grab effective mass
            numEff = EffEnd-EffStart;
            Step(S).Eff = zeros(numEff,7);
            indx = 0;
            for kk = EffStart:EffEnd
                indx = indx + 1;
                Step(S).Eff(indx,:) = sscanf(C{kk},'%i %f %f %f %f %f %f');
            end

            % Do node print or element print data exist?
            % Search between start and end of step
            
            
            % Only one increment since this is a linear analysis step
            if length(Step)>=S+1
                EndPos = Step(S+1).position;
            else
                EndPos = length(C);
            end
            % find all table information
            NPT = 0; EPT = 0; I = 0;
            Step(S).hasNodePrint = 0; Step(S).hasElPrint = 0;
            for LinePos = Step(S).position : EndPos
                
                isInc = contains(C{LinePos},'E I G E N V A L U E    N U M B E R');
                if isInc
                    I = I + 1; NPT = 0; EPT = 0;
                    Step(S).Mode(I).Num = sscanf(C{LinePos},'E I G E N V A L U E    N U M B E R  %d');
                    Step(S).Mode(I).header = C{LinePos};
                    Step(S).Mode(I).position = LinePos;
                end
                
                isNodePrintTable = contains(C{LinePos},'THE FOLLOWING TABLE IS PRINTED FOR NODES BELONGING TO NODE SET');
                if isNodePrintTable
                    Step(S).hasNodePrint = 1;
                    NPT = NPT + 1;
                    Step(S).Mode(I).NodePrintTable(NPT).SetName = sscanf(C{LinePos},'THE FOLLOWING TABLE IS PRINTED FOR NODES BELONGING TO NODE SET %s');
                    Step(S).Mode(I).NodePrintTable(NPT).header = C{LinePos};
                    Step(S).Mode(I).NodePrintTable(NPT).position = LinePos;
                    Variables = textscan(C{LinePos+1},'NODE FOOT- %s %s %s %s %s %s %s %s %s %s');
                    Step(S).Mode(I).NodePrintTable(NPT).Variables = [Variables{:}]; %flatten and remove empty cells
                    innerLinePos = LinePos + 3; indx = 1;

                    while ~(contains(C{innerLinePos},'MAXIMUM') || contains(C{innerLinePos},'INCREMENT') || contains(C{innerLinePos},'ALL VALUES IN THIS TABLE ARE ZERO'))
                        DataLine1 = sscanf(C{innerLinePos},'%i %f %f %f %f %f %f %f %f');
                        Step(S).Mode(I).NodePrintTable(NPT).Nodes(indx) = DataLine1(1);
                        Step(S).Mode(I).NodePrintTable(NPT).Data(:,indx) = DataLine1(2:end);
                        innerLinePos = innerLinePos + 1; indx = indx + 1;
                    end
                end
                
                isElPrintTable = contains(C{LinePos},'THE FOLLOWING TABLE IS PRINTED AT THE INTEGRATION POINTS FOR ELEMENT TYPE');
                if isElPrintTable
                    Step(S).hasElPrint = 1;
                    EPT = EPT + 1;
                    Step(S).Mode(I).ElPrintTable(EPT).Type = sscanf(C{LinePos},'THE FOLLOWING TABLE IS PRINTED AT THE INTEGRATION POINTS FOR ELEMENT TYPE %s');
                    Step(S).Mode(I).ElPrintTable(EPT).SetName = sscanf(C{LinePos+1},'%s');
                    Step(S).Mode(I).ElPrintTable(EPT).header = C{LinePos};
                    Step(S).Mode(I).ElPrintTable(EPT).position = LinePos;
                    Variables = textscan(C{LinePos+2},'ELEMENT  PT FOOT- %s %s %s %s %s %s %s %s %s %s');
                    Step(S).Inc(I).ElPrintTable(EPT).Variables = [Variables{:}]; %flatten and remove empty cells
                    innerLinePos = LinePos + 4; indx = 1;

                    while ~(contains(C{innerLinePos},'MAXIMUM') || contains(C{innerLinePos},'INCREMENT') || ...
                            contains(C{innerLinePos},'ALL VALUES IN THIS TABLE ARE ZERO') || ...
                            contains(C{innerLinePos},'N O D E   O U T P U T') || ...
                            contains(C{innerLinePos},'THE FOLLOWING TABLE IS PRINTED AT THE INTEGRATION POINTS FOR ELEMENT TYPE '))
                        DataLine1 = sscanf(C{innerLinePos},'%i %f %f %f %f %f %f %f %f');
                        Step(S).Mode(I).ElPrintTable(EPT).Element(indx) = DataLine1(1);
                        Step(S).Mode(I).ElPrintTable(EPT).Point(indx) = DataLine1(2);
                        Step(S).Mode(I).ElPrintTable(EPT).Data(:,indx) = DataLine1(3:end);
                        innerLinePos = innerLinePos + 1; indx = indx + 1;
                    end
                end
            end
        end
    catch ME
        if (strcmp(ME.identifier,'MATLAB:badsubscript')) && (LinePos > length(C))
            fprintf('warning: end of file encountered...')
            % erase the last step and attempt to pass back the rest
            Step = Step(1:S-1);
            EigenSteps = EigenSteps(1:end-1);
        else
            rethrow(ME);
        end
    end

return
function Results = CondenseEigenSteps(Step,EigenSteps,Results)
    
    % condense data and prepare to initialize 3-D array for storage
    % do not condense node print or el print data right now
    
    % get envelope size first
    numfreq = length(EigenSteps);
    nummodes = zeros(numfreq,1);
    indx = 0;
    for S = EigenSteps
        indx = indx+1;
        nummodes(indx) = length(Step(S).Eig(:,3));
    end
    maxmodes = max(nummodes);
    Results(1).EigFreq = zeros(maxmodes,numfreq);
    Results(1).PartFactorX = zeros(maxmodes,numfreq);
    Results(1).PartFactorY = zeros(maxmodes,numfreq);
    Results(1).PartFactorZ = zeros(maxmodes,numfreq);
    Results(1).PartFactorXR = zeros(maxmodes,numfreq);
    Results(1).PartFactorYR = zeros(maxmodes,numfreq);
    Results(1).PartFactorZR = zeros(maxmodes,numfreq);
    Results(1).EffMassX = zeros(maxmodes,numfreq);
    Results(1).EffMassY = zeros(maxmodes,numfreq);
    Results(1).EffMassZ = zeros(maxmodes,numfreq);
    Results(1).EffMassXR = zeros(maxmodes,numfreq);
    Results(1).EffMassYR = zeros(maxmodes,numfreq);
    Results(1).EffMassZR = zeros(maxmodes,numfreq);
    indx = 0;
    for S = EigenSteps
        indx = indx+1;
        Results(1).EigFreq(1:nummodes(indx),indx) = Step(S).Eig(:,3);
        Results(1).PartFactorX(1:nummodes(indx),indx) = Step(S).Part(:,2);
        Results(1).PartFactorY(1:nummodes(indx),indx) = Step(S).Part(:,3);
        Results(1).PartFactorZ(1:nummodes(indx),indx) = Step(S).Part(:,4);
        Results(1).PartFactorXR(1:nummodes(indx),indx) = Step(S).Part(:,5);
        Results(1).PartFactorYR(1:nummodes(indx),indx) = Step(S).Part(:,6);
        Results(1).PartFactorZR(1:nummodes(indx),indx) = Step(S).Part(:,7);
        Results(1).EffMassX(1:nummodes(indx),indx) = Step(S).Eff(:,2);
        Results(1).EffMassY(1:nummodes(indx),indx) = Step(S).Eff(:,3);
        Results(1).EffMassZ(1:nummodes(indx),indx) = Step(S).Eff(:,4);
        Results(1).EffMassXR(1:nummodes(indx),indx) = Step(S).Eff(:,5);
        Results(1).EffMassYR(1:nummodes(indx),indx) = Step(S).Eff(:,6);
        Results(1).EffMassZR(1:nummodes(indx),indx) = Step(S).Eff(:,7);
    end

return
