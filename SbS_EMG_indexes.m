function [Indexes, GaitPhases] = SbS_EMG_indexes(FILENAME, CHforIdxFilename, Nexus, toBePlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                                                                                       %%%%%%%%%%%%
%%%%%%%%%%%%  [Indexes, GaitPhases] = SbS_EMG_indexes(FILENAME, CHforIdxFilename, Nexus, toBePlot) %%%%%%%%%%%%
%%%%%%%%%%%%                                                                                       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                                                         %%%
%%%         Eurobench STEPbySTEP routine to calculate SbS EMG indexes                                       %%%
%%%                                                                                                         %%%
%%%         IN:     filename            --->    STRING                                                      %%%
%%%                                             name of the CSV file with EMG and                           %%%
%%%                                             acceleration data from the acquisition system               %%%
%%%                                                                                                         %%%
%%%                 CHforIdxFilename    --->    YAML file with couples of numbers                           %%%
%%%                                             CHforINDEXES to calculate cocontraction                     %%%
%%%                                             e.g. CHforINDEXES:                                          %%%
%%%                                                                 - [1 2]                                 %%%
%%%                                                                 - [1 3]                                 %%%
%%%                                                                 - [2 3]                                 %%%
%%%                                                                                                         %%%
%%%                 Nexus               --->    BOOL                                                        %%%
%%%                                             0 acquired with Delsys software                             %%%
%%%                                             1 acquired with Nexus software                              %%%
%%%                                                                                                         %%%
%%%                 toBePlot            --->    BOOL                                                        %%%
%%%                                             0 do not plot the data                                      %%%
%%%                                             1 plot the data                                             %%%
%%%                                                                                                         %%%
%%%         OUT:    Indexes             --->    STRUCT                                                      %%%
%%%                                             Fields:                                                     %%%
%%%                                                     ActivationIndexes                                   %%%
%%%                                                     CocontractionIndexes                                %%%
%%%                 GaitPhases          --->    STRUCT                                                      %%%
%%%                                             Fields:                                                     %%%
%%%                                                     Time                                                %%%
%%%                                                     Perc                                                %%%
%%%                                                                                                         %%%
%%%                                                                                                         %%%
%%%                                                                                                         %%%
%%%                 SUBROUTINES:    ReadYaml,LoadDelsysData, FiltButterLBH, EnvelopeHilbert, ...            %%%
%%%                                 discrete_integrate, cocontraction_winter, WriteYaml                     %%%
%%%                                                                                                         %%%
%%%                                                                                                         %%%
%%%         Author:     Marco Caimmi                                                                        %%%
%%%                     STIIMA Nationl Research Council of Italy                                            %%%
%%%                     marco.caimmi@stiima.cnr.it                                                          %%%
%%%                     marco.caimmi@gmail.com                                                              %%%
%%%                                                                                                         %%%
%%%         Year:       2021                                                                                %%%
%%%                                                                                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global GlobalVar
GlobalVar.ROOT = pwd;
GlobalVar.DATAdir = strcat(GlobalVar.ROOT,'\input');
GlobalVar.OUTPUTdir = strcat(GlobalVar.ROOT,'\output');
GlobalVar.TOOLSdir = strcat(GlobalVar.ROOT,'\tools');


% if nargin == 2
PHASES_filename = strcat(FILENAME(1:end-8),'phases','.yml');
% end

%%%%%%%% Load data and parameters
%
cd(GlobalVar.DATAdir)
    Phases = ReadYaml(PHASES_filename);
    CHforIdx = ReadYaml(CHforIdxFilename);
cd(GlobalVar.ROOT)
    Acq = LoadDelsysData(FILENAME, Nexus, toBePlot);
cd(GlobalVar.TOOLSdir)
EMG = Acq.EMG;



%
%%%%%%%%

%%%% EMG Elaboration
%
%% 50Hz filtering
EMG(1).Data = EMG(1).Data - mean(EMG(1).Data);
EMG(2) = EMG(1);  %Initialize EMG(2)
EMG(3) = EMG(1);  %Initialize EMG(3)
NumEMGCh = size( EMG(1).Data, 2 );
for iCh = 1: NumEMGCh
    EMG(2).Name{iCh} = strcat(EMG(1).Name{iCh}, 'Filtered');
    EMG(2).Name{iCh} = strcat(EMG(1).Name{iCh}, 'HilbEnvelop');
end

EMG(2).Data = FiltButterLBH(EMG(1).Data,5,[45 55],2000, 'Bs');

%%%% Envelope by Hilbert
%EMG(2).Data = EnvelopeHilbertButterworthFilt(DataFilt, EMG(1).sR);
EMG(3).Data = EnvelopeHilbert(EMG(2).Data, EMG(1).sR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ACTIVATION INDEX calculation
%
%%% Channels for activation index computation
%
%
Couples = [];
for iCouple = 1 : numel(CHforIdx.CHANNELS)
    Couples = [Couples str2num(cell2mat( CHforIdx.CHANNELS{iCouple}) )];
end
Channels = unique( Couples );
NumChElab = numel( Channels );

%%%% Structure GaitPhases
GaitPhases = struct();
%%%%
% GaitPhases.Time
% LEFT FOOT
GaitPhases.Time.LEFT.STANCE.type = 'vector';
GaitPhases.Time.LEFT.STANCE.value.data = [];
GaitPhases.Time.LEFT.STANCE.value.mean = [];
GaitPhases.Time.LEFT.STANCE.value.stdev = [];
%
GaitPhases.Time.LEFT.SWING.type = 'vector';
GaitPhases.Time.LEFT.SWING.value.data = [];
GaitPhases.Time.LEFT.SWING.value.mean = [];
GaitPhases.Time.LEFT.SWING.value.stdev = [];
%
GaitPhases.Time.LEFT.LOADRESPONSE.type = 'vector';
GaitPhases.Time.LEFT.LOADRESPONSE.value.data = [];
GaitPhases.Time.LEFT.LOADRESPONSE.value.mean = [];
GaitPhases.Time.LEFT.LOADRESPONSE.value.stdev = [];
%
GaitPhases.Time.LEFT.PRESWING.type = 'vector';
GaitPhases.Time.LEFT.PRESWING.value.data = [];
GaitPhases.Time.LEFT.PRESWING.value.mean = [];
GaitPhases.Time.LEFT.PRESWING.value.stdev = [];
%
GaitPhases.Time.LEFT.GAITCYCLE.type = 'vector';
GaitPhases.Time.LEFT.GAITCYCLE.value.data = [];
GaitPhases.Time.LEFT.GAITCYCLE.value.mean = [];
GaitPhases.Time.LEFT.GAITCYCLE.value.stdev = [];
%%%%
% RIGHT FOOT
GaitPhases.Time.RIGHT.STANCE.type = 'vector';
GaitPhases.Time.RIGHT.STANCE.value.data = [];
GaitPhases.Time.RIGHT.STANCE.value.mean = [];
GaitPhases.Time.RIGHT.STANCE.value.stdev = [];
%
GaitPhases.Time.RIGHT.SWING.type = 'vector';
GaitPhases.Time.RIGHT.SWING.value.data = [];
GaitPhases.Time.RIGHT.SWING.value.mean = [];
GaitPhases.Time.RIGHT.SWING.value.stdev = [];
%
GaitPhases.Time.RIGHT.LOADRESPONSE.type = 'vector';
GaitPhases.Time.RIGHT.LOADRESPONSE.value.data = [];
GaitPhases.Time.RIGHT.LOADRESPONSE.value.mean = [];
GaitPhases.Time.RIGHT.LOADRESPONSE.value.stdev = [];
%
GaitPhases.Time.RIGHT.PRESWING.type = 'vector';
GaitPhases.Time.RIGHT.PRESWING.value.data = [];
GaitPhases.Time.RIGHT.PRESWING.value.mean = [];
GaitPhases.Time.RIGHT.PRESWING.value.stdev = [];
%
GaitPhases.Time.RIGHT.GAITCYCLE.type = 'vector';
GaitPhases.Time.RIGHT.GAITCYCLE.value.data = [];
GaitPhases.Time.RIGHT.GAITCYCLE.value.mean = [];
GaitPhases.Time.RIGHT.GAITCYCLE.value.stdev = [];

% GaitPhases.Perc
% LEFT FOOT
GaitPhases.Perc.LEFT.STANCE.type = 'vector';
GaitPhases.Perc.LEFT.STANCE.value.data = [];
GaitPhases.Perc.LEFT.STANCE.value.mean = [];
GaitPhases.Perc.LEFT.STANCE.value.stdev = [];
%
GaitPhases.Perc.LEFT.SWING.type = 'vector';
GaitPhases.Perc.LEFT.SWING.value.data = [];
GaitPhases.Perc.LEFT.SWING.value.mean = [];
GaitPhases.Perc.LEFT.SWING.value.stdev = [];
%
GaitPhases.Perc.LEFT.LOADRESPONSE.type = 'vector';
GaitPhases.Perc.LEFT.LOADRESPONSE.value.data = [];
GaitPhases.Perc.LEFT.LOADRESPONSE.value.mean = [];
GaitPhases.Perc.LEFT.LOADRESPONSE.value.stdev = [];
%
GaitPhases.Perc.LEFT.PRESWING.type = 'vector';
GaitPhases.Perc.LEFT.PRESWING.value.data = [];
GaitPhases.Perc.LEFT.PRESWING.value.mean = [];
GaitPhases.Perc.LEFT.PRESWING.value.stdev = [];
%
GaitPhases.Perc.LEFT.GAITCYCLE.type = 'vector';
GaitPhases.Perc.LEFT.GAITCYCLE.value.data = [];
GaitPhases.Perc.LEFT.GAITCYCLE.value.mean = [];
GaitPhases.Perc.LEFT.GAITCYCLE.value.stdev = [];
%%%%
% RIGHT FOOT
GaitPhases.Perc.RIGHT.STANCE.type = 'vector';
GaitPhases.Perc.RIGHT.STANCE.value.data = [];
GaitPhases.Perc.RIGHT.STANCE.value.mean = [];
GaitPhases.Perc.RIGHT.STANCE.value.stdev = [];
%
GaitPhases.Perc.RIGHT.SWING.type = 'vector';
GaitPhases.Perc.RIGHT.SWING.value.data = [];
GaitPhases.Perc.RIGHT.SWING.value.mean = [];
GaitPhases.Perc.RIGHT.SWING.value.stdev = [];
%
GaitPhases.Perc.RIGHT.LOADRESPONSE.type = 'vector';
GaitPhases.Perc.RIGHT.LOADRESPONSE.value.data = [];
GaitPhases.Perc.RIGHT.LOADRESPONSE.value.mean = [];
GaitPhases.Perc.RIGHT.LOADRESPONSE.value.stdev = [];
%
GaitPhases.Perc.RIGHT.PRESWING.type = 'vector';
GaitPhases.Perc.RIGHT.PRESWING.value.data = [];
GaitPhases.Perc.RIGHT.PRESWING.value.mean = [];
GaitPhases.Perc.RIGHT.PRESWING.value.stdev = [];
%
GaitPhases.Perc.RIGHT.GAITCYCLE.type = 'vector';
GaitPhases.Perc.RIGHT.GAITCYCLE.value.data = [];
GaitPhases.Perc.RIGHT.GAITCYCLE.value.mean = [];
GaitPhases.Perc.RIGHT.GAITCYCLE.value.stdev = [];

%%%% Structure ActivationIndexes
ActivationIndexes = struct();
%%%%
% LEFT FOOT
ActivationIndexes.LEFT.STANCE.type = 'vector';
ActivationIndexes.LEFT.STANCE.value.data = [];
ActivationIndexes.LEFT.STANCE.value.mean = [];
ActivationIndexes.LEFT.STANCE.value.stdev = [];
%
ActivationIndexes.LEFT.SWING.type = 'vector';
ActivationIndexes.LEFT.SWING.value.data = [];
ActivationIndexes.LEFT.SWING.value.mean = [];
ActivationIndexes.LEFT.SWING.value.stdev = [];
%
ActivationIndexes.LEFT.LOADRESPONSE.type = 'vector';
ActivationIndexes.LEFT.LOADRESPONSE.value.data = [];
ActivationIndexes.LEFT.LOADRESPONSE.value.mean = [];
ActivationIndexes.LEFT.LOADRESPONSE.value.stdev = [];
%
ActivationIndexes.LEFT.PRESWING.type = 'vector';
ActivationIndexes.LEFT.PRESWING.value.data = [];
ActivationIndexes.LEFT.PRESWING.value.mean = [];
ActivationIndexes.LEFT.PRESWING.value.stdev = [];
%
ActivationIndexes.LEFT.GAITCYCLE.type = 'vector';
ActivationIndexes.LEFT.GAITCYCLE.value.data = [];
ActivationIndexes.LEFT.GAITCYCLE.value.mean = [];
ActivationIndexes.LEFT.GAITCYCLE.value.stdev = [];
%%%%
% RIGHT FOOT
ActivationIndexes.RIGHT.STANCE.type = 'vector';
ActivationIndexes.RIGHT.STANCE.value.data = [];
ActivationIndexes.RIGHT.STANCE.value.mean = [];
ActivationIndexes.RIGHT.STANCE.value.stdev = [];
%
ActivationIndexes.RIGHT.SWING.type = 'vector';
ActivationIndexes.RIGHT.SWING.value.data = [];
ActivationIndexes.RIGHT.SWING.value.mean = [];
ActivationIndexes.RIGHT.SWING.value.stdev = [];
%
ActivationIndexes.RIGHT.LOADRESPONSE.type = 'vector';
ActivationIndexes.RIGHT.LOADRESPONSE.value.data = [];
ActivationIndexes.RIGHT.LOADRESPONSE.value.mean = [];
ActivationIndexes.RIGHT.LOADRESPONSE.value.stdev = [];
%
ActivationIndexes.RIGHT.PRESWING.type = 'vector';
ActivationIndexes.RIGHT.PRESWING.value.data = [];
ActivationIndexes.RIGHT.PRESWING.value.mean = [];
ActivationIndexes.RIGHT.PRESWING.value.stdev = [];
%
ActivationIndexes.RIGHT.GAITCYCLE.type = 'vector';
ActivationIndexes.RIGHT.GAITCYCLE.value.data = [];
ActivationIndexes.RIGHT.GAITCYCLE.value.mean = [];
ActivationIndexes.RIGHT.GAITCYCLE.value.stdev = [];

%%%% Structure CocontractionIndexes
%
CocontractionIndexes = struct();
%%%%
% LEFT FOOT
CocontractionIndexes.LEFT.STANCE.type = 'vector';
CocontractionIndexes.LEFT.STANCE.value.data = [];
CocontractionIndexes.LEFT.STANCE.value.mean = [];
CocontractionIndexes.LEFT.STANCE.value.stdev = [];
%
CocontractionIndexes.LEFT.SWING.type = 'vector';
CocontractionIndexes.LEFT.SWING.value.data = [];
CocontractionIndexes.LEFT.SWING.value.mean = [];
CocontractionIndexes.LEFT.SWING.value.stdev = [];
%
CocontractionIndexes.LEFT.LOADRESPONSE.type = 'vector';
CocontractionIndexes.LEFT.LOADRESPONSE.value.data = [];
CocontractionIndexes.LEFT.LOADRESPONSE.value.mean = [];
CocontractionIndexes.LEFT.LOADRESPONSE.value.stdev = [];
%
CocontractionIndexes.LEFT.PRESWING.type = 'vector';
CocontractionIndexes.LEFT.PRESWING.value.data = [];
CocontractionIndexes.LEFT.PRESWING.value.mean = [];
CocontractionIndexes.LEFT.PRESWING.value.stdev = [];
%
CocontractionIndexes.LEFT.GAITCYCLE.type = 'vector';
CocontractionIndexes.LEFT.GAITCYCLE.value.data = [];
CocontractionIndexes.LEFT.GAITCYCLE.value.mean = [];
CocontractionIndexes.LEFT.GAITCYCLE.value.stdev = [];
%%%%
% RIGHT FOOT
CocontractionIndexes.RIGHT.STANCE.type = 'vector';
CocontractionIndexes.RIGHT.STANCE.value.data = [];
CocontractionIndexes.RIGHT.STANCE.value.mean = [];
CocontractionIndexes.RIGHT.STANCE.value.stdev = [];
%
CocontractionIndexes.RIGHT.SWING.type = 'vector';
CocontractionIndexes.RIGHT.SWING.value.data = [];
CocontractionIndexes.RIGHT.SWING.value.mean = [];
CocontractionIndexes.RIGHT.SWING.value.stdev = [];
%
CocontractionIndexes.RIGHT.LOADRESPONSE.type = 'vector';
CocontractionIndexes.RIGHT.LOADRESPONSE.value.data = [];
CocontractionIndexes.RIGHT.LOADRESPONSE.value.mean = [];
CocontractionIndexes.RIGHT.LOADRESPONSE.value.stdev = [];
%
CocontractionIndexes.RIGHT.PRESWING.type = 'vector';
CocontractionIndexes.RIGHT.PRESWING.value.data = [];
CocontractionIndexes.RIGHT.PRESWING.value.mean = [];
CocontractionIndexes.RIGHT.PRESWING.value.stdev = [];
%
CocontractionIndexes.RIGHT.GAITCYCLE.type = 'vector';
CocontractionIndexes.RIGHT.GAITCYCLE.value.data = [];
CocontractionIndexes.RIGHT.GAITCYCLE.value.mean = [];
CocontractionIndexes.RIGHT.GAITCYCLE.value.stdev = [];


Time = EMG(3).Time;
Data = EMG(3).Data';

%%%%%
%%% Define Begin and End samples for each gait phaase
%%
%  L_ or R_ BST stands for ... Left/Right BeginStance
%  L_ or R_ EST stands for ... Left/Right EndStance
%  L_ or R_ BSW stands for ... Left/Right BeginSwing
%  L_ or R_ ESW stands for ... Left/Right EndSwing

NumLeftSteps = numel( str2num(Phases.LeftBeginStance) );
L_BST = str2num( Phases.LeftBeginStance );
L_EST = str2num( Phases.LeftEndStance ); % same as ... left BeginSwing
L_BSW = L_EST;
L_ESW = str2num( Phases.LeftEndSwing );

NumRightSteps = numel( str2num(Phases.RightBeginStance) );
R_BST = str2num( Phases.RightBeginStance );
R_EST = str2num( Phases.RightEndStance ); % same as ... right BeginSwing
R_BSW = R_EST;
R_ESW = str2num( Phases.RightEndSwing );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%         GAIT PHASES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LEFT FOOT
% stance
GaitPhases.Time.LEFT.STANCE.value.data = [L_EST - L_BST];
GaitPhases.Time.LEFT.STANCE.value.mean = mean(GaitPhases.Time.LEFT.STANCE.value.data);
GaitPhases.Time.LEFT.STANCE.value.stdev = std(GaitPhases.Time.LEFT.STANCE.value.data);
% swing
GaitPhases.Time.LEFT.SWING.value.data = [L_ESW - L_BSW];
GaitPhases.Time.LEFT.SWING.value.mean = mean(GaitPhases.Time.LEFT.SWING.value.data);
GaitPhases.Time.LEFT.SWING.value.stdev = std(GaitPhases.Time.LEFT.SWING.value.data);
% load response
if L_BST(1) > R_BST(1)
    BT = L_BST; % BS stands for Begin Time
    ET = R_EST; % ES stands for End Time
else
    BT = L_BST(2:end); % BS stands for Begin Time
    ET = R_EST(1:end-1); % ES stands for End Time
end
GaitPhases.Time.LEFT.LOADRESPONSE.value.data = [ET - BT];
GaitPhases.Time.LEFT.LOADRESPONSE.value.mean = mean(GaitPhases.Time.LEFT.LOADRESPONSE.value.data);
GaitPhases.Time.LEFT.LOADRESPONSE.value.stdev = std(GaitPhases.Time.LEFT.LOADRESPONSE.value.data);
% Preswing
if L_BST(1) > R_BST(1)
    BT = R_ESW; % BS stands for Begin Time
    ET = L_EST; % ES stands for End Time
else
    BT = R_BST; % BS stands for Begin Time
    ET = L_EST; % ES stands for End Time
end
GaitPhases.Time.LEFT.PRESWING.value.data = [ET - BT];
GaitPhases.Time.LEFT.PRESWING.value.mean = mean(GaitPhases.Time.LEFT.PRESWING.value.data);
GaitPhases.Time.LEFT.PRESWING.value.stdev = std(GaitPhases.Time.LEFT.PRESWING.value.data);
% gait cycle
GaitPhases.Time.LEFT.GAITCYCLE.value.data = [L_ESW - L_BST];
GaitPhases.Time.LEFT.GAITCYCLE.value.mean = mean(GaitPhases.Time.LEFT.GAITCYCLE.value.data);
GaitPhases.Time.LEFT.GAITCYCLE.value.stdev = std(GaitPhases.Time.LEFT.GAITCYCLE.value.data);
%%%%
% RIGHT FOOT
% Stance
GaitPhases.Time.RIGHT.STANCE.value.data = [R_EST - R_BST];
GaitPhases.Time.RIGHT.STANCE.value.mean = mean(GaitPhases.Time.RIGHT.STANCE.value.data);
GaitPhases.Time.RIGHT.STANCE.value.stdev = std(GaitPhases.Time.RIGHT.STANCE.value.data);
% swing
GaitPhases.Time.RIGHT.SWING.value.data = [R_ESW - R_BSW];
GaitPhases.Time.RIGHT.SWING.value.mean = mean(GaitPhases.Time.RIGHT.SWING.value.data);
GaitPhases.Time.RIGHT.SWING.value.stdev = std(GaitPhases.Time.RIGHT.SWING.value.data);
% load response
if R_BST(1) > L_BST(1)
    BT = R_BST; % BS stands for Begin Time
    ET = L_EST; % ES stands for End Time
else
    BT = R_BST(2:end); % BS stands for Begin Time
    ET = L_EST(1:end-1); % ES stands for End Time
end
GaitPhases.Time.RIGHT.LOADRESPONSE.value.data = [ET - BT];
GaitPhases.Time.RIGHT.LOADRESPONSE.value.mean = mean(GaitPhases.Time.RIGHT.LOADRESPONSE.value.data);
GaitPhases.Time.RIGHT.LOADRESPONSE.value.stdev = std(GaitPhases.Time.RIGHT.LOADRESPONSE.value.data);
% Preswing
if R_BST(1) > L_BST(1)
    BT = L_ESW; % BS stands for Begin Time
    ET = R_EST; % ES stands for End Time
else
    BT = L_BST; % BS stands for Begin Time
    ET = R_EST; % ES stands for End Time
end
GaitPhases.Time.RIGHT.PRESWING.value.data = [ET - BT];
GaitPhases.Time.RIGHT.PRESWING.value.mean = mean(GaitPhases.Time.RIGHT.PRESWING.value.data);
GaitPhases.Time.RIGHT.PRESWING.value.stdev = std(GaitPhases.Time.RIGHT.PRESWING.value.data);
% gait cycle
GaitPhases.Time.RIGHT.GAITCYCLE.value.data = [R_ESW - R_BST];
GaitPhases.Time.RIGHT.GAITCYCLE.value.mean = mean(GaitPhases.Time.RIGHT.GAITCYCLE.value.data);
GaitPhases.Time.RIGHT.GAITCYCLE.value.stdev = std(GaitPhases.Time.RIGHT.GAITCYCLE.value.data);
%%%%
cd(GlobalVar.OUTPUTdir)
    Outputfilename = strcat( FILENAME(1:end-8), '_GaitPhases');
    WriteYaml(Outputfilename, GaitPhases.Time,0);
cd(GlobalVar.TOOLSdir)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%         ACTIVATION INDEXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iCh = 1 : NumChElab
    for iStp = 1 : NumLeftSteps
        % LeftStance
        BS = floor( L_BST( iStp ) * EMG(3).sR ); % BS stands for Begin Sample
        ES = floor( L_EST( iStp ) * EMG(3).sR ) ; % ES stands for End Sample
        ActivationIndexes(iCh).LEFT.STANCE.value.data( iStp )  = discrete_integrate( Time(BS:ES), Data(iCh, BS:ES) );
        
        % LeftSwing
        BS = floor( L_BSW( iStp ) * EMG(3).sR ); % BS stands for Begin Sample
        ES = floor( L_ESW( iStp ) * EMG(3).sR ) ; % ES stands for End Sample
        ActivationIndexes(iCh).LEFT.SWING.value.data( iStp )  = discrete_integrate( Time(BS:ES), Data(iCh, BS:ES) );
        
        % LeftLoadResponse
        if L_BST(1) > R_BST(1)
            BS = floor( L_BST( iStp ) * EMG(3).sR ); % BS stands for Begin Sample
            ES = floor( R_EST( iStp ) * EMG(3).sR ) ; % ES stands for End Sample
        ActivationIndexes(iCh).LEFT.LOADRESPONSE.value.data( iStp )  = discrete_integrate( Time(BS:ES), Data(iCh, BS:ES) );
        else
            if iStp ~= NumLeftSteps
                BS = floor( L_BST( iStp + 1) * EMG(3).sR ); % BS stands for Begin Sample
                ES = floor( R_EST( iStp ) * EMG(3).sR ) ; % ES stands for End Sample
                ActivationIndexes(iCh).LEFT.LOADRESPONSE.value.data( iStp )  = discrete_integrate( Time(BS:ES), Data(iCh, BS:ES) );            
            end
        end
        
        % LeftPreSwing
        if L_BST(1) > R_BST(1)
            BS = floor( R_ESW( iStp ) * EMG(3).sR ); % BS stands for Begin Sample
            ES = floor( L_EST( iStp ) * EMG(3).sR ) ; % ES stands for End Sample
        else
            BS = floor( R_BST( iStp ) * EMG(3).sR ); % BS stands for Begin Sample
            ES = floor( L_EST( iStp ) * EMG(3).sR ) ; % ES stands for End Sample
        end
        ActivationIndexes(iCh).LEFT.PRESWING.value.data( iStp )  = discrete_integrate( Time(BS:ES), Data(iCh, BS:ES) );           
        
        % leftGaitCycle
        BS = floor( L_BST( iStp ) * EMG(3).sR ); % BS stands for Begin Sample
        ES = floor( L_ESW( iStp ) * EMG(3).sR ) ; % ES stands for End Sample
        ActivationIndexes(iCh).LEFT.GAITCYCLE.value.data( iStp )  = discrete_integrate( Time(BS:ES), Data(iCh, BS:ES) );
    end
    % left mean and standard deviations
    ActivationIndexes(iCh).LEFT.STANCE.value.mean = mean( ActivationIndexes(iCh).LEFT.STANCE.value.data );
    ActivationIndexes(iCh).LEFT.STANCE.value.stdev = std( ActivationIndexes(iCh).LEFT.STANCE.value.data );
    ActivationIndexes(iCh).LEFT.SWING.value.mean = mean( ActivationIndexes(iCh).LEFT.SWING.value.data );
    ActivationIndexes(iCh).LEFT.SWING.value.stdev = std( ActivationIndexes(iCh).LEFT.SWING.value.data );
    ActivationIndexes(iCh).LEFT.LOADRESPONSE.value.mean = mean( ActivationIndexes(iCh).LEFT.LOADRESPONSE.value.data );
    ActivationIndexes(iCh).LEFT.LOADRESPONSE.value.stdev = std( ActivationIndexes(iCh).LEFT.LOADRESPONSE.value.data );
    ActivationIndexes(iCh).LEFT.PRESWING.value.mean = mean( ActivationIndexes(iCh).LEFT.PRESWING.value.data );
    ActivationIndexes(iCh).LEFT.PRESWING.value.stdev = std( ActivationIndexes(iCh).LEFT.PRESWING.value.data );
    ActivationIndexes(iCh).LEFT.GAITCYCLE.value.mean = mean( ActivationIndexes(iCh).LEFT.GAITCYCLE.value.data );
    ActivationIndexes(iCh).LEFT.GAITCYCLE.value.stdev = std( ActivationIndexes(iCh).LEFT.GAITCYCLE.value.data );
    
    %%%
       
    for iStp = 1 : NumRightSteps
        % RightStance
        BS = floor( R_BST( iStp ) * EMG(3).sR ); % BS stands for Begin Sample
        ES = floor( R_EST( iStp ) * EMG(3).sR ) ; % ES stands for End Sample
        ActivationIndexes(iCh).RIGHT.STANCE.value.data( iStp )  = discrete_integrate( Time(BS:ES), Data(iCh, BS:ES) );
        
        % RightSwing
        BS = floor( R_BSW( iStp ) * EMG(3).sR ); % BS stands for Begin sample
        ES = floor( R_ESW( iStp ) * EMG(3).sR ) ; % ES stands for End sample
        ActivationIndexes(iCh).RIGHT.SWING.value.data( iStp )  = discrete_integrate( Time(BS:ES), Data(iCh, BS:ES) );
        
        % RightLoadResponse
        if R_BST(1) > L_BST(1)
            BS = floor( R_BST( iStp ) * EMG(3).sR ); % BS stands for Begin sample
            ES = floor( L_EST( iStp ) * EMG(3).sR ) ; % ES stands for End sample
        ActivationIndexes(iCh).RIGHT.LOADRESPONSE.value.data( iStp )  = discrete_integrate( Time(BS:ES), Data(iCh, BS:ES) );    
        else
            if iStp ~= NumRightSteps
                BS = floor( R_BST( iStp + 1) * EMG(3).sR ); % BS stands for Begin sample
                ES = floor( L_EST( iStp ) * EMG(3).sR ) ; % ES stands for End sample
                ActivationIndexes(iCh).RIGHT.LOADRESPONSE.value.data( iStp )  = discrete_integrate( Time(BS:ES), Data(iCh, BS:ES) );
            end
        end
        
        % RightPreSwing
        if R_BST(1) > L_BST(1)
            BS = floor( L_ESW( iStp ) * EMG(3).sR ); % BS stands for Begin sample
            ES = floor( R_EST( iStp ) * EMG(3).sR ) ; % ES stands for End sample
        else
            BS = floor( L_BST( iStp ) * EMG(3).sR ); % BS stands for Begin sample
            ES = floor( R_EST( iStp ) * EMG(3).sR ) ; % ES stands for End sample
        end
        ActivationIndexes(iCh).RIGHT.PRESWING.value.data( iStp )  = discrete_integrate( Time(BS:ES), Data(iCh, BS:ES) );

        % leftGaitCycle
        BS = floor( R_BST( iStp ) * EMG(3).sR ); % BS stands for Begin sample
        ES = floor( R_ESW( iStp ) * EMG(3).sR ) ; % ES stands for End sample
        ActivationIndexes(iCh).RIGHT.GAITCYCLE.value.data( iStp )  = discrete_integrate( Time(BS:ES), Data(iCh, BS:ES) );
    end
    % right mean and standard deviations
    ActivationIndexes(iCh).RIGHT.STANCE.value.mean = mean( ActivationIndexes(iCh).RIGHT.STANCE.value.data );
    ActivationIndexes(iCh).RIGHT.STANCE.value.stdev = std( ActivationIndexes(iCh).RIGHT.STANCE.value.data );
    ActivationIndexes(iCh).RIGHT.SWING.value.mean = mean( ActivationIndexes(iCh).RIGHT.SWING.value.data );
    ActivationIndexes(iCh).RIGHT.SWING.value.stdev = std( ActivationIndexes(iCh).RIGHT.SWING.value.data );
    ActivationIndexes(iCh).RIGHT.LOADRESPONSE.value.mean = mean( ActivationIndexes(iCh).RIGHT.LOADRESPONSE.value.data );
    ActivationIndexes(iCh).RIGHT.LOADRESPONSE.value.stdev = std( ActivationIndexes(iCh).RIGHT.LOADRESPONSE.value.data );
    ActivationIndexes(iCh).RIGHT.PRESWING.value.mean = mean( ActivationIndexes(iCh).RIGHT.PRESWING.value.data );
    ActivationIndexes(iCh).RIGHT.PRESWING.value.stdev = std( ActivationIndexes(iCh).RIGHT.PRESWING.value.data );
    ActivationIndexes(iCh).RIGHT.GAITCYCLE.value.mean = mean( ActivationIndexes(iCh).RIGHT.GAITCYCLE.value.data );
    ActivationIndexes(iCh).RIGHT.GAITCYCLE.value.stdev = std( ActivationIndexes(iCh).RIGHT.GAITCYCLE.value.data );
    
    cd(GlobalVar.OUTPUTdir)
        Outputfilename = strcat( FILENAME(1:end-8), EMG(1).Name{iCh},'_ActIndex');
        WriteYaml(Outputfilename, ActivationIndexes(iCh),0);
    cd(GlobalVar.TOOLSdir)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%         COCONTRACTION INDEXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumCouples = numel(CHforIdx.CHANNELS);

for iCp = 1 : NumCouples
    Couple = str2num( cell2mat(CHforIdx.CHANNELS{iCp} ) );
    Ch_1 = Couple(1);
    Ch_2 = Couple(2);
    
    for iStp = 1 : NumLeftSteps
        % LeftStance
        BS = floor( L_BST( iStp ) * EMG(3).sR ); % BS stands for Begin Sample
        ES = floor( L_EST( iStp ) * EMG(3).sR ) ; % ES stands for End Sample
        CocontractionIndexes(iCp).LEFT.STANCE.value.data( iStp )  = cocontraction_winter( Time(BS:ES), Data(Ch_1,BS:ES), Data(Ch_2,BS:ES)  );
        
        % LeftSwing
        BS = floor( L_BSW( iStp ) * EMG(3).sR ); % BS stands for Begin Sample
        ES = floor( L_ESW( iStp ) * EMG(3).sR ) ; % ES stands for End Sample
        CocontractionIndexes(iCp).LEFT.SWING.value.data( iStp )  = cocontraction_winter( Time(BS:ES), Data(Ch_1,BS:ES), Data(Ch_2,BS:ES)  );
        
        % LeftLoadResponse
        if L_BST(1) > R_BST(1)
            BS = floor( L_BST( iStp ) * EMG(3).sR ); % BS stands for Begin Sample
            ES = floor( R_EST( iStp ) * EMG(3).sR ) ; % ES stands for End Sample
        CocontractionIndexes(iCp).LEFT.LOADRESPONSE.value.data( iStp )  = cocontraction_winter( Time(BS:ES), Data(Ch_1,BS:ES), Data(Ch_2,BS:ES)  );
        else
            if iStp ~= NumLeftSteps
                BS = floor( L_BST( iStp + 1) * EMG(3).sR ); % BS stands for Begin Sample
                ES = floor( R_EST( iStp ) * EMG(3).sR ) ; % ES stands for End Sample
                CocontractionIndexes(iCp).LEFT.LOADRESPONSE.value.data( iStp )  = cocontraction_winter( Time(BS:ES), Data(Ch_1,BS:ES), Data(Ch_2,BS:ES)  );            
            end
        end
        
        % LeftPreSwing
        if L_BST(1) > R_BST(1)
            BS = floor( R_ESW( iStp ) * EMG(3).sR ); % BS stands for Begin Sample
            ES = floor( L_EST( iStp ) * EMG(3).sR ) ; % ES stands for End Sample
        else
            BS = floor( R_BST( iStp ) * EMG(3).sR ); % BS stands for Begin Sample
            ES = floor( L_EST( iStp ) * EMG(3).sR ) ; % ES stands for End Sample
        end
        CocontractionIndexes(iCp).LEFT.PRESWING.value.data( iStp )  = cocontraction_winter( Time(BS:ES), Data(Ch_1,BS:ES), Data(Ch_2,BS:ES)  );           
        
        % leftGaitCycle
        BS = floor( L_BST( iStp ) * EMG(3).sR ); % BS stands for Begin Sample
        ES = floor( L_ESW( iStp ) * EMG(3).sR ) ; % ES stands for End Sample
        CocontractionIndexes(iCp).LEFT.GAITCYCLE.value.data( iStp )  = cocontraction_winter( Time(BS:ES), Data(Ch_1,BS:ES), Data(Ch_2,BS:ES)  );
    end
    % left mean and standard deviations
    CocontractionIndexes(iCp).LEFT.STANCE.value.mean = mean( CocontractionIndexes(iCp).LEFT.STANCE.value.data );
    CocontractionIndexes(iCp).LEFT.STANCE.value.stdev = std( CocontractionIndexes(iCp).LEFT.STANCE.value.data );
    CocontractionIndexes(iCp).LEFT.SWING.value.mean = mean( CocontractionIndexes(iCp).LEFT.SWING.value.data );
    CocontractionIndexes(iCp).LEFT.SWING.value.stdev = std( CocontractionIndexes(iCp).LEFT.SWING.value.data );
    CocontractionIndexes(iCp).LEFT.LOADRESPONSE.value.mean = mean( CocontractionIndexes(iCp).LEFT.LOADRESPONSE.value.data );
    CocontractionIndexes(iCp).LEFT.LOADRESPONSE.value.stdev = std( CocontractionIndexes(iCp).LEFT.LOADRESPONSE.value.data );
    CocontractionIndexes(iCp).LEFT.PRESWING.value.mean = mean( CocontractionIndexes(iCp).LEFT.PRESWING.value.data );
    CocontractionIndexes(iCp).LEFT.PRESWING.value.stdev = std( CocontractionIndexes(iCp).LEFT.PRESWING.value.data );
    CocontractionIndexes(iCp).LEFT.GAITCYCLE.value.mean = mean( CocontractionIndexes(iCp).LEFT.GAITCYCLE.value.data );
    CocontractionIndexes(iCp).LEFT.GAITCYCLE.value.stdev = std( CocontractionIndexes(iCp).LEFT.GAITCYCLE.value.data );
    
    %%%
       
    for iStp = 1 : NumRightSteps
        % RightStance
        BS = floor( R_BST( iStp ) * EMG(3).sR ); % BS stands for Begin Sample
        ES = floor( R_EST( iStp ) * EMG(3).sR ) ; % ES stands for End Sample
        CocontractionIndexes(iCp).RIGHT.STANCE.value.data( iStp )  = cocontraction_winter( Time(BS:ES), Data(Ch_1,BS:ES), Data(Ch_2,BS:ES)  );
        
        % RightSwing
        BS = floor( R_BSW( iStp ) * EMG(3).sR ); % BS stands for Begin sample
        ES = floor( R_ESW( iStp ) * EMG(3).sR ) ; % ES stands for End sample
        CocontractionIndexes(iCp).RIGHT.SWING.value.data( iStp )  = cocontraction_winter( Time(BS:ES), Data(Ch_1,BS:ES), Data(Ch_2,BS:ES)  );
        
        % RightLoadResponse
        if R_BST(1) > L_BST(1)
            BS = floor( R_BST( iStp ) * EMG(3).sR ); % BS stands for Begin sample
            ES = floor( L_EST( iStp ) * EMG(3).sR ) ; % ES stands for End sample
        CocontractionIndexes(iCp).RIGHT.LOADRESPONSE.value.data( iStp )  = cocontraction_winter( Time(BS:ES), Data(Ch_1,BS:ES), Data(Ch_2,BS:ES)  );    
        else
            if iStp ~= NumRightSteps
                BS = floor( R_BST( iStp + 1) * EMG(3).sR ); % BS stands for Begin sample
                ES = floor( L_EST( iStp ) * EMG(3).sR ) ; % ES stands for End sample
                CocontractionIndexes(iCp).RIGHT.LOADRESPONSE.value.data( iStp )  = cocontraction_winter( Time(BS:ES), Data(Ch_1,BS:ES), Data(Ch_2,BS:ES)  );
            end
        end
        
        % RightPreSwing
        if R_BST(1) > L_BST(1)
            BS = floor( L_ESW( iStp ) * EMG(3).sR ); % BS stands for Begin sample
            ES = floor( R_EST( iStp ) * EMG(3).sR ) ; % ES stands for End sample
        else
            BS = floor( L_BST( iStp ) * EMG(3).sR ); % BS stands for Begin sample
            ES = floor( R_EST( iStp ) * EMG(3).sR ) ; % ES stands for End sample
        end
        CocontractionIndexes(iCp).RIGHT.PRESWING.value.data( iStp )  = cocontraction_winter( Time(BS:ES), Data(Ch_1,BS:ES), Data(Ch_2,BS:ES)  );

        % leftGaitCycle
        BS = floor( R_BST( iStp ) * EMG(3).sR ); % BS stands for Begin sample
        ES = floor( R_ESW( iStp ) * EMG(3).sR ) ; % ES stands for End sample
        CocontractionIndexes(iCp).RIGHT.GAITCYCLE.value.data( iStp )  = cocontraction_winter( Time(BS:ES), Data(Ch_1,BS:ES), Data(Ch_2,BS:ES)  );
    end
    % right mean and standard deviations
    CocontractionIndexes(iCp).RIGHT.STANCE.value.mean = mean( CocontractionIndexes(iCp).RIGHT.STANCE.value.data );
    CocontractionIndexes(iCp).RIGHT.STANCE.value.stdev = std( CocontractionIndexes(iCp).RIGHT.STANCE.value.data );
    CocontractionIndexes(iCp).RIGHT.SWING.value.mean = mean( CocontractionIndexes(iCp).RIGHT.SWING.value.data );
    CocontractionIndexes(iCp).RIGHT.SWING.value.stdev = std( CocontractionIndexes(iCp).RIGHT.SWING.value.data );
    CocontractionIndexes(iCp).RIGHT.LOADRESPONSE.value.mean = mean( CocontractionIndexes(iCp).RIGHT.LOADRESPONSE.value.data );
    CocontractionIndexes(iCp).RIGHT.LOADRESPONSE.value.stdev = std( CocontractionIndexes(iCp).RIGHT.LOADRESPONSE.value.data );
    CocontractionIndexes(iCp).RIGHT.PRESWING.value.mean = mean( CocontractionIndexes(iCp).RIGHT.PRESWING.value.data );
    CocontractionIndexes(iCp).RIGHT.PRESWING.value.stdev = std( CocontractionIndexes(iCp).RIGHT.PRESWING.value.data );
    CocontractionIndexes(iCp).RIGHT.GAITCYCLE.value.mean = mean( CocontractionIndexes(iCp).RIGHT.GAITCYCLE.value.data );
    CocontractionIndexes(iCp).RIGHT.GAITCYCLE.value.stdev = std( CocontractionIndexes(iCp).RIGHT.GAITCYCLE.value.data );
    
    cd(GlobalVar.OUTPUTdir)
        Outputfilename = strcat( FILENAME(1:end-8), EMG(1).Name{Ch_1} ,'-vs-', EMG(1).Name{Ch_2},'_CoCoIndex');
        WriteYaml(Outputfilename, CocontractionIndexes(iCp),0);
    cd(GlobalVar.TOOLSdir)
   
end

cd(GlobalVar.ROOT)

Indexes = struct();
Indexes.ActivationIndexes = ActivationIndexes;
Indexes.CocontractionIndexes = CocontractionIndexes;

end





