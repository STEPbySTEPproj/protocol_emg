function [Delsys] = LoadDelsysData(filename, Nexus, toBePlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                                           %%%%%%%%%%%%%%
%%%%%%%%%%%%%%    [Delsys] = LoadDelsysData(filename, toBePlot)          %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                                           %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                                 %%%
%%%     To load, elaborate (interpolation to eliminate zero signal) and             %%%
%%%     plot he EMG and acceleratoin signals                                        %%%
%%%                                                                                 %%%
%%%                                                                                 %%%
%%%     IN:     filename        --->    STRING                                      %%%
%%%                                     name of the CSV file with EMG and           %%%
%%%                                     acceleration data from the Delsys system    %%%
%%%                                                                                 %%%
%%%             (CHANNELS_NAME) --->    CELL with strings                           %%%
%%%                                     Name of CHANNELS Loaded from Ch_NAME.yml    %%%
%%%                                     e.g. CH_NAME: {' LEFT RF', ' LEFT TA',      %%%
%%%                                     ...}                                        %%%
%%%                                                                                 %%%
%%%             Nexus           --->    BOOL                                        %%%
%%%                                     0 acquired wih Delsys software              %%%
%%%                                     1 acquired with Nexus software              %%%
%%%                                                                                 %%%
%%%             toBePlot        --->    BOOL                                        %%%
%%%                                     0 do not plot thedata                       %%%
%%%                                     1 plot the data                             %%%
%%%                                                                                 %%%
%%%     OUT:    Delsys          --->    STRUCT                                      %%%
%%%                                     Delsys.ACC.Name acceleration channels name  %%%
%%%                                     Delsys.ACC.Data acceleration channels data  %%%
%%%                                     Delsys.ACC.sR ...                           %%%
%%%                                     (see signalElab structure)                  %%%
%%%                                     Delsys.EMG.Name EMG channels name           %%%
%%%                                     Delsys.EMG.Data EMG channels data           %%%
%%%                                                                                 %%%
%%%                                                                                 %%%
%%%                 SUBROUTINES:    ReadYaml,WriteYaml,                             %%%
%%%                                 discrete_integrate, cocontraction_winter,       %%%
%%%                                                                                 %%%
%%%                                                                                 %%%
%%%   Author:   Marco Caimmi                                                        %%%
%%%             STIIMA Nationl Research Council of Italy                            %%%
%%%             marco.caimmi@stiima.cnr.it                                          %%%
%%%             marco.caimmi@gmail.com                                              %%%
%%%                                                                                 %%%
%%%   Year:     2021                                                                %%%
%%%                                                                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	File:		View3D.m
%	Author:		Marco Caimmi
%	Release:	1.0
%	Date:		02.03.2000
%	Copyright:	(c) Department of Bioengineering - Villa Beretta
%			H Valduce Como - Italy 
%	Platform:	MATLAB 5.2.0.3084


global GlobalVar
Delsys = struct();
EMGChFilename = strcat(filename(1:end-8) , 'EMGchNames.yml');cd(GlobalVar.DATAdir)
iSFile = dir(EMGChFilename);
if isempty(iSFile)
    NAMES = ReadYaml('EMGChFile.yml');
    WriteYaml(EMGChFilename, NAMES, 1);
else
    NAMES = ReadYaml(EMGChFilename);
end

cd(GlobalVar.ROOT)

EMGChNames = NAMES.CH_NAMES;   clear NAMES;
NumEMGCh = numel(EMGChNames);
Delsys.ACC.NumCh = 6 * NumEMGCh; % num EMG  * 6 deg of an IMU

if Nexus
    cd(GlobalVar.DATAdir)
    EMG_FILENAME = strcat(filename(1:end-8) , 'emgs.csv');
    ACC_FILENAME = strcat(filename(1:end-8) , 'imus.csv');
    %
    %%%%%%%%% EMG
    %
    EMGTabData = readtable(EMG_FILENAME);
    Probes = [];
    EMG = [];
    for iProbe = 1:16
        if EMGTabData{1,iProbe+2} ~= 0
            Probes = [Probes,iProbe];
            EMG = [EMG EMGTabData{:,iProbe+2}];
        end
    end
    
    %
    %%%%%%%%% ACC
    %
    ACCTabData = readtable(ACC_FILENAME);
    cd(GlobalVar.ROOT)
    ACC =[];
    for iProbe = 1:numel(Probes)
        Begin = 3+(Probes(iProbe)-1)*9;
        End = Begin+5;
        ACC = [ACC ACCTabData{ : , Begin : End } ];
    end
    
%     ACCX = ACC(:,1:6:end);
%     ACCY = ACC(:,2:6:end)   ;
%     ACCZ = ACC(:,3:6:end);
%     GYROX1 = ACC(:,4:6:end);
%     GYROX2 = ACC(:,5:6:end);
%     GYROX3 = ACC(:,6:6:end);
    %
    %%%%%%%% TIMES
    %
    Delsys.EMG.sR = 2000;
    Delsys.ACC.sR = 200;
    
    tEMG = [0:size(EMG,1)-1] / Delsys.EMG.sR;
    tACC = [0:size(ACC,1)-1] / Delsys.ACC.sR;
    tGYR = tACC;
    
   
    Delsys.ACC.Data = ACC;
    
else % Delsys software used to acquire the data
    
    cd(GlobalVar.DATAdir)
    TabData = readtable(filename);
    cd(GlobalVar.ROOT)
    
    Data = TabData{:,:};  % Da tabella a matrice;
    
    if NumEMGCh ~= size(Data,2)/14
        msgbox('Incorrect number of channels! Verify names in CHANNELS');
        return
    end
    
    EMG  = TabData{:,2:14:end};
    ACCX = TabData{:,4:14:end};
    ACCY = TabData{:,6:14:end};
    ACCZ = TabData{:,8:14:end};
    GYROX1 = TabData{:,10:14:end};
    GYROX2 = TabData{:,12:14:end};
    GYROX3 = TabData{:,14:14:end};
    
    tEMG = TabData{:,1};
    tACC = TabData{:,3};
    tGYR = TabData{:,9};
    
    
    %%%%% The final part of the ACC signals are NaN
    %%%
    %%% Find last valid value different from NaN
    %%%
    EndEMG = min( find( isnan( tEMG(:,1) ) == 1) ) - 1;
    EndACC = min( find( isnan( tACC(:,1) ) == 1) ) - 1;
    EndGYR = min( find( isnan( tGYR(:,1) ) == 1) ) - 1;
    
    % Resizing
    tEMG = tEMG(1:EndEMG,1);
    tACC = tACC(1:EndACC,1);
    tGYR = tGYR(1:EndGYR,1);
    
    EMG = EMG(1:EndEMG,:);
    ACCX = ACCX(1:EndACC,:);
    ACCY = ACCY(1:EndACC,:);
    ACCZ = ACCZ(1:EndACC,:);
    GYROX1 = GYROX1(1:EndGYR,:);
    GYROX2 = GYROX2(1:EndGYR,:);
    GYROX3 = GYROX3(1:EndGYR,:);
    DataRaw = [ACCX ACCY ACCZ GYROX1 GYROX2 GYROX3];
    
    
    
    %%%%%% Interpolation (needed to eliminate NaN inside the signals but not
    %%%%%% the NaN at the beginning of the signals
    
    NumACCsamples = numel(tACC);
    xData = [ 1:NumACCsamples ]';
    Interpolated_data = zeros(NumACCsamples,8);
    
    NumColumns = size(DataRaw,2);
    
    for i = 1:NumColumns
        if sum(DataRaw(:,i)) ~= 0 % signal is not empty
            iZeros = find(DataRaw(:,i) == 0);
            DataRaw(iZeros,i) = NaN;
            yData = DataRaw(:,i);
            yData(iZeros) = interp1( xData( ~isnan(yData) ), yData( ~isnan(yData) ), iZeros );
            Interpolated_data(:,i) = yData;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%% Elaboration of the signals to find first valid
    %%%%%%%%%%%%%%%%%%%%%% samples
    %%%
    %%%% Find first NaNs in ACC signals
    %
    %
    ACCBeginFile = 0;
    for i = 1:Delsys.ACC.NumCh
        maxiCh = find( isnan( Interpolated_data(:,i)) == 1, 1, 'last' )+1;
        %     plot( isnan( Interpolated_data( :,Delsys.ACC.NumCh) ),'bx' )
        %     axis([0 50 0 1])
        %     pause
        ACCBeginFile = max(ACCBeginFile,maxiCh);
    end
    
    %%%% Find first zeros in EMG signals
    %
    %
    EMGBeginFile = 1;
    for iCh =1:NumEMGCh
        iBeginEMG = find( EMG(:,iCh) == 0, 1, 'last' );
        EMGBeginFile = max(EMGBeginFile,iBeginEMG);
    end
    
    TBeginEMG = tEMG(EMGBeginFile);
    TBeginACC = tACC(ACCBeginFile);
    TBegin = max(TBeginEMG,TBeginACC);
    EMGBeginFile = ceil(TBegin * numel(tEMG)/tEMG(end) );
    ACCBeginFile = ceil(TBegin * numel(tACC)/tACC(end) );
    
    %%%%%%%%%%%%% Resizing
    %
    tEMG = tEMG(EMGBeginFile:end) - tEMG(EMGBeginFile);
    EMG = EMG(EMGBeginFile:end,:);
    %
    tACC = tACC(ACCBeginFile:end) - tACC(ACCBeginFile);
    tGYR = tACC;
    Interpolated_data = Interpolated_data(ACCBeginFile:end,:);
    %%%%
    
    % Elaboration to eliminate the offet
    yy = Interpolated_data;
    NumNaN = isnan(yy);
    if sum(sum((NumNaN))) ~= 0
        msgbox('NaN in signals. Check LoadDelsysData')
        return
    end
    %%%%% substract offset
    %
    pp = fft(yy);
    n = numel(tACC);
    yy = yy - sign(pp(1,:)).*abs(pp(1,:))/n;
    %yy = yy-pp(1)/n;
    ElabData = yy;

    
    Delsys.ACC.sR = round(Delsys.ACC.NumSample/tACC(end));
    Delsys.ACC.Data = ElabData;
end %%% if Nexus

%%%% ACC DATA
%
Name = {'ACCX', 'ACCY', 'ACCZ', 'GYROX1', 'GYROX2', 'GYROX3'};
IndexCh = 0;
for  i = 1: NumEMGCh
    for j = 1: 6
        IndexCh = IndexCh + 1;
        Delsys.ACC.Name{IndexCh} = {strcat(Name{j},'-', num2str(i))};
    end
end

Units= {};
for iCh = 1:Delsys.ACC.NumCh
    Units{iCh}= 'mms-2';
end
Delsys.ACC.NumSample = length(tACC);
Delsys.ACC.Time = tACC;
Delsys.ACC.NZoom =[1 Delsys.ACC.NumSample];

%%%%%%% EMG DATA
%
%

%%%% EMG interpolation
xData = [ 1:length(tEMG) ]';
EMGRaw = EMG;
for i = NumEMGCh
    if sum(EMG(:,i)) ~= 0 % signal is not empty
        iZeros = find( abs(EMG(:,i)) < 1e-6 );
        if ~isempty(iZeros)
            EMGRaw(iZeros,i) = NaN;
            yData = EMGRaw(:,i);
            yData(iZeros) = interp1( xData( ~isnan(yData) ), yData( ~isnan(yData) ), iZeros );
            EMG(:,i) = yData;
        end
    end
end


Delsys.EMG.Name = EMGChNames;
Delsys.EMG.Data = EMG;
Delsys.EMG.NumFootSw = 0;
Delsys.EMG.NumCh = NumEMGCh;
Units= {};
for iCh = 1:NumEMGCh
    Units{iCh}= 'V';
end

Delsys.EMG.NumSample = length(tEMG);
Delsys.EMG.sR = round(Delsys.EMG.NumSample/tEMG(end));
Delsys.EMG.Time = tEMG;
Delsys.EMG.NZoom = [1 Delsys.EMG.NumSample];


if toBePlot
    
    % EMG
    figure(68)
    ChPlotNames = Delsys.EMG.Name;
    set(gcf,'Name','EMG')
    OrderCh = {1 3 5 7 2 4 6 8};
    for jCh = 1:NumEMGCh
        subplot(4,2,OrderCh{jCh});
        plot(tEMG,Delsys.EMG.Data(:, jCh) );
        title( ChPlotNames{jCh});
    end
    
    % IMUS
    ChPlotNames = Name;
    for iPlot = 1:6
        figure(iPlot +10)
        clf
            set(gcf,'Name',ChPlotNames{iPlot})
        OrderCh = {1 3 5 7 2 4 6 8};
        for jCh = 1:NumEMGCh
            subplot(4,2,OrderCh{jCh});
            titleName = strcat( ChPlotNames{iPlot},' -  ', EMGChNames{jCh} );
            plot(tACC,Delsys.ACC.Data(:,(iPlot-1) + jCh));
            title(titleName);
        end
    end
    
    
    
end


end


