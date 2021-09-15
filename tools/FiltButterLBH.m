function [DataFilt,Time,A,B] = FiltButterLBH(Data,Order,Fc,Fs,Type);

%%%function [DataFilt,Time,A,B] = FiltButterLBH(Data,Order,Fc,Fs,Type);
%
%  Digital Butterworth filter
%
%
% INPUT:    Data -->  Data matrix 
%           Order --> Filter order 
%           Fstop --> cuto-off frequency (for BandPass Fc=[Fc1 F2])
%           Fs --> sampling frequency
%           Type --> 'L' ---> 'Low'
%                    'Bs' ---> 'Band Stop'        note  Fc=[Fc1 Fc2]
%                    'Bp' ---> 'Band Pass'        note  Fc=[Fc1 Fc2]
%                    'H' --->  'High'
%
%%%
%         
% OUTPUT:   DataFilt --> Filtered Data (column)
%           Time     --> Time vector
%
%
%%%    By Marco Caimmi September 2008   %%%
%
%
OriginalData=Data;

if size(Data,1) < size(Data,2)
    Data = Data';
end

len = max(size(Data));
Time = [0:len-1]/Fs; 
Time = Time(:);

if upper(Type) == 'L'
    [B,A] = butter(Order,Fc*2/Fs);

elseif upper(Type) == 'BS'
    [B,A] = butter(Order,Fc*2/Fs,'stop');
    
elseif upper(Type) == 'BP'
    [B,A] = butter(Order,Fc(1)*2/Fs,'high');
    Data=filtfilt(B,A,Data); % trasfomazione in colonna e filtraggio
    [B,A] = butter(Order,Fc(2)*2/Fs);
    
elseif upper(Type) == 'H'
    [B,A] = butter(Order,Fc*2/Fs,'high');
  
else
    msgbox('Invalid Type');
end

DataFilt=filtfilt(B,A,Data); % trasfomazione in colonna e filtraggio


if size(OriginalData,1) ~= size(Data,1)
DataFilt = DataFilt';
Time=Time(:)';
end
