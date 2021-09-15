function EnvepData = EnvelopeHilbert(Data,sR)

%%%%%%%%%%%%%%%%%% EnvepData = function EnvelopeHilbert(Data,sR)
%%%
%%%
%%%%%%%%% Caculates the envelope of data using the Hilbert method  
%%%
%%%
%%%%%%%%%   IN:  ----> Data         ----> Array 1 x N 
%%%
%%%
%%%%%%%%    OUT:  ----> EnvepData  ---->  Vector with evenvelope
%%
%%
%%%%%%%    By Marco Caimmi october 2020  %%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% commento:
% low pass filter design : FIR- EQUIRIPPLE
% All frequency values are in Hz.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters have been calculated using fdatool

Fpass = 5;               % Passband Frequency
Fstop = 20;              % Stopband Frequency
Dpass = 0.057501127785;  % Passband Ripple
Dstop = 0.0001;          % Stopband Attenuation
dens  = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(sR/2), [1 0], [Dpass, Dstop]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);

OriginalData = Data;
S = size(Data);
if S(1) < S(2)
    Data = Data';   
end

%%%%% Elaboration
Data = Data - mean(Data);
Data = abs( hilbert(Data) );
Data1 = filter(Hd, Data);
Dataf1Flip = flipud(Data1);
Dataf2Filp = filter(Hd, Dataf1Flip);

EnvepData = ( flipud(Dataf2Filp) );

if size(EnvepData) ~= size(OriginalData)
    EnvepData = EnvepData';   
end

