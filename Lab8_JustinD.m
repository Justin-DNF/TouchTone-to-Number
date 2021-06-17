% File: Lab8_JustinD.m
% Name: Justin D
% Class: ECEGR 3711 01 LAB 3 Signals and Systems
% Lab: #08
% Lab Date: 05/26/2020
% Due Date: 06/02/2020

%% LAB 8

clear all; % removes all variables, globals, functions and MEX links
clc; % clear Command Window

% VARIABLES

% Define parameters of a signal
frequency = 100.1; % Hz
amplitude = 1; %  Voltage
phase = 0; % (0-2*pi) radians
offset = 0;

%Define Plot limits for x-axis
xlimits = [90 110]; 

% Define sample rate 
sampleRate = 1000; % Hz

% Define length of sample
totaltime = 2; %Seconds

% Sample Sinusoidal wave 
% Must have even number of values
t = 0:(1/(sampleRate)):(totaltime -(1/(sampleRate)));

% Create signal
signal = amplitude*cos((frequency*t*2*pi)+phase)+ offset;

% Apply Windowing
%signalWindow = signal.*boxcar(length(t))';
%signalWindow = signal.*hanning(length(t))';
%signalWindow = signal.*flattopwin(length(t))';

%{

% Setup 'Sinusoidal Plot'
f1 = figure(1);
f1.Name = 'Sinusoidal Plot';
clf % Clear current figure

% Subplot 1:  Individual signals
subplot(3,1,1)
hold on % Hold the current figure on
ylabel('x(t)')
grid on
box on
%Plot the sinusoidal 
plot(t,signal,'b-')
% Plot the zero line
plot(t,(t.*0),'k')
hold off

% Calculate Fast Fourier Transform of Signal
Y = fft(signal);
YWindow = fft(signalWindow);
L = length(Y);

% Calculate Frequency Array
freq = linspace(sampleRate/L,sampleRate/2,L/2);

% Calculate Amplitude Spectra for Pure Signal
Amp2 = abs(Y/L);
Amp1 = 2*Amp2(2:L/2+1);
dBAmp1=20*log10(Amp1);

% Calculate Amplitude Spectra for Windowed Signal
Amp2Window = abs(YWindow/L);
Amp1Window = 2*Amp2Window(2:L/2+1);
dBAmp1Window=20*log10(Amp1Window);


% Subplot 2:  Amplitude Spectra
subplot(3,1,2)
hold on
stem(freq,Amp1,'r-')
stem(freq,Amp1Window,'b-')
xlabel('Frequency')
ylabel('Amplitude')
xlim(xlimits)
hold off

% Subplot 3:  Amplitude Spectra Decibel Scale
subplot(3,1,3)
hold on
plot(freq,dBAmp1,'r-')
plot(freq,dBAmp1Window,'b-')
xlabel('Frequency')
ylabel('Decibels')
xlim(xlimits)
hold off

Pos = [800 200 800 800];
set(gcf, 'Position',  Pos)
%}


% Read infomation for audio file
information = audioinfo('justin.wav');
SampleRate = information.SampleRate;
NumChannels = information.NumChannels;
[Y,Fs] = audioread('justin.wav');
size(Y);
Y = Y';

% The 8th Signal 
signal8 = Y(1, 1+SampleRate*2.1 : SampleRate*2.4);
t8 = 2.1:(1/SampleRate):(2.4 - (1/SampleRate));
signalWindow = signal8.*hanning(length(t8))';


% Setup 'Sinusoidal Plot'
f1 = figure(1);
f1.Name = 'Sinusoidal Plot';
clf % Clear current figure

% Subplot 1:  Individual signals
subplot(3,1,1)
hold on % Hold the current figure on
ylabel('x(t)')
grid on
box on
%Plot the sinusoidal 
plot(t8,signal8,'b-')
% Plot the zero line
plot(t8,(t8.*0),'k')
hold off

% Calculate Fast Fourier Transform of Signal
Y8 = fft(signal8);
YWindow_8 = fft(signalWindow);
L_8 = length(Y8);

% Calculate Frequency Array
freq_8 = linspace(SampleRate/L_8,SampleRate/2,L_8/2);

% Calculate Amplitude Spectra
Amp2_8 = abs(Y8/L_8);  % Divide the L and find the magnitude
Amp1_8 = 2*Amp2_8(2:L_8/2+1); % Parse the array and multiply by 2
dBAmp1_8=20*log10(Amp1_8);

% Calculate Amplitude Spectra for Windowed Signal
Amp2Window_8 = abs(YWindow_8/L_8);
Amp1Window_8 = 2*Amp2Window_8(2:L_8/2+1);
dBAmp1Window_8 =20*log10(Amp1Window_8);


% Subplot 2:  Amplitude Spectra
subplot(3,1,2)
hold on
stem(freq_8,Amp1_8,'r-')
stem(freq_8,Amp1Window_8,'b-')
xlabel('Frequency')
ylabel('Amplitude')
hold off

% Subplot 3:  Amplitude Spectra Decibel Scale
subplot(3,1,3)
hold on
plot(freq_8,dBAmp1_8,'r-')
plot(freq_8,dBAmp1Window_8,'b-')
xlabel('Frequency')
ylabel('Decibels')
hold off

Pos = [800 200 800 800];
set(gcf, 'Position',  Pos)