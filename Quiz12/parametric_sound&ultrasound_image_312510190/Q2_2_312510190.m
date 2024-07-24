%% Quiz 2: Time-Frequency Analysis of Music. One is the original music and the other one is the sound recorded by parametric speaker (Google it!)
clear all 

%%% Read out two music (fs is the sampling frequency)
[org_music fs] = audioread('org.mp4'); %original music
[par_music fs] = audioread('par.mp4'); % parametric sound

N1 = length(org_music);  % length of the data
N2 = length(par_music);  %  N1 is the same as N2
t_indx = [0:1:N1-1]*1/fs;  % time index in second


%%% Question 2.2
%%% Time-frequency analysis (Complete all codes and display the STFT results)  
%%% Note1: You can use MATLAB "stft" function or write
%%% STFT by yourself, but you need to specific the following four
%%% parameters: Window size, Windowing function, FFT length, Overlapping
%%% length.
       
%%% Note2: Because the analysis is time comsuming, you can crop your data
%%% to be a length of 220000 (roughly 5 seconds).

%%% Note3: Display your STFT result in dB within the frequency range of 0-4000 Hz. 
org_part = org_music(1:220000);  % crop the data to be 5 seconds
par_part = par_music(1:220000);  % crop the data to be 5 seconds

w_len = 0;  % window size for STFT
win = 0 ;  % Windowing function  
fft_len = 0;  % FFT length for STFT (can be different from the window size).
OverlapLength =0;  % Overlapping window length between two STFT analses.  

[s,f,t] = stft(org_part,fs,Window=hamming(1024,'periodic'),OverlapLength=512,FFTLength=1024);
figure
sdb = mag2db(abs(s));
mesh(t,f/1000,sdb);

cc = max(sdb(:))+[-60 0];
ax = gca;
ax.CLim = cc;
view(2)

xlim([0 5])  % confine the display time range
ylim([0 4])  % confine the display frequency range
xlabel('time (s)')
ylabel('Frequency (kHz)')
title('Spectrogram of orginal music')
colorbar

[s1,f1,t1] = stft(par_part,fs,Window=hamming(1024,'periodic'),OverlapLength=512,FFTLength=1024);
figure
sdb1 = mag2db(abs(s1));
mesh(t1,f1/1000,sdb1);

cc1 = max(sdb(:))+[-60 0];
ax = gca;
ax.CLim = cc1;
view(2)

xlim([0 5])  % confine the display time range
ylim([0 4])  % confine the display frequency range
xlabel('time (s)')
ylabel('Frequency (kHz)')
title('Spectrogram of parametric music')
colorbar

%%% Question 2.3.  Improve the parametric sound by filtering. Complete your
%%% design below, and output the sound as a mp4 file using the following
%%% function.

%audiowrite('par_improvement.mp4',par_music,fs);


