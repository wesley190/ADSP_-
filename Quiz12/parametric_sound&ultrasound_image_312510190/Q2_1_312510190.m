%% Quiz 2: Time-Frequency Analysis of Music. One is the original music and the other one is the sound recorded by parametric speaker (Google it!)
clear all 

%%% Read out two music (fs is the sampling frequency)
[org_music fs] = audioread('org.mp4'); %original music
[par_music fs] = audioread('par.mp4'); % parametric sound

N1 = length(org_music);  % length of the data
N2 = length(par_music);  %  N1 is the same as N2
t_indx = [0:1:N1-1]*1/fs;  % time index in second

%%% Question 2.1: Display waveform and spectrum
figure
plot(t_indx, org_music)
hold on
plot(t_indx, par_music)
xlabel('Time (s)')
ylabel('Magnitude')
title('Time-domain waveform comparison')
legend('Original music','Parametric sound')

%%% FFT of both music 
f_org = fft(org_music);
f_par = fft(par_music);

figure
plot(([0:1:N1-1]/N1-0.5)*fs,20*log10(abs(fftshift(f_org))))
hold on
plot(([0:1:N1-1]/N1-0.5)*fs,20*log10(abs(fftshift(f_par))))
xlim([0 2000])  % confine the display frequency range
ylim([-20 80])  % confine the display magnitude
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Spectrum comparision')
legend('Original music','Parametric sound')


%%% Question 2.2
%%% Time-frequency analysis (Complete all codes and display the STFT results)  
%%% Note1: You can use MATLAB "stft" function or write
%%% STFT by yourself, but you need to specific the following four
%%% parameters: Window size, Windowing function, FFT length, Overlapping
%%% length.
       
%%% Note2: Because the analysis is time comsuming, you can crop your data
%%% to be a length of 220000 (roughly 5 seconds).

%%% Note3: Display your STFT result in dB within the frequency range of 0-4000 Hz. 


w_len = 0;  % window size for STFT
win = 0 ;  % Windowing function  
fft_len = 0;  % FFT length for STFT (can be different from the window size).
OverlapLength =0;  % Overlapping window length between two STFT analses.  


%%% Question 2.3.  Improve the parametric sound by filtering. Complete your
%%% design below, and output the sound as a mp4 file using the following
%%% function.

audiowrite('par_improvement.mp4',par_music,fs);


