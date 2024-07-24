%% Quiz 2: Time-Frequency Analysis of Music. One is the original music and the other one is the sound recorded by parametric speaker (Google it!)
clear all 

%%% Read out two music (fs is the sampling frequency)
[org_music fs] = audioread('org.mp4'); %original music
[par_music fs] = audioread('par.mp4'); % parametric sound

N1 = length(org_music);  % length of the data
N2 = length(par_music);  %  N1 is the same as N2
t_indx = [0:1:N1-1]*1/fs;  % time index in second
%PI = 3.14159265358979323846;

%%% Question 2.1: Display waveform and spectrum
%figure
%plot(t_indx, org_music)
%hold on
%plot(t_indx, par_music)
%xlabel('Time (s)')
%ylabel('Magnitude')
%title('Time-domain waveform comparison')
%legend('Original music','Parametric sound')
%
%%%% FFT of both music 
%f_org = fft(org_music);
%f_par = fft(par_music);
%
%figure
%plot(([0:1:N1-1]/N1-0.5)*fs,20*log10(abs(fftshift(f_org))))
%hold on
%plot(([0:1:N1-1]/N1-0.5)*fs,20*log10(abs(fftshift(f_par))))
%xlim([0 8000])  % confine the display frequency range
%ylim([-20 80])  % confine the display magnitude
%xlabel('Frequency (Hz)')
%ylabel('Magnitude (dB)')
%title('Spectrum comparision')
%legend('Original music','Parametric sound')


%%% Question 2.2
%%% Time-frequency analysis (Complete all codes and display the STFT results)  
%%% Note1: You can use MATLAB "stft" function or write
%%% STFT by yourself, but you need to specific the following four
%%% parameters: Window size, Windowing function, FFT length, Overlapping
%%% length.
       
%%% Note2: Because the analysis is time comsuming, you can crop your data
%%% to be a length of 220000 (roughly 5 seconds).

%%% Note3: Display your STFT result in dB within the frequency range of 0-4000 Hz. 
%org_part = org_music(1:220000);  % crop the data to be 5 seconds
%par_part = par_music(1:220000);  % crop the data to be 5 seconds

%w_len = 0;  % window size for STFT
%win = 0 ;  % Windowing function  
%fft_len = 0;  % FFT length for STFT (can be different from the window size).
%OverlapLength =0;  % Overlapping window length between two STFT analses.  
%
%[s,f,t] = stft(org_music,fs,Window=hamming(1024,'periodic'),OverlapLength=512,FFTLength=1024);
%figure
%sdb = mag2db(abs(s));
%mesh(t,f/1000,sdb);
%
%cc = max(sdb(:))+[-60 0];
%ax = gca;
%ax.CLim = cc;
%view(2)
%
%%xlim([0 8])  % confine the display time range
%ylim([0 20])  % confine the display frequency range
%xlabel('time (s)')
%ylabel('Frequency (kHz)')
%title('Spectrogram of orginal music')
%colorbar

%[s1,f1,t1] = stft(par_part,fs,Window=hamming(1024,'periodic'),OverlapLength=512,FFTLength=1024);
%figure
%sdb1 = mag2db(abs(s1));
%mesh(t1,f1/1000,sdb1);
%
%cc1 = max(sdb(:))+[-60 0];
%ax = gca;
%ax.CLim = cc1;
%view(2)
%
%xlim([0 5])  % confine the display time range
%ylim([0 4])  % confine the display frequency range
%xlabel('time (s)')
%ylabel('Frequency (kHz)')
%title('Spectrogram of parametric music')
%colorbar

%%% Question 2.3.  Improve the parametric sound by filtering. Complete your
%%% design below, and output the sound as a mp4 file using the following
%%% function.
% Design a band-pass filter



%function out = band_pass_filter(m, n)
%
%    wh = 2*PI*4000/fs;
%    wl = 2*PI*1/fs;
%    if n == m
%        out = 2*(wh/PI - wl/PI);
%    else
%        out = 2*(sin(wh*(n-m))-sin(wl*(n-m)))/PI/(n-m) * hamming(2*m+1, n);
%    end
%end
%
%function out = hamming(N, n)
%
%    out = 0.54 - 0.46*cos(2*PI*n/(N-1));
%end

%Fl = 1; 
%Fh = 8000; 
%Wc = [Fl Fh].*(2/fs);  
%bbp = fir1(1024, Wc, 'band');
%%figure
%%plot (1:length(bbp), bbp)
%%figure
%%freqz(bbp,1)
%par_improved = filter(bbp,1,par_music);
%
%[s,f,t] = stft(par_improved,fs,Window=hamming(1024,'periodic'),OverlapLength=512,FFTLength=1024);
%figure
%sdb = mag2db(abs(s));
%mesh(t,f/1000,sdb);
%
%cc = max(sdb(:))+[-60 0];
%ax = gca;
%ax.CLim = cc;
%view(2)
%
%%xlim([0 8])  % confine the display time range
%ylim([0 15])  % confine the display frequency range
%xlabel('time (s)')
%ylabel('Frequency (kHz)')
%title('Spectrogram of parametric improvement music')
%colorbar

%-----------------
% Design a high-pass filter by using chebwin window(length=1025, 60dB attenuation)
%-----------------
Fl = 50; 
Fh = 1000; 
%Wc = [Fl Fh].*(2/fs);
Wc = Fl .*(2/fs);    
bbpL = fir1(1024, Wc, 'high', chebwin(1025, 60));
%figure
%plot (1:length(bbp), bbp)
%figure
%freqz(bbp,1)
par_improved1 = filter(bbpL,1,par_music);

%-----------------
% Design a band-pass filter to increase gain by 4 times
%-----------------
Fl = 100; 
Fh = 800; 
Wc = [Fl Fh].*(2/fs);
bbpL = fir1(1024, Wc, 'band');
%figure
%plot (1:length(bbp), bbp)
%figure
%freqz(bbp,1)
par_improvedL = filter(bbpL*4,1,par_improved1);

%-----------------
% Design a band-stop filter to remove the noise
%-----------------
bbpR = fir1(1024, Wc, 'stop');
%figure
%plot (1:length(bbp), bbp)
%figure
%freqz(bbp,1)
par_improvedR = filter(bbpR,1,par_improved1);

%-----------------
% Display the spectrogram of the improved sound
%-----------------
[s,f,t] = stft(par_improvedL,fs,Window=hamming(1024,'periodic'),OverlapLength=512,FFTLength=1024);
figure
sdb = mag2db(abs(s));
mesh(t,f/1000,sdb);

cc = max(sdb(:))+[-60 0];
ax = gca;
ax.CLim = cc;
view(2)

xlim([0 8])  % confine the display time range
ylim([0 8])  % confine the display frequency range
xlabel('time (s)')
ylabel('Frequency (kHz)')
title('Spectrogram of parametricL music')
colorbar

[s,f,t] = stft(par_improvedR,fs,Window=hamming(1024,'periodic'),OverlapLength=512,FFTLength=1024);
figure
sdb = mag2db(abs(s));
mesh(t,f/1000,sdb);

cc = max(sdb(:))+[-60 0];
ax = gca;
ax.CLim = cc;
view(2)

xlim([0 8])  % confine the display time range
ylim([0 8])  % confine the display frequency range
xlabel('time (s)')
ylabel('Frequency (kHz)')
title('Spectrogram of parametricR music')
colorbar

par_improved = par_improvedL + par_improvedR;

[s,f,t] = stft(par_improved,fs,Window=hamming(1024,'periodic'),OverlapLength=512,FFTLength=1024);
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
title('Spectrogram of parametric improvement music')
colorbar

figure
plot(t_indx, org_music)
hold on
plot(t_indx, par_improved)
xlabel('Time (s)')
ylabel('Magnitude')
title('Time-domain waveform comparison')
legend('Original music','Parametric improved sound')


%%% FFT of both music 
f_org = fft(org_music);
f_par = fft(par_improved);

figure
plot(([0:1:N1-1]/N1-0.5)*fs,20*log10(abs(fftshift(f_org))))
hold on
plot(([0:1:N1-1]/N1-0.5)*fs,20*log10(abs(fftshift(f_par))))
xlim([0 2000])  % confine the display frequency range
ylim([-20 80])  % confine the display magnitude
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Spectrum comparision')
legend('Original music','Parametric improved sound')

% Output the improved sound as an mp4 file
audiowrite('par_improvement.mp4',par_improved,fs);



