%%%% Quiz 1 (ADSP 2024)
%%%% We will make an image from the raw 96-channel data collected by a 96-element phased array. The beamforming technique used here is synthetic aperture.  
%%%% You need to complete the beamforming and envelope detection to obtain the final B-mode image.
%%%% Created by Geng-Shi Jeng 5/21/2024

close all
clear all

load PA_chdata      %% Load the channel data and the associated parameters. 
%==============================

%data;  %Ultrasound full dataset collected by each Tx and Rx combination.
       %4096 x 96 x 96 matrix size. 4096 means number of the range samples. The first 96 means the received channel index whereas the 2nd 96 means the transmitted channel index.
[Nz, Nx, Nt] = size (data);  % Nz=4096 (Samples), Nx = 96 (Rx), Nt = 96 (Tx)

fc = B_fo;               %  center frequency [Hz]
c = SV;                  %  Speed of sound [m/s]
fs;                      % sampling rate or frequency  (Hz)
toff = B_Toffset;        % the starting time for channel data  (s)

%%%% Array coordinate
chx =-1*(Nx-1)/2*pitch:pitch:(Nx-1)/2*pitch; % x-coordinate of a phased array. The element spacing (i.e., pitch) is the variable called "pitch". The number of elements is Nx (or Nt)  (m)
chz = zeros(1,Nx);   %z-coordinate of a phased array (m)

%%%% Paramters for Field of View (FOV) and scan lines
FOV_angle = 90*pi/180; % The angle of FOV (in radian).  
d_th = sin(0.9/1.5/180*pi);  % the angle increment between the successive scan lines.
no_lines = round(2*sin(FOV_angle/2)/d_th);  % determining the number of scan lines
if mod(no_lines,2) == 0  % Forcing no_lines to be odd.
   no_lines = no_lines+1;             
end

%%% Depth information for B-mode
dz = 1/(4*fc)*c/2;        % range sampling interval (m)
range = 0:dz:80*10^(-3);  % The sampling range points for beamforming  (m)
Nr = max(size(range));  % number of the sampling range points

% ======================== Synthetic aperture beamforming =============================
% interpolation on acquired channel data
interp_num = 1;     %%%-------------------> You can play with this value to see if the imaging quality varies.

%%% high-pass filter to remove the unwanted signals of the channel data around DC.  
hpf = fir1(32,0.2,'high');

%%%% Synthetic aperture beamforming
bf = zeros(Nr,no_lines); % initialization of final beamformation result 

% Reduce the number of Tx events by using every alternate Tx channel
Tx_channels = 1:6:Nt;

for ii= Tx_channels  % for every ii transmit element
ii
    %%% Read the raw data for each ii transmit element
    bb_data1 =data(:,:,ii);   
    [Nz, Nx] = size(bb_data1);

    bb_data = conv2(hpf.',1,bb_data1,'same');  % high-pass filter to remove the unwanted signals of the channel data around DC.  
    bb_data = resample(bb_data,interp_num,1);  %interpolation of the channel data in the range direction by a factor of interp_num
    
    %%% Initialization for each sub-image associated with ii-th Tx emission.
    subbf = zeros(Nr, no_lines);

    for jj = 1: no_lines  %% for each scan line
        
        % The scan line position
         sin_th = -1*(no_lines-1)/2*d_th  + (jj-1)*d_th;
         cos_th = cos(asin(sin_th));

         %%%% Delay and Sum (DAS) to realize the beamforming
         % Rx distance  (You should calculate the distance between each image pixel to ALL receive channels)
          Rx_dist = zeros(Nr,Nx);           % initialization of the Rx distance
         % Rx_dist =                                                ----------------------------------------------->   You  need complete the code here 
         [R, C] = meshgrid(range, chx); % Create a meshgrid for range and channel positions
         Rx_dist = sqrt((R*cos_th).^2 + (C - R * sin_th).^2).'; % Euclidean distance from each point to the Rx elements

         % Tx distance (You should calculate the distance between each image pixel to ii-th transmit channel)
          Tx_dist = zeros(Nr,Nx);           % initialization of the Tx distance
         % Tx_dist =                                                ----------------------------------------------->   You  need complete the code here 
         Tx_dist = sqrt((R*cos_th).^2 + (chx(ii) - R * sin_th).^2).'; % Euclidean distance from each point to the Tx element

         % Calculate the total delay for ii-th Tx channel to ALL Rx channels.
         totdelay = (Rx_dist + Tx_dist)/c-toff;  % total delay 
         delayindx = round(totdelay*fs*interp_num)+1;  % convert the total delay to the sample index.
          
         %%% Find the sample according to "delayindx" and do the summation 
         % subbf(:,jj) =   % sub-image for ii-th Tx emission       ----------------------------------------------->   You  need complete the code here
         for k = 1:Nx
            valid_idx = delayindx(:, k) <= size(bb_data, 1); % Ensure indices are within bounds
            subbf(valid_idx, jj) = subbf(valid_idx, jj) + bb_data(delayindx(valid_idx, k), k);
        end
         
    end
    bf = bf + subbf;   % summation among all sub-images.
    
end

%%% Envelope detection 
bf_env = zeros(Nr,no_lines); % initialization of envelope detected result
%bf_env =             % envelope detected result               ----------------------------------------------->   You  need complete the code here
bf_env =  abs(hilbert(bf)); % envelope detected result 

%%% Log compression
bf_env_nor = abs(bf_env)/max(max(abs(bf_env)));  % normalized to 1.
srcInt = 20*log10(bf_env_nor);

%%% Log compression of the beamformed image
bf_nor = abs(bf)/max(max(abs(bf)));  % normalized to 1.
bf_srcInt = 20*log10(bf_nor);

%%% Display the result
DR=60;   % dynamic range

%%% scan conversion (From rectangle to polar format)
sin_th = -1*(no_lines-1)/2*d_th + ([1:1:no_lines]-1)*d_th; 
cos_th = cos(asin(sin_th));
Z=range.'*cos_th;
X=range.'*sin_th;

figure
pcolor(X*1e3, Z*1e3, bf_srcInt+DR )
axis('ij');
shading interp;
caxis([0 DR])
colorbar
colormap gray
xlabel('Lateral distance [mm]')
ylabel('Axial distance [mm]')
title('Phased array imaging without envelope detection')
axis equal

figure
pcolor(X*1e3, Z*1e3, srcInt+DR )
axis('ij');
shading interp;
caxis([0 DR])
colorbar
colormap gray
xlabel('Lateral distance [mm]')
ylabel('Axial distance [mm]')
title('Phased array imaging(emmision every 6 Tx channels) ')
axis equal



















