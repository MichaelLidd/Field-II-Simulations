%% field.m
%  Start up the Field simulation system
%
%  Version 1.0, April 2, 1998, JAJ

path(path, '/Users/michaelliddelow/Documents/MATLAB/VeinTech/Field_II_Simulation/Field_II_ver_3_30_mac/cfm_image_example')

field_init(0)

%% make_sct.m
%  Make a simulation of the received rf signal from
%  flow with a parabolic velocity profile
%  The result is stored as scatterer positions in a file
%
%  Version 2.0, 2/4-98, JAJ

%  Set the seed of the random number generator

randn('seed',sum(100*clock))

%  Initialize the ranges for the scatteres
%  Notice that the coordinates are in meters

R=0.005;         %  Radius of blood vessel [m]
x_range=0.08;    %  x range for the scatterers  [m]
y_range=2*R;     %  y range for the scatterers  [m]
z_range=2*R;     %  z range for the scatterers  [m]
z_offset=0.06;   %  Offset of the mid-point of the scatterers [m]

%  Set the number of scatterers. It should be roughly
%  10 scatterers per cubic wavelength

c=1540;    %  Ultrasound propagation velocity [m/s]
f0=3e6;    %  Center frequency of transducer  [Hz]
lambda=c/f0;
N=round(10*x_range/(5*lambda)*y_range/(5*lambda)*z_range/(lambda*2));
disp([num2str(N),' Scatterers'])

%  Generate the coordinates and amplitude
%  Coordinates are rectangular within the range.
%  The amplitude has a Gaussian distribution.

x=x_range*(rand(1,N)-0.5);
y=y_range*(rand(1,N)-0.5);
z=z_range*(rand(1,N)-0.5);

%  Find which scatterers that lie within the blood vessel

r=(y.^2+z.^2).^0.5;
within_vessel= r < R;

%  Assign an amplitude and a velocity for each scatterer

v0=1;   %  Largest velocity of scatterers [m/s]
velocity=v0*(1-(r/R).^2).*within_vessel;

blood_to_stationary= 10;   %  Ratio between amplitude of blood to stationary tissue
amp=randn(1,N).*((1-within_vessel) + within_vessel*blood_to_stationary);
amp=amp';

%  Generate files for the scatteres over a number of pulse emissions

Tprf=1/10e3;  %  Time between pulse emissions  [s]
Nshoots=10;   %  Number of shoots

for i=1:Nshoots

  %  Generate the rotated and offset block of sample

  theta=45/180*pi;
  xnew=x*cos(theta)+z*sin(theta);
  znew=z*cos(theta)-x*sin(theta) + z_offset;
  positions=[xnew; y; znew;]';

  %   Save the matrix with the values

  cmd = ['save sim_flow/scat_',num2str(i),'.mat positions amp'];
  eval(cmd)

  %  Propagate the scatterers and aliaze them
  %  to lie within the correct range
  
  x1=x;
  x=x + velocity*Tprf;
  outside_range= (x > x_range/2);
  x=x - x_range*outside_range;
  end
  

%% make_image.m

%  Compress the data to show 40 dB of
%  dynamic range for the CFM phantom image
%
%  version 1.2 by Joergen Arendt Jensen, April 6, 1997.

f0=3.5e6;                 %  Transducer center frequency [Hz]
fs=100e6;                 %  Sampling frequency [Hz]
c=1540;                   %  Speed of sound [m/s]
no_lines=50;              %  Number of lines in image
d_x=40/1000/no_lines;     %  Increment for image

%  Read the data and adjust it in time 

min_sample=0;
for i=1:no_lines

  %  Load the result

  cmd=['load sim_bmd/rf_ln',num2str(i),'.mat'];
  eval(cmd)
  
  %  Find the envelope
  
  if (tstart>0)
    rf_env=abs(hilbert([zeros(fix(tstart*fs-min_sample),1); rf_data]));
  else
    rf_env=abs(hilbert( rf_data( abs(tstart*fs):max(size(rf_data)) ) ));
    end
  env(1:max(size(rf_env)),i)=rf_env;
  end

%  Do logarithmic compression

D=20;   %  Sampling frequency decimation factor

log_env=env(1:D:max(size(env)),:)/max(max(env));
log_env=log(log_env+0.01);
log_env=log_env-min(min(log_env));
log_env=64*log_env/max(max(log_env));

%  Make an interpolated image

ID_bmode=10;
[n,m]=size(log_env);
new_env=zeros(n,m*ID_bmode);
for i=1:n
  if(rem(i,100) == 0)
    i;
    end
  new_env(i,:)=abs(interp(log_env(i,:),ID_bmode));
  end
[n,m]=size(new_env);
  
fn_bmode=fs/D;
clf
image(((1:(ID_bmode*no_lines-1))*d_x/ID_bmode-no_lines*d_x/2)*1000,((1:n)/fn_bmode+min_sample/fs)*1540/2*1000,new_env)
xlabel('Lateral distance [mm]')
ylabel('Axial distance [mm]')
colormap(gray(64))
%brighten(-0.35)
%axis([-20 20 30 90])
axis('image')

%% sim_flow.m
% 
% %  Example of use of the new Field II program running under Matlab
% %
% %  This example shows how a linear array B-mode system scans an image
% %  when doing color flow mapping
% %
% %  This script assumes that the field_init procedure has been called
% %  Here the field simulation is performed and the data is stored
% %  in rf-files; one for each rf-line done. The data must then
% %  subsequently be processed to yield the image. The data for the
% %  scatteres are read from the file pht_data.mat, so that the procedure
% %  can be started again or run for a number of workstations.
% %
% %  Example by Joergen Arendt Jensen and Peter Munk, March 14, 1997.
% %  Version 2.1 by Joergen Arendt Jensen, August 14, 1998 for Matlab 5.
% 
% %  Generate the transducer apertures for send and receive
% 
% f0=5e6;                  %  Transducer center frequency [Hz]
% fs=100e6;                %  Sampling frequency [Hz]
% c=1540;                  %  Speed of sound [m/s]
% lambda=c/f0;             %  Wavelength [m]
% width=lambda;            %  Width of element
% element_height=5/1000;   %  Height of element [m]
% kerf=0.05/1000;          %  Kerf [m]
% focus=[0 0 70]/1000;     %  Fixed focal point [m]
% N_elements=196;          %  Number of physical elements
% rec_N_active=64;         %  Number of active elements in receive
% xmit_N_active=64;        %  Number of active elements in transmit
% 
% % Use triangles
% 
% set_field('use_triangles',0);
% 
% %  Set the sampling frequency
% 
% set_sampling(fs);
% 
% %  Generate aperture for emission
% 
% xmit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 10,focus);
% 
% %  Set the impulse response and excitation of the xmit aperture
% 
% impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
% impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
% xdc_impulse (xmit_aperture, impulse_response);
% 
% excitation=sin(2*pi*f0*(0:1/fs:2/f0));
% xdc_excitation (xmit_aperture, excitation);
% 
% %  Generate aperture for reception
% 
% receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 10,focus);
% 
% %  Set the impulse response for the receive aperture
% 
% xdc_impulse (receive_aperture, impulse_response);
% 
% %  Do for the number of CFM lines
% 
% Ncfm=10;
% for k=1:Ncfm
% 
% %   Load the computer phantom
% 
%   cmd=['load sim_flow/scat_',num2str(k),'.mat']
%   eval(cmd)
% 
%   %   Do linear array imaging
% 
%   no_lines=20;                    %  Number of lines in image
%   image_width=40/1000;           %  Size of image sector
%   d_x=image_width/(no_lines-1);   %  Increment for image
% 
%   %  Set the different focal zones for reception
% 
%   rec_zone_start=30/1000;
%   rec_zone_stop=100/1000;
%   rec_zone_size=10/1000;
% 
%   focal_zones_center=[rec_zone_start:rec_zone_size:rec_zone_stop]';
%   focal_zones=focal_zones_center-0.5*rec_zone_size;
%   Nf=max(size(focal_zones));
%   focus_times=focal_zones/1540;
% 
%   %  Set a Hanning apodization on the receive aperture
%   %  Dynamic opening aperture is used.
% 
%   Fnumber=2.0;
%   rec_N_active_dyn=round(focal_zones_center./(Fnumber*(width+kerf)));
% 
%   for ii=1:Nf
%     if rec_N_active_dyn(ii)>rec_N_active 
%       rec_N_active_dyn(ii)=rec_N_active; 
%       end
%     rec_N_pre_dyn(ii) = ceil(rec_N_active/2  - rec_N_active_dyn(ii)/2);
%     rec_N_post_dyn(ii) = rec_N_active - rec_N_pre_dyn(ii) - rec_N_active_dyn(ii);
%     rec_apo=(ones(1,rec_N_active_dyn(ii)));
%     rec_apo_matrix_sub(ii,:)=[zeros(1,rec_N_pre_dyn(ii)) rec_apo zeros(1,rec_N_post_dyn(ii))];
%     end
% 
%   %  Transmit focus
% 
%   z_focus=40/1000;   
% 
%   %   Set a Hanning apodization on the xmit aperture
% 
%   xmit_apo=hanning(xmit_N_active)';
%   %xmit_apo=ones(1,xmit_N_active);
% 
%   % Do imaging line by line
% 
%   i_start=1;
%   x= -image_width/2 +(i_start-1)*d_x;
% 
%   for i=i_start:no_lines
%   i
%     %   Set the focus for this direction
% 
%     xdc_center_focus (xmit_aperture, [x 0 0]);
%     xdc_focus (xmit_aperture, 0, [x 0 z_focus]);
%     xdc_center_focus (receive_aperture, [x 0 0]);
%     xdc_focus (receive_aperture, focus_times, [x*ones(Nf,1), zeros(Nf,1), focal_zones]);
% 
%     %  Calculate the apodization 
%    
%     xmit_N_pre  = round(x/(width+kerf) + N_elements/2 - xmit_N_active/2);
%     xmit_N_post = N_elements - xmit_N_pre - xmit_N_active;
%     xmit_apo_vector=[zeros(1,xmit_N_pre) xmit_apo zeros(1,xmit_N_post)];
% 
%     rec_N_pre(i) = round(x/(width+kerf) + N_elements/2 - rec_N_active/2);
%     rec_N_post(i) = N_elements - rec_N_pre(i) - rec_N_active;
%  
%     rec_apo_matrix=[zeros(size(focus_times,1),rec_N_pre(i)) rec_apo_matrix_sub zeros(size(focus_times,1),rec_N_post(i))];
% 
%     xdc_apodization (xmit_aperture, 0, xmit_apo_vector);
%     xdc_apodization (receive_aperture, focus_times , rec_apo_matrix);
%    
%     %   Calculate the received response
% 
%     [rf_data, tstart]=calc_scat(xmit_aperture, receive_aperture, positions, amp);
% 
%     %  Store the result
% 
%     cmd=['save sim_flow/rft',num2str(k),'l',num2str(i),'.mat rf_data tstart']
%     eval(cmd)
% 
%     %  Steer in another direction
% 
%     x = x + d_x;
%     
%     end  %  Loop for lines
% 
%   end  %  CFM loop
%   
% %   Free space for apertures
% 
% xdc_free (xmit_aperture)
% xdc_free (receive_aperture)


%% cfm_image.m

%  Read data to make a CFM image from simulated data
%
%  Version 1.0, March 24, 1996, Joergen Arendt Jensen

%  Physical data

f0=3.5e6;                %  Transducer center frequency [Hz]
fs=100e6;                %  Sampling frequency [Hz]
c=1540;                  %  Speed of sound [m/s]
Ncfm=10;                 %  Number of pulses in one direction
fprf=10e3;               %  Pulse emissions frequency  [Hz]
D=5;                     %  Sampling frequency decimation rate
est_dist=1/1000;         %  Distance between velocity estimates [m]
no_lines=20;             %  Number of imaging directions
dy=2/1000;               %  Distance between imaging lines [m]
v_est=0;                 %  Set estimates to zero [m/s]

%  Load the data for one image line and a number of
%  pulse emissions

for i=1:no_lines
  min_sample=0;
  data=0;

  for k=1:Ncfm
    cmd=['load sim_flow/rft',num2str(k),'l',num2str(i),'.mat']
    eval(cmd)
  
    %  Decimate the data and store it in data
  
    if (tstart>0)
      rf_sig = [zeros(fix(tstart*fs-min_sample),1); rf_data];
    else
      rf_sig = rf_data( abs(tstart*fs):max(size(rf_data)) ) ;
      end
     
    rf_sig=hilbert(rf_sig(1:D:max(size(rf_sig))));
    data(1:max(size(rf_sig)),k)=rf_sig;
    end

  %  Make the velocity estimation

  Ndist=floor(2*est_dist/c*fs/D);  %  Rf samples between velocity estimates
  [Nsamples,M]=size(data);
  index=1;

  RMS=std(data(:,1));
  for k=1:Ndist:Nsamples

    %  Find the proper data
  
    vdata=diff(data(k,:));

    %  Calculate the autocorrelation and the velocity

    if (std(vdata) > RMS/5)
      auto  = vdata(2:(M-1)) * vdata(1:(M-2))' ;
      v_est(index,i) = c*fprf/(4*pi*f0) * atan2(imag(auto),real(auto));
    else
      v_est(index,i)=0;
      end
    index=index+1;
    end
  end

[Nx,Ny]=size(v_est);  
imagesc(((0:Ny-1)-Ny/2)*dy*1000,(1:Nx)*Ndist*D/fs*c/2*1000,v_est)
map=[1:64; zeros(2,64)]/64;
colormap(map')
colorbar
ylabel('Depth in tissue [mm]')
xlabel('Lateral distance [mm]')
drawnow

%  Make an interpolated image

ID=25;
[n,m]=size(v_est)
new_est1=zeros(n,m*ID);
for i=1:n
i
  new_est1(i,:)=abs(interp(v_est(i,:),ID));
  end
[n,m]=size(new_est1)
new_est=zeros(n*5,m);
Ndist=Ndist/5;
for i=1:m
i
  new_est(:,i)=abs(interp(new_est1(:,i),5));
  end

[Nx,Ny]=size(new_est);  
new_est=new_est/max(max(new_est))*64;
imagesc(((0:Ny-1)-Ny/2)*dy/ID*1000,(1:Nx)*Ndist*D/fs*c/2*1000,new_est)
map=[1:64; zeros(2,64)]/64;
colormap(map')
colorbar
ylabel('Depth in tissue [mm]')
xlabel('Lateral distance [mm]')
axis('image')



%% cfm_bmode.m


%  Make a combined CFM and B-mode image from
%  simulated data
%
%  Version 1.0, 3/4-97, JAJ

%  Do the gray scale data first

make_image

%  Then do the cfm image

cfm_image

%  Combine the two images
%  Ndist is the number of samples between two velocity estimates

Ndist_new=D*Ndist*fn_bmode/fs;  %  New conversion factor between sampling frequencies

[Nsamples,Nlines]=size(new_env)
[Nest,Nlest]=size(new_est);

for line=1:Nlines
line
  for sample=1:Nsamples
    if (sample/Ndist_new+1 < Nest) 
      if (new_est(fix(sample/Ndist_new+1),line) > 0.25*64)
        comb_image(sample,line)=new_est(fix(sample/Ndist_new+1),line) + 64;
      else
        comb_image(sample,line)=new_env(sample,line)/2;
        end
    else
      comb_image(sample,line)=new_env(sample,line)/2;
      end
    end
  end
  
%  Make a color map for the combined image

map=[[0:63 0:63]; [0:63 zeros(1,64)]; [0:63 zeros(1,64)]]/63;
colormap(map')

%  Display the combined result

clf
dx=40/500;
image( ((1:Nlines)-Nlines/2)*dx,((1:Nsamples)/fn_bmode)*c/2*1000,comb_image)
xlabel('Lateral distance [mm]')
ylabel('Axial distance [mm]')
colormap(map')
%brighten(-0.35)
axis([-20 20 30 90])
axis('image')

print -depsc cfm_image.eps


%% tissue_pht.m

%  Create a computer model of the tissue surrounding a
%  vessel
%
%  Calling: [positions, amp] = tissue_pht (N);
%
%  Parameters:  N - Number of scatterers in the phantom
%
%  Output:      positions  - Positions of the scatterers.
%               amp        - amplitude of the scatterers.
%
%  Version 1.0, March 26, 1997 by Joergen Arendt Jensen

N = 100000;

%  Make the scatteres for a simulation and store
%  it in a file for later simulation use

%   Joergen Arendt Jensen, April 2, 1998

[phantom_positions, phantom_amplitudes] = tissue_phtm (100000);
save pht_data.mat phantom_positions phantom_amplitudes

function [positions, amp] = tissue_phtm (N)

%  Dimensions of the phantom

x_size = 40/1000;   %  Width of phantom [mm]
y_size = 10/1000;   %  Transverse width of phantom [mm]
z_size = 90/1000;   %  Height of phantom [mm]
z_plus = z_size/2 + 10/1000;  %  Start of phantom surface [mm];

%  Initialize the ranges for the vessel

R=0.005;         %  Radius of blood vessel [m]
y_range=2*R;     %  y range for the scatterers  [m]
z_range=2*R;     %  z range for the scatterers  [m]
z_offset=0.06;   %  Offset of the mid-point of the scatterers [m]

%  Creat the general scatterers

x = (rand(N,1)-0.5)*x_size;
y = (rand(N,1)-0.5)*y_size;
z = (rand(N,1)-0.5)*z_size + z_plus - z_offset;

%  Generate the amplitudes with a Gaussian distribution

amp=randn(N,1);

%  Generate the rotated and offset block of sample

theta=-45/180*pi;
xnew=x*cos(theta)+z*sin(theta);
znew=z*cos(theta)-x*sin(theta);

%  Make the vessel and set the amplitudes to -40 dB below inside

inside = (( y.^2 + (znew).^2) < R^2);
amp = amp .* (1-inside) + amp .* inside/100*0; 

%  Generate the rotated and offset block of sample

theta=45/180*pi;
x=xnew*cos(theta)+znew*sin(theta);
z=znew*cos(theta)-xnew*sin(theta) + z_offset;

%  Return the variables

positions=[x y z];
end

