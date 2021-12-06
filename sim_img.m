


%  Example of use of the new Field II program running under Matlab
%
%  This example shows how a linear array B-mode system scans an image
%
%  This script assumes that the field_init procedure has been called
%  Here the field simulation is performed and the data is stored
%  in rf-files; one for each rf-line done. The data must then
%  subsequently be processed to yield the image. The data for the
%  scatteres are read from the file pht_data.mat, so that the procedure
%  can be started again or run for a number of workstations.
%
%  Version 2.1 by Joergen Arendt Jensen, August 14, 1998 for Matlab 5.

%  Generate the transducer apertures for send and receive

f0=3e6;                  %  Transducer center frequency [Hz]
fs=100e6;                %  Sampling frequency [Hz]
c=1540;                  %  Speed of sound [m/s]
lambda=c/f0;             %  Wavelength [m]
width=lambda;            %  Width of element
element_height=5/1000;   %  Height of element [m]
kerf=0.05/1000;          %  Kerf [m]
focus=[0 0 70]/1000;     %  Fixed focal point [m]
N_elements=196;          %  Number of physical elements
N_active=64;             %  Number of active elements

% Do not use triangles

set_field('use_triangles',0);

%  Set the sampling frequency

set_sampling(fs);

%  Generate aperture for emission

emit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 1,focus);

%  Set the impulse response and excitation of the emit aperture

impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (emit_aperture, impulse_response);

excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (emit_aperture, excitation);

%  Generate aperture for reception

receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 1,focus);

%  Set the impulse response for the receive aperture

xdc_impulse (receive_aperture, impulse_response);

%   Set a Hanning apodization on the apertures

apo=hanning(N_active)';

%   Load the computer phantom data

load pht_data

%  Set the different focal zones for reception

focal_zones=[30:20:200]'/1000;
Nf=max(size(focal_zones));
focus_times=(focal_zones-10/1000)/1540;
z_focus=70/1000;          %  Transmit focus

%   Do linear array imaging

no_lines=50;                            %  Number of lines in image
image_width=40/1000;                    %  Size of image sector
d_x=image_width/no_lines;               %  Increment for image
for i=[1:no_lines]
i
   %  Find the direction for the imaging
  
  x= -image_width/2 + (i-1)*d_x;

   %   Set the focus for this direction

  xdc_center_focus (emit_aperture, [x 0 0]);
  xdc_focus (emit_aperture, 0, [x 0 z_focus]);
  xdc_center_focus (receive_aperture, [x 0 0]);
  xdc_focus (receive_aperture, focus_times, [x*ones(Nf,1), zeros(Nf,1), focal_zones]);

   %  Calculate the apodization 
   
  N_pre  = round(x/(width+kerf) + N_elements/2 - N_active/2);
  N_post = N_elements - N_pre - N_active;
  apo_vector=[zeros(1,N_pre) apo zeros(1,N_post)];
  xdc_apodization (emit_aperture, 0, apo_vector);
  xdc_apodization (receive_aperture, 0, apo_vector);
  
  %   Calculate the received response

  [rf_data, tstart]=calc_scat(emit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);

  %  Store the result

  cmd=['save sim_bmd/rf_ln',num2str(i),'.mat rf_data tstart']
  eval(cmd)

  %  Steer in another angle

  x = x + d_x;
  end

%   Free space for apertures

xdc_free (emit_aperture)
xdc_free (receive_aperture)
