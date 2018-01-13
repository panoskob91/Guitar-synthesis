function [ out,Fs] = Synthesis_6GuitarString( f0_E,f0_A,f0_D,f0_G,f0_B,f0_E2 )
%Synthesis_6GuitarString: This function applies an algorithm for synthesising a six string guitar.  
%       This function has 6 inputs representing the fundamental frequency of
%   each string of the notes played.
%   f0_E : The frequency of the note played in the low E string. 
%   f0_A : The frequency of the note played in the A string.
%   f0_D : The frequency of the note played in the D string.
%   f0_G : The frequency of the note played in the G string.
%   f0_B : The frequency of the note played in the B string.
%   f0_E2 : The frequency of the note played in the high E string.
%   
%   Type (0,0,0,0,0,0) to strum open strings.
%
%       PRESETS: Chords
%   Am : (82.41,110,164.81,220,261.63,329.63) (input arguments)
%   Bm : (92.5,123.47,185,246.94,293.66,369.99)
%   Cm : (98,130.81,196,261.63,311.13,392)
%   Dm : (82.41,110,146.83,220,293.66,349.23)
%   Em : (82.41,123.47,164.81,196,246.94,329.63)
%   Fm : (87.31,130.81,174.61,207.65,261.63,349.23)
%   Gm : (98,146.83,196,233.08,293.66,392)
%
%   A  : (82.41,110,196,220,277.18,329.63)
%   B  : (92.5,123.47,185,246.94,311.13,369.99)
%   C  : (82.41,130.81,164.81,196,261.63,329.63)
%   D  : (82.41,110,146.83,220,293.66,369.99)
%   E  : (82.41,123.47,164.81,207.65,246.94,329.63)
%   F  : (87.31,130.81,174.61,220,261.63,349.23)
%   G  : (98,123.47,146.83,196,246.94,392)
%   
%   A#m: (87.31,116.54,174.61,233.08,277.18,349.23)
%   C#m: (103.83,155.56,207.65,246.94,311.13,415.30)
%   D#m: (82.41,116.54,155.56,233.08,311.13,369.99)
%   F#m: (87.31,138.59,185,220,277.18,369.99)
%   G#m: (103.83,155.56,207.65,246.94,311.13,415.30)
% 
%   A#: (87.31,116.54,174.61,233.08,293.66,349.23)
%   C#: (103.83,138.59,207.65,277.18,349.23,415.30)
%   D#: (98,155.56,196,233.08,311.13,392)
%   F#: (92.5,138.59,185,253.08,277.18,369.99)
%   G#: (103.83,138.59,207.65,261.63,311.13,415.30)
%
%%

%Fs = 176400;         % Sample Rate for audio output
Fs = 44100;
%T = 1/Fs;           % Sample period
N = 10*Fs;          % Number of samples in output
%f0 = 441;         % Fundamental Frequency of string in Hz.


% Code for playing a default sound vector.
if (f0_E==0 && f0_A==0 && f0_D==0 && f0_G==0 && f0_B==0 && f0_E2==0)

%if (nargin <= 6)  %Default values if input arguments are blank. 
f0_E = 82.41;
f0_A = 110;
f0_D = 146.83;
f0_G = 196;
f0_B =246.94;
f0_E2 = 329.63;
else             % User's values.
f0_E;  
f0_A;
f0_D;
f0_G;
f0_B;
f0_E2;
end

%%
% Band Pass Filters simulating the body of a guitar
bpFilterE1 = designfilt('bandpassfir','filterOrder',200,...
    'CutoffFrequency1',223,'CutoffFrequency2',280,...
    'SampleRate',Fs);
bpFilterA = designfilt('bandpassfir','filterOrder',200,...
    'CutoffFrequency1',217,'CutoffFrequency2',260,...
    'SampleRate',Fs);
bpFilterD = designfilt('bandpassfir','filterOrder',200,...
    'CutoffFrequency1',440,'CutoffFrequency2',470,...
    'SampleRate',Fs);
bpFilterG = designfilt('bandpassfir','filterOrder',200,...
    'CutoffFrequency1',146,'CutoffFrequency2',246,...
    'SampleRate',Fs);
bpFilterB = designfilt('bandpassfir','filterOrder',200,...
    'CutoffFrequency1',239,'CutoffFrequency2',259,...
    'SampleRate',Fs);
bpFilterE2 = designfilt('bandpassfir','filterOrder',200,...
    'CutoffFrequency1',279.63,'CutoffFrequency2',379.63,...
    'SampleRate',Fs);
%%

L_E = floor(0.5*Fs/f0_E);   % String Length in samples
L_A = floor(0.5*Fs/f0_A);
L_D = floor(0.5*Fs/f0_D);
L_G = floor(0.5*Fs/f0_G);
L_B = floor(0.5*Fs/f0_B);
L_E2 = floor(0.5*Fs/f0_E2);

x = 0.5;                % Pluck Position as a proportion of string length
pickup_E = floor(L_E/2);    % Output Position
pickup_A = floor(L_A/2);
pickup_D = floor(L_D/2);
pickup_G = floor(L_G/2);
pickup_B = floor(L_B/2);
pickup_E2 = floor(L_E2/2);



r = -0.99;                 % Bridge Reflection Coefficient
%r = -0.8;

% Right-going delay line, defined by L 
right_E = zeros(1,L_E);
right_A = zeros(1,L_A);
right_D = zeros(1,L_D);
right_G = zeros(1,L_G);
right_B = zeros(1,L_B);
right_E2 = zeros(1,L_E2);

% Left-going delay line, defined by L
left_E = zeros(1,L_E);
left_A = zeros(1,L_A);
left_D = zeros(1,L_D);
left_G = zeros(1,L_G);
left_B = zeros(1,L_B);
left_E2 = zeros(1,L_E2);

% Define initial string shape from pluck position x
pluck_E = x*(L_E-1);% Find point on string corresponding to pluck position
pluck_A = x*(L_A-1);
pluck_D = x*(L_D-1);
pluck_G = x*(L_G-1);
pluck_B = x*(L_B-1);
pluck_E2 = x*(L_E2-1);

x_E = [ ((0:floor(pluck_E))/pluck_E),(L_E-1-((floor(pluck_E)+1):(L_E-1)))/(L_E-1-pluck_E)];
x_A = [ ((0:floor(pluck_A))/pluck_A),(L_A-1-((floor(pluck_A)+1):(L_A-1)))/(L_A-1-pluck_A)];
x_D = [ ((0:floor(pluck_D))/pluck_D),(L_D-1-((floor(pluck_D)+1):(L_D-1)))/(L_D-1-pluck_D)];
x_G = [ ((0:floor(pluck_G))/pluck_G),(L_G-1-((floor(pluck_G)+1):(L_G-1)))/(L_G-1-pluck_G)];
x_B = [ ((0:floor(pluck_B))/pluck_B),(L_B-1-((floor(pluck_B)+1):(L_B-1)))/(L_B-1-pluck_B)];
x_E2 = [ ((0:floor(pluck_E2))/pluck_E2),(L_E2-1-((floor(pluck_E2)+1):(L_E2-1)))/(L_E2-1-pluck_E2)];


% Initial displacement for each delay line is equivalent to plucked string
% excitation shape divided by 2.
% Delay line representing left going wave 
left_E(1:L_E) = x_E(1:L_E)/2;
left_A(1:L_A) = x_A(1:L_A)/2;
left_D(1:L_D) = x_D(1:L_D)/2;
left_G(1:L_G) = x_G(1:L_G)/2;
left_B(1:L_B) = x_B(1:L_B)/2;
left_E2(1:L_E2) = x_E2(1:L_E2)/2;

% Delay line representing right going wave 
right_E(1:L_E) = x_E(1:L_E)/2;
right_A(1:L_A) = x_A(1:L_A)/2;
right_D(1:L_D) = x_D(1:L_D)/2;
right_G(1:L_G) = x_G(1:L_G)/2;
right_B(1:L_B) = x_B(1:L_B)/2;
right_E2(1:L_E2) = x_E2(1:L_E2)/2;


% Initialize output
out_E = zeros(1,N);
out_A = zeros(1,N);
out_D = zeros(1,N);
out_G = zeros(1,N);
out_B = zeros(1,N);
out_E2 = zeros(1,N);
% Initialize variables for display
%pkval = max(abs(x));
%string_pos = 1:L;

% Main digital waveguide loop
for n = 1:N
  
  % Shift left-going wave one step left; append dummy value for now
  left_E = [left_E(2:L_E),0];
  left_A = [left_A(2:L_A),0];
  left_D = [left_D(2:L_D),0];
  left_G = [left_G(2:L_G),0];
  left_B = [left_B(2:L_B),0];
  left_E2 = [left_E2(2:L_E2),0];
 
  % New right-going value is negative of new value at nut of left-going
  
  nut_E = -left_E(1);
  nut_A = -left_A(1);
  nut_D = -left_D(1);
  nut_G = -left_G(1);
  nut_B = -left_B(1);
  nut_E2 = -left_E2(1);
  
  % Add reflection from nut into first element of right-going delay line;
  % Shift right-going wave one step
  right_E = [nut_E, right_E(1:L_E-1)];
  right_A = [nut_A, right_A(1:L_A-1)];
  right_D = [nut_D, right_D(1:L_D-1)];
  right_G = [nut_G, right_G(1:L_G-1)];
  right_B = [nut_B, right_B(1:L_B-1)];
  right_E2 = [nut_E2, right_E2(1:L_E2-1)];
  
  % At the 'bridge' (right-hand end), assume perfect reflection (* -1).
  % New left-going value is negative of new value at bridge of right-going
  
  bridge_E = (r)*right_E(L_E);
  bridge_A = (r)*right_A(L_A);
  bridge_D = (r)*right_D(L_D);
  bridge_G = (r)*right_G(L_G);
  bridge_B = (r)*right_B(L_B);
  bridge_E2 = (r)*right_E2(L_E2);
  
  % Add new bridge value to end of left-going delay line, replacing dummy
  % value from above:
  
  left_E(L_E) = bridge_E;
  left_A(L_A) = bridge_A;
  left_D(L_D) = bridge_D;
  left_G(L_G) = bridge_G;
  left_B(L_B) = bridge_B;
  left_E2(L_E2) = bridge_E2;

  % Output is sum of left and right going delay lines at pickup point.
  % Calculate output:
  out_E(n) = (left_E(pickup_E) + right_E(pickup_E));
  out_A(n) = (left_A(pickup_A) + right_A(pickup_A))./2;
  out_D(n) = (left_D(pickup_D) + right_D(pickup_D))./3;
  out_G(n) = (left_G(pickup_G) + right_G(pickup_G))./4;
  out_B(n) = (left_B(pickup_B) + right_B(pickup_B))./5;
  out_E2(n) = (left_E2(pickup_E2) + right_E2(pickup_E2))./6;
  
end

out_E = filter(bpFilterE1,out_E);
out_E = out_E./max(abs(out_E));
out_A = filter(bpFilterA,out_A);
out_A = out_A./max(abs(out_A));
out_D = filter(bpFilterD,out_D);
out_D = out_D./max(abs(out_D));
out_G = filter(bpFilterG,out_G);
out_G = out_G./max(abs(out_G));
out_B = filter(bpFilterB,out_B);
out_B = out_B./max(abs(out_B));
out_E2 = filter(bpFilterE2,out_E2);
out_E2 = out_E2./max(abs(out_E2));
 
out = out_E + out_A + out_D + out_G + out_B + out_E2;
out = out./max(abs(out));

%Plot time domain response of output
figure(2);
t = linspace(0,N/Fs,N);
plot(t,out);
xlabel('Time (s)');
ylabel('Amplitude');

fftSize = 8192;
f = (0:fftSize-1)*(Fs/fftSize);

% Plot frequency domain response of output
figure(3);
semilogx(f,20*log10(abs(fft(out,fftSize))/max(abs(fft(out,fftSize)))));
axis([0 16000 -40 0]);
xlabel('Frequency (Hz)');
ylabel('Magnitude Response (dB)');

sound(out, Fs);


end

