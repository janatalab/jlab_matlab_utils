function F0 = midi2f0(mn)
% F0 = midi2f0(midi_notes)
%
% Converts a vector of MIDI note numbers to corresponding fundamental
% frequencies (F0s)
%
% Assumes that middle-C is midi note 60, and that A440 tuning is used.
%

% 09/13/01 PJ
% 11/03/04 PJ Corrected to have MIDI note 69 = A440

a4 = 69;

F0 = 440 * (2 .^ ((mn-a4)/12));

return
