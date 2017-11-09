# CoupledFDPlateAndString
A collection of Finite Difference models coupling a Kirchhoff Thin Plate and Stiff String along with related functions 

NOTE: The included files were created with MATLAB 2016b on macOS. Some behaviour may be slightly different, particularly with
vector .* element-wise multiplication. I've routed out most of these but there may still be a few lurking in bushes.

The files that should be run are:

thin_plate_stiff_string_MIDI.m
thin_plate_stiff_string_loss.m
thin_plate_reverb.m
thin_plate_loss_xy.m
reverbStringPlate.m
multiStringPianoModel.m


There are respective instrument files which control the parameters stated in the first few lines of each script.

These are:

plate_and_string_vars.m
multiStringPianoVars.m
plateVerbVars.m
reverbParams.m
