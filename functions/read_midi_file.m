% MIDI Read in Function
%----------------------------------------------%
%----------------------------------------------%
function midi_notes = read_midi_file(filename,BPM)

  % The following resources were invaluable.
  % http://www.ccarh.org/courses/253/handout/smf/
  % http://sites.uci.edu/camp2014/2014/05/19/timing-in-midi-files/
  % https://ccrma.stanford.edu/~craig/articles/linuxmidi/misc/essenmidi.html
  % http://www.indiana.edu/~emusic/etext/MIDI/chapter3_MIDI4.shtml
  % http://computermusicresource.com/MIDI.Commands.html

  % ---------------------------------------------------- %
  % Start Open File
  FID = fopen(filename);
  [DATA count] = fread(FID,'uint8'); %read in 8 bit unsigned
  fclose(FID);
  % End Open File
  % ---------------------------------------------------- %


  % ---------------------------------------------------- %
  % MIDI Info
  midi.track_count = uint82dec(DATA(11:12));
  midi.ticks = uint82dec(DATA(13:14));

  if nargin<2
    midi.tempo = 500000;
  else
    midi.tempo = 60/(BPM*1e-6);
  end
  % tempo in BPM = 60/(500000* 1e-6)
  % midi.clock = 60/(BPM*1e-6)
  % ---------------------------------------------------- %


  % ----------------------------------------------------%
  % Split up MIDI Tracks
  byte_cnt = 15; %skip to 15th byte as the rest is header data

  % track_chunk = "MTrk" + <length> + <track_event> [+ <track_event> ...]
  for i=1:midi.track_count

    byte_cnt = byte_cnt+4;
    track_length = uint82dec(DATA(byte_cnt:byte_cnt+3));
    byte_cnt = byte_cnt+4;

    % track_bytes holds initial 8 bytes + rest of track data
    track_bytes{i} = DATA((byte_cnt-8):(byte_cnt + track_length-1));
    byte_cnt = byte_cnt+track_length;
  end
  % ----------------------------------------------------%



  % ----------------------------------------------------%
  % Iterate through MIDI tracks
  for i=1:midi.track_count

    track = track_bytes{i};

    message_count = 1;

    % Skip the firt 8 Bytes which are all Header
    % Byte_cnt increments the relevant number of bytes
    % given the MIDI data type
    byte_cnt=9;

    % ----------------------------------------------------%
    % Track Data While Loop
    while (byte_cnt < length(track_bytes{i}))
      %  track_event = <v_time> + <midi_event> | <meta_event> | <sysex_event>
      clear current_message;
      [delta_time,byte_cnt] = calculate_v_time(track, byte_cnt);
      % ----------------------------------------------------%
      % Meta-event beginning with 'FF': 0xFF
      % meta_event = 0xFF + <meta_type> + <v_length> + <event_data_bytes>
      if track(byte_cnt)==255
        meta_type = track(byte_cnt+1); %<meta_type> 1 byte
        byte_cnt = byte_cnt+2;
        % get variable length 'length' field
        [v_length,byte_cnt] = calculate_v_time(track, byte_cnt); % <v_length> variable
        event_data_bytes = track(byte_cnt:byte_cnt+v_length-1);
        byte_cnt = byte_cnt + v_length;
        meta_event = 0;
        % ----------------------------------------------------%


        % ----------------------------------------------------%
      else
        meta_event = 1;
        % keep reading 16 bytes
        if (track(byte_cnt)<128)
          current_byte = last_byte;
        else
          current_byte  = track(byte_cnt);
          byte_cnt = byte_cnt + 1;
        end

        if (bitshift(current_byte,-4)>=8 && bitshift(current_byte,-4)<=14) % Channel Voice Message

          meta_type = bitshift(bitshift(current_byte,-4),4);
          v_length = channel_command_size(meta_type); % var length data:
          event_data_bytes = track(byte_cnt:byte_cnt+v_length-1);%the actual event data.
          byte_cnt = byte_cnt + v_length;

        end

        last_byte = bitand(current_byte,15) + bitshift(bitshift(current_byte,-4),4);

      end % end midi event 'if'

      % compile message elements into a structure
      current_message.delta_time = delta_time;
      current_message.meta_event = meta_event;
      current_message.meta_type = meta_type;
      current_message.data = event_data_bytes;

      %Message structure is placed in the approriate message index of the
      % MIDI data structure.
      midi.track(i).messages(message_count) = current_message;
      message_count = message_count + 1;

    end % end while loop
    % ----------------------------------------------------%
  end % end track loop

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % GET MIDI notes
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  track_index = 1:midi.track_count;

  % set tempo MIDI tempo

  notes = zeros(0,3); % Note Number, Velocity, start time

  for i=1:length(track_index)
    tracknum = track_index(i);
    IncrementDeltatime=0; %accumulator for delta time values
    note_time=0; % note_time in seconds

    for MIDI_Msg_number=1:length(midi.track(tracknum).messages)

      current_message = midi.track(tracknum).messages(MIDI_Msg_number);

      % get vairables from current MIDI message
      meta_event  = current_message.meta_event;
      delta_time   = current_message.delta_time;
      data        = current_message.data;
      meta_type   = current_message.meta_type;

      IncrementDeltatime = IncrementDeltatime + delta_time;

      % since the clock is in microseconds, the time between events is
      % the clock_delta mulitplied by tempo and divided by 1 million and
      % the number of ticks
      note_time = note_time + delta_time*1e-6*midi.tempo/midi.ticks;

      % from the data stream, note on messages will begin will
      % have a meta_type == 144 and velocity of more than 1
      % data(1) in this case will represent the note number
      % from http://computermusicresource.com/MIDI.Commands.html
      if (meta_event==1 && meta_type==144 && data(2)>0)
        % note on:
        notes(end+1,:) = [data(1) data(2) note_time];
      end

    end %% end track.

  end

  % Index rows by note start time
  [Y,order] = sort(notes(:,3));
  midi_notes = notes(order,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIDI BASED Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% http://computermusicresource.com/MIDI.Commands.html
function size=channel_command_size(chan_msg)

  switch chan_msg

    case {128, 144, 160, 176, 224}

    size=2;

    % Program Number and pressure only take up 1 byte
    case {192,208}

    size=1;

  end
end

% caluclate variable length
function [value,pointer] = calculate_v_time(bytes, pointer)

  value = 0;
  flag=1;
  while (flag)
    % if MSB=1, then delta-time continues into next byte...
    if(~bitand(bytes(pointer),128))
      flag=0;
    end
    % append last 7 bits from each byte in the delta_time:
    value = value*128 + rem(bytes(pointer), 128);
    pointer=pointer+1;
  end
end

function num=uint82dec(array)
  % change from and 8bit array to decimal
  % need to read bytes in reverse apparently, which was a bugger
  % I'll tell you
  num = (2.^(8*(0:length(array)-1)))*array(end:-1:1);
end
%----------------------------------------------%
