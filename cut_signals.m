%%% This code is supposed to be used for cut the extra parts of raw
%%% recordings


%% Cut using the trigger

%% Cut using the given time
clc
clear
close all

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
ALLEEG = struct([]);

data_path = "/Users/amin/Documents/Syntropic/HVS/Recordings";
list_of_subjects = {% "SU1_24_RDj",...
                      % "SU1_24_gTo",...
                      % "SU1_24_Q0R",...
                      % "SU1_24_nub",...
                      % "SU1_24_u1W",...
                      "SU1_24_jp6"
                      }
recording_days = {"day_1", "day_5", "day_19"};

fs = 250;
no_light_interval = [1/fs+5, 60*5-5]; % minutes
constant_light_interval = [60*5+1/fs+5, 2*60*5-5]; % minutes
stim_light_interval = [2*60*5+1/fs+5, 30*60-5]; % minutes
stim_no_frames = 30*60*fs;

dataset_no = 0;

for subject = list_of_subjects
    for recording_day = recording_days
        subject = string(subject(1));
        recording_day = string(recording_day(1));
        recording_path = fullfile(data_path, subject, strcat(recording_day, ".txt"));
        data = read_openbci_txt_files(recording_path);
        data_mat = table2array(data(:, 2:9))';
        data_mat_cut = data_mat(:, end-stim_no_frames+1:end);
        recordings_struct.(sprintf(subject)).(sprintf(recording_day)).raw_data = data_mat;
        recordings_struct.(sprintf(subject)).(sprintf(recording_day)).raw_data_cut = data_mat_cut;

        EEG = pop_importdata('dataformat','array','nbchan',0,'data','data_mat_cut','srate',fs,'pnts',0,'xmin',0);
        recording_id = char(strcat(subject, "_", recording_day));
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, dataset_no,'setname',recording_id,'gui','off');
        recordings_struct.(sprintf(recording_id)) = ALLEEG(CURRENTSET).data;

        EEG=pop_chanedit(EEG, {'lookup','sample_locs/Standard-10-20-Cap8.locs'},'load',{'sample_locs/Standard-10-20-Cap8.locs','filetype','autodetect'});
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

        dataset_no = dataset_no + 1;

        EEG = pop_eegfiltnew(EEG, 'locutoff',47,'hicutoff',53,'revfilt',1,'plotfreqz',0);
        recording_id = char(strcat(subject, "_", recording_day, "_NotchFilt"));
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, dataset_no,'setname',recording_id,'gui','off');
        recordings_struct.(sprintf(subject)).(sprintf(recording_day)).data_notch = EEG.data;

        dataset_no = dataset_no + 1;

        EEG = pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',80,'plotfreqz',0);
        recording_id = char(strcat(subject, "_", recording_day, "_NotchFilt", "_BPfilt"));
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, dataset_no,'setname',recording_id,'gui','off');
        recordings_struct.(sprintf(subject)).(sprintf(recording_day)).data_notch_BPF = EEG.data;

        dataset_no = dataset_no + 1;

        EEG = pop_reref( EEG, []);
        recording_id = char(strcat(subject, "_", recording_day, "_NotchFilt", "_BPfilt", "_AvgRR"));
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, dataset_no,'setname',recording_id,'gui','off');
        recordings_struct.(sprintf(subject)).(sprintf(recording_day)).data_notch_BPF_AvgRR = EEG.data;

        dataset_no = dataset_no + 1;

        EEG_full = EEG;

        EEG = pop_select(EEG_full, 'time', no_light_interval);
        recording_id = char(strcat(subject, "_", recording_day, "_NotchFilt", "_BPfilt", "_AvgRR", "_no_ligth"));
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, dataset_no,'setname',recording_id,'gui','off');
        recordings_struct.(sprintf(subject)).(sprintf(recording_day)).data_notch_BPF_AvgRR_noLight = EEG.data;
        dataset_no = dataset_no + 1; 
    
        EEG = pop_select(EEG_full, 'time', constant_light_interval);
        recording_id = char(strcat(subject, "_", recording_day, "_NotchFilt", "_BPfilt", "_AvgRR", "_constant_ligth"));
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, dataset_no,'setname',recording_id,'gui','off');
        recordings_struct.(sprintf(subject)).(sprintf(recording_day)).data_notch_BPF_AvgRR_constantLight = EEG.data;
        dataset_no = dataset_no + 1;

        EEG = pop_select(EEG_full, 'time', stim_light_interval);
        recording_id = char(strcat(subject, "_", recording_day, "_NotchFilt", "_BPfilt", "_AvgRR", "_stimulus_ligth"));
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, dataset_no,'setname',recording_id,'gui','off');
        recordings_struct.(sprintf(subject)).(sprintf(recording_day)).data_notch_BPF_AvgRR_stimulusLight = EEG.data;
        dataset_no = dataset_no + 1;
    end
end
%% Load noise rejected EEG signals from files to run ICA
clc
clear
close all

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
ALLEEG = struct([]);
data_path = "/Users/amin/Documents/Syntropic/HVS/Recordings";
list_of_subjects = {"SU1_24_RDj",...
                      "SU1_24_gTo",...
                      "SU1_24_Q0R",...
                      "SU1_24_nub",...
                      "SU1_24_u1W",...
                      "SU1_24_jp6"
                      }
recording_days = {"day_1", "day_5", "day_19"};

stim_types = {"no_light", ...
            "constant_light", ...
            "stimulus_light"};

fs = 250;

dataset_no = 1;

for subject = list_of_subjects
    for recording_day = recording_days
        for stim_type = stim_types
            subject = string(subject(1));
            recording_day = string(recording_day(1));
            stim_type = string(stim_type(1));

            filepath = char(fullfile(data_path, subject));
            filename = char(strcat(recording_day, "_", stim_type, "_rejected", ".set"));

            EEG = pop_loadset('filename',filename,'filepath',filepath);
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, dataset_no);
            EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'rndreset','yes','interrupt','on','pca',8);
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            pop_selectcomps(EEG, [1:8] );
            EEG = pop_iclabel(EEG, 'default');
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            
            dataset_no = dataset_no + 1;
        end
    end
end