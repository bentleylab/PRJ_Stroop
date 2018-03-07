"""
Functions for preprocessing data for the omission task.
Modified from 'Exploring a single subject ipynb'.

@author: Christine Tseng

20161129: updates by Yvonne Fonken according to omission.py by Arjun
	  Updated load_data, 
"""


### Setup environment
import os, glob, sys
import numpy as np
import pandas as pd
import mne
import tqdm
import h5io
import ecogtools
import csv


### Load data and create Raw Array object ###

def load_data(f):
    """
    Load raw data, bad channels, and return an MNE Raw Array object.

    Input
    -----
    f : string
        File name of data file, usually .fif. Should contain entire
        path, e.g. '/home/knight/ECOGprediction/.../ST44.fif'

    fbad : string
        File name of CSV file with list of bad channels. Should not be entire
        path and is by default set to 'BadChannels.csv'. This file should be
        in the same folder as f.

    Output
    ------
    raw_array : MNE Raw Array object
        Raw Array object with raw data, list of bad channels, etc.
    """
    sub_dir = os.path.dirname(f)
    # Get the list of bad channels
    bads = pd.read_csv(sub_dir + '/BadChannels.csv')["badCh"].tolist()
    print "Bad channels are: %s" %bads

    # Load data
    raw = mne.io.Raw( f, add_eeg_ref=False, preload=True)
    info = mne.create_info(raw.ch_names, raw.info['sfreq'], ch_types='misc')
    print 'There are %d timepoints.'%raw.n_times	

    # Creating raw array object
    raw_array = mne.io.RawArray(raw._data, info)

    # Bad electrode names
    ch_names = raw_array.info['ch_names']
    bad_idx = np.array(bads) - 1
    raw_array.info['bads'] = [ch_names[i] for i in bad_idx]
    print ch_names[bad_idx[0]]
    print ('Bad chans are ')
    print [ch_names[i] for i in bad_idx]
    #print 'Bad Channels: ',', '.join(raw.info['bads']) # This doesn't return content for some reason??

    return raw_array


### Referencing ###

def separate_by_electrode(raw_array):
    """
    Separate raw_array data by depth electrode.
    Supposedly to track which electrodes are on the same shaft? YF

    Input
    -----
    raw_array : MNE Raw array object
        Raw array w/ preprocessed data.

    Output
    ------
    separated : 1 x (# electrodes) list
        List containing the data in raw_array separated by electrode.
        Each entry is structured as follows:
        ['Name of electrode' (e.g. 'RAM'), [List of channel names on electrode],
        np.array(data from electrode)]
    """
    separated = []
    current_electrode_names = []
    start = 0
    prev_name = raw_array.info['ch_names'][0][0:3]

    for ch in range(raw_array.info['nchan']+1):
        if ch == raw_array.info['nchan']: # last iteration
            current_name = ''
        else:
            current_name = raw_array.info['ch_names'][ch][0:3]

        if (current_name == prev_name):
            current_electrode_names.append(raw_array.info['ch_names'][ch])
        else:
            separated.append([prev_name, current_electrode_names,
               raw_array._data[start:ch, :]])
            if ch != raw_array.info['nchan']:
                start = ch
                prev_name = current_name
                current_electrode_names = [raw_array.info['ch_names'][ch]]

    return separated

def wm_csv_to_dict(wmfile):
    """
    Convert csv file of wm channels into dictionary.

    Input
    -----
    wmfile : csv file
        File listing white matter electrodes. Each line in the csv file
        represents one electrode and should contain:
        Col 1: name of electrode
        Col >=2: names of white matter electrodes, one per column

    Output
    ------
    mydict : dictionary (size = # electrodes)
        Dictionary of electrode name and white matter channels.
        Keys: electrode names. Values: names of white matter
        channels. E.g.: {'RAM': ['RAM3', 'RAM4']}
    """
    mydict = {}
    with open(wmfile, mode='r') as infile:
        reader = csv.reader(infile)
        for row in reader:
            channels = [c for _,c in enumerate(row) if c != '']
            mydict[row[0]] = channels[1:]

    return mydict

def wm_indices(separated, wm_electrodes):
    """
    Returns indices of white matter channels in electrodes.
    separated : list of data/channel names separated by electrode
    wm_electrodes : dictionary of white matter electrodes
    """
    inds = [] # Indices to return

    for electrode in separated:
        electrode_name = electrode[0]
        wm_chans = wm_electrodes[electrode_name] # white matter channels
        inds.append([i for i,v in enumerate(electrode[1]) if v in wm_chans])

    return inds

def reference(raw_array, f_car, data, times, method='bipolar', f=None, direction='out', rm_id = None):
    """
    Reference an array using bipolar ('bipolar') or white matter referencing
    ('white') or common average referencing ('car')

    Input
    -----
    raw_array : MNE raw array object
        Preprocessed array with data to be referenced.

    method : string
        Referencing method to use.
        'bipolar': bipolar (default)
        'white': white matter
        'car': common average (only for grid array)

    f : string
        File name of white matter channels if using white matter referencing.

    direction : string
        If using bipolar, direction of referencing.
        'out': starting from small # channels.
        'in': starting from large # channels.

    Output
    ------
    raw_array : MNE raw array object
        Same as input raw_array but with referenced data.

    """
    # Separate the data by electrode
    if method != 'car':
        separated_array = separate_by_electrode(raw_array)

    # Do referencing on separated channels
    if method == 'bipolar':
        for electrode in separated_array:
            data = electrode[2]

            if direction == 'out':
                ch_subtract = np.r_['0, 2', data[1:, :],
                    np.zeros((1, raw_array._data.shape[1]))]
            elif direction == 'in':
                ch_subtract = np.r_['0, 2',
                    np.zeros((1, raw_array._data.shape[1])), data[:-1, :]]
            else:
                raise ValueError('Need to specify direction')

            electrode[2] = data - ch_subtract
    elif method == 'white':
        wm_dict = wm_csv_to_dict(f)
        inds = wm_indices(separated_array, wm_dict)
        print inds

        for i, electrode in enumerate(separated_array):
            if not inds[i]: # if inds[i] is empty
                print electrode[0]
                continue
            data = electrode[2]
            wm_data = data[inds[i], :]
            ref = data - wm_data.mean(axis=0)
            ref[inds[i], :] = 0

            electrode[2] = ref
    elif method == 'car':
	print(f_car )
        groups = pd.read_csv(f_car)["Grouping"].tolist()
	groups = np.array(groups)
	mask = np.ones(len(groups), dtype = bool) # remove channels that may have been removed using the sd method
	mask[rm_id] = False
	
        car_data = ecogtools.car(raw_array._data.T, grouping=groups[mask]).T

    # Move data back into raw_array
    if method != 'car':
        ref_data = [d[2] for _,d in enumerate(separated_array)]
        ref_data = np.vstack(tuple(ref_data))
        raw_array._data = ref_data
    else:
        raw_array._data = car_data

    return raw_array


### Preprocessing ###

def discard_ends(t_start, t_end, raw_array):
    """
    Discard first and last few seconds of data
    """
    n_secs_total = raw_array.n_times/raw_array.info['sfreq']
    start, stop = raw_array.time_as_index([t_start, n_secs_total-t_end])
    data, times = raw_array[:, start:stop]

    return data, times


def demean(raw_array, data):
    """
    Demean all data using mean of data w/discarded ends
    """
    means = np.mean(data, axis=1)
    means = means[:, np.newaxis]
    print raw_array._data[0,0] - raw_array._data[0,1], np.mean(raw_array._data[0,:])
    raw_array._data = raw_array._data - means
    print raw_array._data[0,0] - raw_array._data[0,1], np.mean(raw_array._data[0,:])

    return raw_array

def remove_bad_channels_by_name(raw_array):
    """
    Remove bad channels using list of bad channel names
    """
    good_chs = [ch for ch in raw_array.info['ch_names'] if ch not in
        raw_array.info['bads']]
    #print('remove_bad_channels_by_name')
    #print(raw_array.info['bads'])
    #print(good_chs)
    raw_array=raw_array.pick_channels(good_chs)

    return raw_array



def remove_variable_channels(raw_array, t_discard):
    """
    Drop channels whose stdev is > 2 * mean stdev across all channels
    """
    t_start, t_end = t_discard
    n_secs_total = raw_array.n_times/raw_array.info['sfreq']
    start, stop = raw_array.time_as_index([t_start, n_secs_total-t_end])

    data, times = raw_array[:, start:stop]
    ch_sds = np.std(data, axis=1)

    mean_ch_sd = np.mean(ch_sds)

    bad_sd_ch_idxs = [idx for idx in range(ch_sds.shape[0]) if ch_sds[idx] > 2*mean_ch_sd]
    bad_sd_ch_names = [raw_array.info['ch_names'][i] for i in bad_sd_ch_idxs]
    print "# Bad channels: ", len(bad_sd_ch_idxs), bad_sd_ch_names

    raw_array = raw_array.drop_channels(bad_sd_ch_names)

    return raw_array, bad_sd_ch_idxs


def preprocess_raw_array(raw_array, d, f_car, bp_freq, t_discard,
ref_method, ref_wm_f, ref_direction, ref_groups):
    # Use **kwargs above instead??

    # Discard ends
    t_start, t_end = t_discard
    data, times = discard_ends(t_start, t_end, raw_array)

    # Demean
    if d:
        raw_array = demean(raw_array, data)

    # Remove bad channels
    raw_array = remove_bad_channels_by_name(raw_array)

    # Filter
    bp_l, bp_h = bp_freq
    raw_array = raw_array.filter(l_freq=bp_l, h_freq=bp_h, picks=range(raw_array._data.shape[0]))
    # Notch filter to remove line noise (60 Hz) and its harmonics at 120, 180, 240
    raw_array = raw_array.notch_filter(np.arange(60,241,60),
        picks=range(raw_array._data.shape[0]), notch_widths=1)

    # Remove channels with SD > 2 * mean channel SD
    raw_array, bad_sd_ch_idxs = remove_variable_channels(raw_array, t_discard)

    # Reference
    if ref_groups == {}:
        raw_array = reference(raw_array, f_car, raw_array._data, times, ref_method,
            ref_wm_f, ref_direction, rm_id = bad_sd_ch_idxs)
    else:
        ref_arrays = []
        for method, chans in ref_groups:
            raw_group = raw_array.pick_channels(chans)
            raw_group = reference(raw_group, f_car, raw_group._data, times, method,
                ref_wm_f, ref_direction)
            ref_arrays.append(raw_group)
        raw_array = ref_arrays[0].append(ref_arrays[1:])

    return raw_array

### Load block ###

def load_block (f, ev_id, d=True, fbad='BadChannels.csv', fcar = None, bp_freq=(1.0, 220.0),
t_discard=(1.0, 1.0), ref_method='bipolar', ref_wm_f=None, ref_direction='out',
ref_groups={}):
    """
    Load the data, return processed raw array.
    """
    
	
    # Specify directory where the data is located (multiple blocks not enabled yet)
    sub_dir = os.path.dirname(f) + '/'
    

    f_bad_ch = sub_dir + fbad
    if fcar is not None:  f_car = sub_dir + fcar
    else: f_car = ''
    print 'Files:\n', f+'\n', f_bad_ch+'\n', f_car+'\n'

    # Load data
    raw_array = load_data(f)
    print '\nLoaded data.'
    print 'bad channels output from load_data'
    print(raw_array.info['bads'])

    # Preprocess data
    raw_array = preprocess_raw_array(raw_array, d, f_car, bp_freq, t_discard,
        ref_method, ref_wm_f, ref_direction, ref_groups)
    print '\nPreprocessed data with parameters: \
        \nBand pass frequency: ', bp_freq, \
        '\nDiscard times: ', t_discard, \
        '\nReference method: ', ref_method, \
        '\nWhite matter ref. file: ', ref_wm_f, \
        '\nBipolar ref. direction: ', ref_direction

    return raw_array


### Make epoch ###

def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def replace_omission_with_code(evs, ev_id):
    """
    Label omissions with expected ba (22) or expected ga (23).
    """
    omission_idxs = np.where(evs[:,2]==14)

    n_omissions = len(omission_idxs[0])
    print evs

    for oidx in omission_idxs[0]:
        t_idx = oidx - 3
        if evs[t_idx,2] == 12:
            new_trigger = 23
        elif evs[t_idx,2] == 13:
            new_trigger = 22
        else:
            continue

        evs[oidx, 2] = new_trigger

    trigger_id = evs[:,2]
   # how many 14s are left?
    n_extra_omissions = sum(evs[:,2]==14)

    if sum(evs[:,2] ==19) > 0:
         if 'omission_control' in ev_id: 
	     ev_id['Misc'] = 19
    if sum(evs[:,2] == 16) >0:
         if 'omission_control' not in ev_id:
    	     ev_id['Misc'] = 16
    
	
    if n_extra_omissions > 0:
	 ev_id['omission'] = 14 # make sure to label leftover trials

    print evs

    return evs, ev_id

def remove_bad_epochs(epochs, raw_array, evs, f_bad_segs, t_pre, t_post):
    bt = np.array(pd.read_csv(f_bad_segs))
    
    intervals = [(i[0], i[1]) for i in bt]
    #print intervals

    # event onsets in ms trial boundaries will be (this - t_pre)
    # and (this + t_post)
    ev_onsets = evs[:,0]/epochs.info['sfreq']
    pre_idx = t_pre * epochs.info['sfreq']
    post_idx = t_post * epochs.info['sfreq']

    ev_idxs = raw_array.time_as_index(ev_onsets.tolist())
    bounds = np.vstack([ev_idxs+pre_idx, ev_idxs+post_idx]).T.astype(int)
    #print bounds
    drops = np.zeros(bounds.shape[0], dtype=bool)
    print 'dropping bad epochs here'
    for i, trial in enumerate(bounds):
        ba = trial.tolist()
        #print 'Evaluating trial {i} with index bounds: {ba}...'.format(i=i, ba=ba)
        for interval in intervals:
            t_beg = interval[0]
            t_end = interval[1]
            i_beg, i_end = raw_array.time_as_index([t_beg, t_end])
            ol = getOverlap([i_beg, i_end], ba)

            if ol > 0:
                drops[i] = True
                print 'Trial %d should be dropped because it falls \
                    within %s'%(i, str(interval))
    
    indFakeTrial = evs[:,0] == 3001 #remove fake trials from control
    drops[indFakeTrial] = True
    #indMisc = evs[:,2] == 19
    #drops[indMisc] = True

    epochs.drop(drops, reason='USER: Bad time interval') # changed 'drop_epochs' to 'drop' according to MNE update
    #bounds = bounds[~drops, :]
    #ev_idxs = ev_idxs[~drops]
    print 'evs shape is'
    print evs.shape
    evs = evs[~drops, :]
    print evs.shape
    print 'epochs shape after drops is'
    print epochs

    return epochs, evs, drops

def make_epochs(f_evs, ev_id, raw_array,resample, t, f_bad_seg='BadSegments.csv',
f_output_name='', omission = True):
    """
    Break data into epochs.

    f_evs : string
        Full path to -eve.fif file.

    f_bad_segs : string
        File name containing bad segments

    f_output_name : string
        String of what the exported epochs data should be called

    ev_id : dict
        Dictionary mapping event numbers to what they're called 
	specify in batch file

    t : tuple
        [Pre, Post] times after event to include in each epoch
    """
    sub_dir = os.path.dirname(f_evs) + '/'
    fbadseg = sub_dir + f_bad_seg
   
    print(sub_dir)
    
    # output file name
    block_str = f_evs.split('/')[-2]
    ofn = 'epochData%s-%s'%(f_output_name,block_str)
    print(ofn)

    print 'reading events'
    print f_evs, ev_id
    
    evs = mne.read_events(f_evs)
    

    # evs has a 3-col shape. first is onset, second is trigger_id(=event type),
    # third isall-0s
    assert(all(evs[:,2]==0))
    # we need the trigger_ids in the last column (swap 2 and 3)
    evs[:,2] = evs[:,1]
    evs[:,1] = 0
    assert(all(evs[:,1]==0) & all(evs[:,2] != 0))

    
    # Change omission names
    if omission:
        evs, ev_id = replace_omission_with_code(evs, ev_id)
    

    # use the provided onsets and event types to slice the data into epochs (=trials)
    # this constructor has many options! detrending can be applied here, for example
    # things to consider: detrending, eeg ref
    epochs = mne.Epochs(raw_array, evs, ev_id, *t)
    print 'size of data in epochs immediately after slicing is '
    print epochs

    epochs, evs, drops = remove_bad_epochs(epochs, raw_array, evs, fbadseg, *t)
    print 'size of data in epochs after remove bad epochs'
    print epochs
    print evs.shape 

    # Resample
    if resample:
        epochs.decimate(resample)

    print 'now saving to...'
    print f_evs.split('/')[-4]
    if f_evs.split('/')[-4] == 'Data':
        f2 = os.path.split(sub_dir)
	print(os.path.split(f2[0])[0])
        sub_dir2 = os.path.split(f2[0])[0]
        print sub_dir2
    else:
        sub_dir2 = sub_dir
    fout = sub_dir2 + '/sge/' + ofn
    os.chdir(sub_dir2 + '/sge')
    print 'epochs right before saving'
    print epochs
    print epochs.info
    # save file - ADD THIS LATER
    np.savez(fout, data=epochs.get_data(), info=[epochs.info],
        events=evs, event_id=ev_id, tmin=[epochs.tmin], drops = drops)
    print 'data saved'

    return epochs

def hg_filter(raw_array):
    # bandpass filter for HG range
    n_good_chs = raw_array._data.shape[0]
    raw_array.filter(70, 150, picks=range(n_good_chs))
    # Apply hilbert transform
    print 'Applying hilbert transform...'
    next_pow2 = int(np.ceil(np.log2(raw_array.n_times))) 
    print raw_array, raw_array.info, np.sum(raw_array._data < 0)
    raw_array.apply_hilbert(picks=range(n_good_chs), envelope=True,n_fft=2**next_pow2)
    print raw_array, raw_array.info, np.sum(raw_array._data < 0)
    raw_array.filter(None, 25, picks=range(n_good_chs)) #lowpass at 25Hz, to remove noise and prevent aliasing
    print raw_array, raw_array.info, np.sum(raw_array._data < 0) 
    print '<0: %d'%np.sum(raw_array._data<0)

    return raw_array

