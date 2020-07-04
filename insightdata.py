
import sys, getopt

arguments = sys.argv
print('Check control: This should be 3 = ' + str(len(arguments)));

event2find = str(arguments[1])
comp = str(arguments[2])

print('Running code for: Event ' + event2find + ' and ' + comp);

# 1. Run this script in **bash** environment.
# 2. Make sure that the following packages are installed on the **Python3** of the working system:
#     1. NumPy
#     2. Obspy
#     3. Matplotlib
# 3. Make sure that the following files are store in the same working directory with this script:
#     1. A1_Marsquake_catalog.txt
#     2. dataless.XB.ELYSE.seed
#
# - - -

# Import the first package (*numpy*) which is going to be used in order to read the variables from the **Marsquake catalog** file

import numpy as np
import obspy
from obspy import read
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy import read_inventory
from obspy.signal.rotate import rotate2zne
import matplotlib.pyplot as plt

# Ask the user for the event
# event2find = input('Enter Event: ')

def data_process_2(event2find):
    event2find=str(event2find)
    print('Event found: ' + event2find) # check
    print(' ')

    # Import the **Marsquake catalog** file
    fname = 'A1_Marsquakes_catalog.txt'
    data = np.loadtxt(fname,dtype=str,usecols=(0,1,2,7)).T

    eventtime = data[0]
    event = data[1]
    eventclass = data[2]
    quality = data[3]

    # Insert event and find the values

    for i in range(0,len(event)):
        event1 = str(event[i])
        if event1 == event2find:
            a=i

    time_1 = str(eventtime[a])
    t1l=len(time_1)-3
    time_2 = time_1
    class_1 = str(eventclass[a])
    c1l = len(class_1)-3
    class_2 = class_1
    quality_1 = str(quality[a])
    q1l = len(quality_1)-3
    quality_2 = quality_1

    print ('Event: ' + event2find)
    print ('Datetime of the event: ' + time_2)
    print ('Class: ' + class_2)
    print ('Quality: ' + quality_2)

    # Now the process of the data starts, using **ObsPy**
    # In order to download the data from **IRIS**, we should define it as the client.

    client = Client('IRIS')

    time = UTCDateTime(time_2)

    # We take data 30 minutes before until 90 minutes after the event.

    starttime = time - 60*15
    print(starttime)
    endtime = time + 60*45
    print(endtime)

    # We download the data for **BHU**, **BHV** and **BHW** channels.

    net = 'XB'
    sta = 'ELYSE'
    loc = '02'
    chan = 'BH*'

    # We add the data in a variable, by making sure that we include the information about the instrument response.

    st = client.get_waveforms(net, sta, loc, chan, starttime, endtime, attach_response = True)
    print(st)
    # st.plot();

    # We perform the instrument response. In this case we take the *ground velocity* output, the user can change it into the displacement or acceleration, with the appropriate inputs *('DISP' and 'ACC' respectively)*.

    st_rem1=st.copy()
    pre_filt = [0.005, 0.01, 8, 10] #for 20 Hz data
    # comp = input('Choose component (DISP / VEL / ACC): ')
    st_rem1.remove_response(output = comp, taper_fraction=0.05, pre_filt = pre_filt);
    # st_rem1.plot();

    # Apply a filter from 0.1 Hz to 9.9 Hz. This is actually raw processed data, in order to see the events further filtering is needed.

    st_filt1=st_rem1.copy()
    st_filt1.filter('bandpass', freqmin=0.1, freqmax=9.9)
    # st_filt1.plot();

    # Read the properties of the U, V and W component in order to perform the inversion. In order to check, the channel data are printed below this part of the script.

    inv = obspy.read_inventory('dataless.XB.ELYSE.seed')

    sta = inv[0][0]
    azs = []
    dips = []
    trs = []

    channels = ['BHU','BHV','BHW']

    for chn in channels:
        chndata = sta.select(channel=chn)[0]
        print ('CHNDATA--------------------------------')
        print (chndata)
        azs.append(chndata.azimuth)
        dips.append(chndata.dip)

    # Generate the rotated data. A check is provided as a timeseries of the vertical component.

    (z, n, e) = rotate2zne(st_filt1[0], azs[0], dips[0], st_filt1[1], azs[1], dips[1], st_filt1[2], azs[2], dips[2])

    #plt.plot(z, color = 'black', linewidth = 0.04);
    #plt.ylim(-2e-9, 2e-9);

    # Apply a Tukey window (alpha = 5%)

    from scipy import signal

    lenz = len(z)
    alp = 5e-2
    window = signal.tukey(len(z), alpha = alp)
    z = z * window
    n = n * window
    e = e * window

    plt.plot(z, color = 'black', linewidth = 0.04);

    # Create new channel properties for the Z, N and E components.

    st_new1=st_filt1.copy()
    st_new1[0].data = z;
    st_new1[0].stats.channel = 'BHZ'
    st_new1[1].data = n;
    st_new1[1].stats.channel = 'BHN'
    st_new1[2].data = e;
    st_new1[2].stats.channel = 'BHE'
    print(st_new1[0].stats)
    print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
    print(st_new1[1].stats)
    print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
    print(st_new1[2].stats)

    # See the treated data. In most cases, mostly in low quality and not broadband events, only noise will be shown in these plots, as further filtering is needed for events that are observed in narrow frequency ranges.

    dataZ = st_new1[0].data
    timesZ = st_new1[0].times()
    plt.figure(figsize=(20,10))
    plt.subplot(311)
    plt.plot(timesZ, dataZ, color = 'black', linewidth = '0.1');
    plt.ylim(-2e-9, 2e-9)
    plt.xlim(1, 3600);
    plt.ylabel('Ground Velocity (m/s)', fontsize = '14')
    plt.suptitle(event2find, fontsize = '14');
    plt.title('BHZ', fontsize = '14');

    dataN = st_new1[1].data
    timesN = st_new1[1].times()
    plt.subplot(312)
    plt.plot(timesN, dataN, color = 'black', linewidth = '0.1');
    plt.ylim(-2e-9, 2e-9)
    plt.xlim(1, 3600);
    plt.ylabel('Ground Velocity (m/s)', fontsize = '14')
    plt.title('BHN', fontsize = '14');

    dataE = st_new1[2].data
    timesE = st_new1[2].times()
    plt.subplot(313)
    plt.plot(timesE, dataE, color = 'black', linewidth = '0.1');
    plt.ylim(-2e-9, 2e-9)
    plt.xlim(1, 3600);
    plt.ylabel('Ground Velocity (m/s)', fontsize = '14')
    plt.xlabel('Time (s)', fontsize = '14');
    plt.title('BHE', fontsize = '14');

    # Write new *.mseed* files for the treated data, with instrument response removal and rotation. In this part the data are categorized under a directory tree following class and quality.

    import os
    path = ('Treated/')
    check = os.path.isdir(path)
    if check == False:
        os.mkdir(path)
    path = ('Treated/' + class_2 + '/')
    check = os.path.isdir(path)
    if check == False:
        os.mkdir(path)
    path = ('Treated/' + class_2 + '/' + quality_2 + '/')
    check = os.path.isdir(path)
    if check == False:
        os.mkdir(path)
    path = ('Treated/' + class_2 + '/' + quality_2 + '/' + event2find + '/')
    check = os.path.isdir(path)
    if check == False:
        os.mkdir(path)
    filename1 = (path + '/' + event2find + '_BHZ_' + comp + '.mseed')
    st_new1[0].write(filename1, format = 'MSEED')
    filename2 = (path + '/' + event2find + '_BHN_' + comp + '.mseed')
    st_new1[1].write(filename2, format = 'MSEED')
    filename3 = (path + '/' + event2find + '_BHE_' + comp + '.mseed')
    st_new1[2].write(filename3, format = 'MSEED')

    return None

data_process_2(event2find)
