import os
import requests

quality2find = input('Select quality: ')
class2find = input('Select class: ')
comp = input('Select DISP / VEL / ACC: ')

from urllib.parse import urlparse

username = input('Give username: ')
import getpass
password = getpass.getpass('Give password: ')

# import sys, getopt
# arguments = sys.argv
#event2find = str(arguments[1])
# comp = str(arguments[1])

import numpy as np

fname = 'A1_Marsquakes_catalog.txt'
data = np.loadtxt(fname,dtype=str,usecols=(0,1,2,7)).T

eventtime = data[0]
event = data[1]
eventclass = data[2]
quality = data[3]

def data_process_2(event2find):

    # A line printed on the screen will indicate that the code is going to perform the calculations for the specific **event** and **output**, both selected by the user.

    print('Running code for: Event ' + event2find + ' and ' + comp);

    # The code is going to look for the event in the **Marsquake catalog**. This catalog is at the moment a *.txt* file, because this is the only format where the event names, classes and qualities are included. This convention is not included in the *.xml* files, available at the portal, yet. Therefore, the file **A1_Marsquakes_catalog.txt** is necessary to be available at the same directory for this step.

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

    # The code has found the requested event, its **name**, the associated **time** for the start of the waveform (which is before the P arrivals), its **class** and its **quality**, according to the *MQS convention*. Now it prints it on the screen for a check.

    print ('Event: ' + event2find)
    print ('Datetime of the event: ' + time_2)
    print ('Class: ' + class_2)
    print ('Quality: ' + quality_2)

    # The input of the **time** is converted into the UTC convention, using *ObsPy* in order to prepare the data request and downloading. As the data of the IPGP portal are sometimes missing in the start time (reason unknown), it is wise to define a relatively early start time of the requested waveforms. Therefore, 30 minutes before the time of the event and 90 minutes after this time, is considered a good time window to be sure that the event will be contained in every channel and any application of window function (*hamming, Tukey, etc*) will not have an effect to the part of the data concerning the event waveforms.

    from obspy import UTCDateTime
    time = UTCDateTime(time_2)

    starttime = time - 60*30
    print(starttime)
    endtime = time + 60*90
    print(endtime)

    # Now we are ready to proceed to the data download. In order to do this, we should define which kind of data we want. If the user is interested for the VBB InSight data, the channels *BHU, BHV and BHW* should be downloaded. Therefore, we define the appropriate, **network**, **station**, **location** and **channels** for the request.

    net = 'XB'
    sta = 'ELYSE'
    loc = '02'
    chan = 'BH*'

    # The next step is to download the data. For this we will need the authorization for the IPGP portal, a username and a password. The user is asked to provide them. After the authorization the requested data are saved in an *.mseed* file, named **fdsnws_msds.mseed**.

    url = ('https://ws.seis-insight.eu/fdsnws/dataselect/1/query?network=' + net + '&station=' + sta + '&startTime=' + str(starttime) + '&endTime=' + str(endtime) + '&location=*&channel=' + chan)

    r = requests.get(url, auth=(username,password))

    filename=('fdsnws_msds.mseed')

    if r.status_code == 200:
        with open(filename, 'wb') as out:
            for bits in r.iter_content():
                out.write(bits)

    # The *.mseed* file is downloaded and the data should be read and added to a *stream*. The stream parameters are printed in the screen in order to check the downloaded channels.

    import obspy
    from obspy import read

    st = read('fdsnws_msds.mseed')
    print(st)

    # Sometimes more than one traces are available for one component, therefore a merge is necessary.

    for tr in st.select(component='U'):
        st.merge(tr)
    print ('----------------------')
    print('Streams after merging')
    print ('----------------------')
    print(st)

    # On can notice that the channel streams do not have the same starttime and endtime. However, in order to perform the inversion, they should be all of the same length, with equal start and end times. Therefore, a trimming is applied, starting at the latest starttime and ending at the earlier endtime. In the end, the new channel properties are printed on the screen, in order to check that the channel streams have equal starrtime, endtime and number of samples.

    timeE1=st[0].stats.starttime;
    timeN1=st[1].stats.starttime;
    timeZ1=st[2].stats.starttime;

    stime=max(timeE1, timeN1, timeZ1)

    timeEe=st[0].stats.endtime;
    timeNe=st[1].stats.endtime;
    timeZe=st[2].stats.endtime;

    etime=min(timeEe, timeNe, timeZe)

    st[0].trim(stime, etime);
    st[1].trim(stime, etime);
    st[2].trim(stime, etime);

    print ('----------------------')
    print ('Streams after trimming')
    print ('----------------------')
    print(st)

    # The obtained stream is going to make new *.mseed* files, following the naming convention **eventname_channel_COMP.mseed** and save it in a directory which follows the naming convention **Treated/Class/Quality/Eventname/**.

    import os
    path = ('RAW/')
    check = os.path.isdir(path)
    if check == False:
        os.mkdir(path)
    path = ('RAW/' + class_2 + '/')
    check = os.path.isdir(path)
    if check == False:
        os.mkdir(path)
    path = ('RAW/' + class_2 + '/' + quality_2 + '/')
    check = os.path.isdir(path)
    if check == False:
        os.mkdir(path)
    path = ('RAW/' + class_2 + '/' + quality_2 + '/' + event2find + '/')
    check = os.path.isdir(path)
    if check == False:
        os.mkdir(path)
    filename1 = (path + '/' + event2find + '_' + comp + '.mseed')
    st.write(filename1, format = 'MSEED')
    filenameU = (path + '/' + event2find + '_BHU_' + comp + '.sac')
    filenameV = (path + '/' + event2find + '_BHV_' + comp + '.sac')
    filenameW = (path + '/' + event2find + '_BHW_' + comp + '.sac')
    st_new1[0].write(filenameU, format = 'SAC')
    st_new1[1].write(filenameV, format = 'SAC')
    st_new1[2].write(filenameW, format = 'SAC')

    # In the end of the code there will be a warning, but this is not of great importance, as the *.mseed* files are found to be created and saved consistently approprietly.

    print('You may ignore this warning')

    return None

# data_process_2(event2find)
for i in range(0,len(event)):
    class_s = str(eventclass[i])
    if class_s == class2find:
        quality_s = str(quality[i])
        if quality_s == quality2find:
            event_input = str(event[i])
            data_process_2(event_input)
