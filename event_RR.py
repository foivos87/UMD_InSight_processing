#!/usr/bin/env python
# coding: utf-8

# # Download InSight data from IPGP portal ##
# ###### This is a script for the UMD group of InSight, for downloading and treating (rotation and instrument response removal) the InSight data that are available on the dedicated portal.
# *Contributors: Foivos Karakostas, Doyeon Kim, Ross Maguire and the UMD InSight group*
# - - -

# The first input is the **Event name**. The events are cataloged with the number of the Sol (Martian Day), followed by a letter (a, b, c, ...) which indicates if it is the 1st, 2nd, 3rd, etc. event of that day.

event2find = input('Enter event: ')

# Then, the user should select the type of the output. When the instrument response removal is performed, the output can be given in terms of ground displacement, velocity or acceleration.

comp = input('Choose output DISP / VEL / ACC: ')

# A line printed on the screen will indicate that the code is going to perform the calculations for the specific **event** and **output**, both selected by the user.

print('Running code for: Event ' + event2find + ' and ' + comp);

# The code is going to look for the event in the **Marsquake catalog**. This catalog is at the moment a *.txt* file, because this is the only format where the event names, classes and qualities are included. This convention is not included in the *.xml* files, available at the portal, yet. Therefore, the file **A1_Marsquakes_catalog.txt** is necessary to be available at the same directory for this step.

import numpy as np

fname = 'A1_Marsquakes_catalog.txt'
data = np.loadtxt(fname,dtype=str,usecols=(0,1,2,7)).T

eventtime = data[0]
event = data[1]
eventclass = data[2]
quality = data[3]

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

import os
import requests

from urllib.parse import urlparse

username = input('Give username: ')
import getpass
password = getpass.getpass('Give password: ')

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

# When the aformentioned statement is true for the set of channels, the instrument response removal is ready to be performed. In the beginning of the code, the user has already chosen which seismogram wants to obtain, **ground displacement** *(DISP)*, **velocity** *(VEL)* or **acceleration** *(ACC)*. The information for the instrument response is found in the **dataless.XB.ELYSE.seed** file. *(Note: the instrument response is performed by applying a prefiltering in the data, a taper window which is the appropriate for 20 Hz data)*.

from obspy import read_inventory

inv = obspy.read_inventory('dataless.XB.ELYSE.seed')

st_rem1=st.copy()
pre_filt = [0.005, 0.01, 8, 10] #for 20 Hz data
st_rem1.remove_response(output = comp, taper_fraction=0.05, pre_filt = pre_filt, inventory = inv);

# As the instrument response removal is performed, the data should be rotated from the **BHU**, **BHV** and **BHW** channels to **BHE**, **BHN** and **BHZ**. In order to do this, the information available in the *dataless* file, which is already loaded as the *inventory* is necessary. More precisely, the *dip* and *azimuth* of the channels is going to be used in order to perform the rotation of the stream data. The code will print in the screen the channel properties that are taken from the *inventory* and then it will perform the rotation.

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

from obspy.signal.rotate import rotate2zne

(z, n, e) = rotate2zne(st_rem1[0], azs[0], dips[0], st_rem1[1], azs[1], dips[1], st_rem1[2], azs[2], dips[2])

# A Tukey window with 5% of data in the sin function, is applied, in order to remove the artifacts of the high amplitudes at the edges of the waveform.

from scipy import signal

lenz = len(z)
alp = 5e-2
window = signal.tukey(len(z), alpha = alp)
z = z * window
n = n * window
e = e * window

# Now the data for the vertical and horizontal components are available and ready to be added into a new stream, which is the final product of this data processing. The code will print in the end the stats of the new channels with the new names introduced.

st_new1=st_rem1.copy()
st_new1[0].data = z;
st_new1[0].stats.channel = 'BHZ'
st_new1[1].data = n;
st_new1[1].stats.channel = 'BHN'
st_new1[2].data = e;
st_new1[2].stats.channel = 'BHE'
print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
print(st_new1[0].stats)
print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
print(st_new1[1].stats)
print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
print(st_new1[2].stats)

# The obtained stream is going to make new *.mseed* files, following the naming convention **eventname_channel_COMP.mseed** and save it in a directory which follows the naming convention **Treated/Class/Quality/Eventname/**.

import os
path = ('RR/')
check = os.path.isdir(path)
if check == False:
    os.mkdir(path)
path = ('RR/' + class_2 + '/')
check = os.path.isdir(path)
if check == False:
    os.mkdir(path)
path = ('RR/' + class_2 + '/' + quality_2 + '/')
check = os.path.isdir(path)
if check == False:
    os.mkdir(path)
path = ('RR/' + class_2 + '/' + quality_2 + '/' + event2find + '/')
check = os.path.isdir(path)
if check == False:
    os.mkdir(path)
filename1 = (path + '/' + event2find + '_' + comp + '.mseed')
st_new1.write(filename1, format = 'MSEED')
st_new1.write(filename1, format = 'SAC')

# In the end of the code there will be a warning, but this is not of great importance, as the *.mseed* files are found to be created and saved consistently approprietly.

print('You may ignore this warning')
