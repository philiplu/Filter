'''
Filter script to separate non-events from variables (usually periodic), candidate microlensing (one magnification), 
and binary+ microlensing events (one main magnification with second peak during event)
'''

import numpy as np

import matplotlib.pyplot as plt

import os, sys

#Change depending on amount of variation during non-variable phase
mag_tol = .002
#Determine how many points are smoothed together
smoothlength = 10
outfilename = 'ulwdc1_1_293_filter.txt'
outfile = open(outfilename,'w+')
first = ['Number','Type','Pspl','Binary','Burst','Pos Binary','Square']
line_new = ''
for word in first:
	line_new += ('{:<15}'.format(word))
line_new += '\r\n'
outfile.write(line_new)
pspls = []
binaries = []
bursts = []
pbinaries = []
squares = []
periodics = []
flats = []
nobases =[]
longpspl = []

for filenumber in range(1,294):
	print filenumber
	if(filenumber > 99):
		infilenumber = str(filenumber)
	elif(filenumber > 9):
		infilenumber = '0'+str(filenumber)
	else:
		infilenumber = '00'+str(filenumber) 
	infile = '../data challenge/lc/ulwdc1_'+infilenumber+'_W149.txt'

	data = np.loadtxt(infile)
	HJD = data[:,0]
	mag = data[:,1]
	err = data[:,2]

	#Make dates easier to read by subtract 245000
	HJD = [date-2450000 for date in HJD]

	#Initiate event tabulation
	events = [0,0,0,0,0,0]
	pspl = 0
	binary = 1
	burst = 2
	pbinary = 3
	square = 4
	longpspl = 5

	length = len(HJD)

	#Average through
	smoothed = np.zeros(length-smoothlength)
	for i in range(0,length-smoothlength):
		for j in range(0,smoothlength):
			smoothed[i] += mag[i+j]
		smoothed[i] = smoothed[i]/smoothlength

	#find background magnitude
	percents = [np.percentile(smoothed,i) for i in range(0,101)]
	counts = np.zeros(100, dtype=int)
	span = np.zeros(100)
	for i in range(1,101):
		for j in range(0,i):
			if(span[j]<mag_tol/1.5):
				span[j] += percents[i]-percents[i-1]
				if(span[j]>mag_tol/1.5):
					counts[j] = i-j
	#In case many of the points lie within mag_tol
	for j in range(0,100):
		if(counts[j]==0):
			counts[j] = 100-j

	startind = np.argmax(counts)
	endind = startind + counts[startind]
	if(counts[startind]<5):
		for i in range(1,101):
			for j in range(0,i):
				if(span[j]<mag_tol*2):
					span[j] += percents[i]-percents[i-1]
					if(span[j]>mag_tol*2):
						counts[j] = i-j
		#In case many of the points lie within mag_tol
		for j in range(0,100):
			if(counts[j]==0):
				counts[j] = 100-j

		startind = np.argmax(counts)
		endind = startind + counts[startind]
		if(counts[startind]<5):
			data = [str(filenumber),'No Baseline',str(events[pspl]),str(events[binary]),str(events[burst]),str(events[pbinary]),str(events[square])]
			nobases.append(str(filenumber))
			data_line = ''
			for value in data:
				data_line += ('{:<15}'.format(value))
			data_line += '\r\n'
			outfile.write(data_line)
			quit()
	#Take average of lower index to upper index
	baseline = (percents[startind] + percents[endind])/2

	#To calculate standard deviation of background
	background = []
	for i in range(0,length-smoothlength):
		if(baseline-smoothed[i]<0):
			background.append(smoothed[i])
	stdev = np.std(background)
	#Correction factor for noisy lightcurves
	corr = max(1,3*stdev/mag_tol)
	#used to identify burst events which would either have a fast incline or fast dropoff
	edge_tol = .1
	peak_tol = 3*mag_tol*corr

	event_tol = mag_tol * corr
	event_tol2 = event_tol * .5

	area_thresh = .05 * max(1,corr/2)

	#Useful when calculating cuts later on
	eventline = baseline - event_tol

	#Find events
	startev = []
	endev = []
	#few strikes
	endnumber = 0
	for i in range(0,length-smoothlength):
		if(baseline-smoothed[i]>event_tol):
			endnumber = 0
			if(len(startev)==len(endev)):
				startev.append(i)
		elif(baseline-smoothed[i]<event_tol2):
			if(len(startev)>len(endev)):
				if(endnumber<5):
					endnumber += 1
				else:
					endev.append(i-5)
					endnumber = 0
		else:
			endnumber = 0
	if(len(startev)>len(endev)):
		if(startev[-1] != length-smoothlength-1):
			endev.append(length-smoothlength-1)
		else:
			del startev[-1]

	#Find max, and set up for trickle down code
	nevents = len(startev)
	start = 0
	end = 0
	for i in range(0,nevents):
		ncuts = 100
		start = startev[i]
		end = endev[i]
		evmin = np.min(smoothed[start:end])
		evlength = end-start
	#	Start cutting from the top
		if(evmin-eventline<-.25):
			ncuts = min(int(ncuts*(evmin-eventline)/-.25),1000)
		percentcut = 100./ncuts
		peaks = []
		area = 0
		noburst = False
		nobinary = False
	#	Take into account events that are widely separated (i.e. 002)
		if(np.max([HJD[ind+1]-HJD[ind] for ind in range(start,end)])>50):
			noburst = True
		elif(start==0 or end==length-smoothlength-1):
			noburst = True
		elif(end<length-2*smoothlength-1):
			if(HJD[end+smoothlength]-HJD[end]>1):
				noburst = True
				nobinary = True
		for j in range(1,ncuts+1):
			base_tol = (event_tol/3)
			cutbase = eventline + (evmin-eventline)*0.01*(100-j*percentcut)
			cutbase2 = eventline + (evmin-eventline)*0.01*(100-j*percentcut)+base_tol
	#		Repeat similar process as before
			startcut = []
			endcut = []
			endnumber = 0
			endgoal = 1
			for k in range(start,end):
				if(len(startcut)==len(endcut)):
					endnumber = 0
					if(cutbase-smoothed[k]>0):
						startcut.append(k)
				elif(len(startcut)>len(endcut)):
					cuttime = HJD[k]-HJD[startcut[-1]]
					endgoal = min(1+int(np.log(k-startcut[-1])),5)
					if(cuttime > 2):
	#					So that long runs can be more easily kept
						cutbase2 = eventline + (evmin-eventline)*0.01*(100-j*percentcut)+base_tol*max(cuttime/2,3)
					if(cutbase2-smoothed[k]<0):
						if(endnumber<endgoal):
							endnumber += 1
						else:
							endnumber = 0
							endcut.append(k-endgoal)
							cutbase2 = eventline + (evmin-eventline)*0.01*(100-j*percentcut)+base_tol
							endgoal = 1
				else:
					endnumber = 0

			if(len(startcut)>len(endcut)):
	#			So will not have a zero length cut
				if(startcut[-1] != end-1):
					endcut.append(end-1)
				else:
	#				Get rid of extraneous cut
					del startcut[-1]
	#		Find area of each event in the cut
			for l in range(0,len(startcut)):
				cutstart = startcut[l]
				cutend = endcut[l]
				if(cutstart == cutend):
					print 'Error: cutend=cutstart, HJD = ', HJD[cutend]
				peaked = False
				area = 0
				for peak in peaks:
					if(peak>=cutstart and peak<=cutend):
						peaked = True
				if(not peaked):
					logtime = 4*np.log10(min(HJD[cutend]-HJD[cutstart],20))
	#				To create a peak at 4 days or so, to better identify secondary peaks
					seccorr = max(4-.25*((HJD[cutend]-HJD[cutstart]-4)**2),8-4*((HJD[cutend]-HJD[cutstart]-.8)**2))
					timecorr = max(1,seccorr,logtime)
					for m in range(cutstart,cutend):
	#					To prevent very large separation (i.e. between epochs) to have large weight
						timestep = min(HJD[m+1]-HJD[m],.05)
						area += max(np.expm1(60*(cutbase-smoothed[m])),0)*timecorr*timestep/60
					if(np.max(cutbase-smoothed[cutstart:cutend])>peak_tol):
						peaks.append(np.argmin(smoothed[cutstart:cutend])+cutstart)
						noburst = True
					elif(area>area_thresh):
						peaks.append(np.argmin(smoothed[cutstart:cutend])+cutstart)
	#	Interpret the event based on the peaks
		npeak = len(peaks)

		if(npeak>1):
			events[binary] += 1
		elif(npeak==1):
			peakind = peaks[0]
			starttime = HJD[start]
			endtime = HJD[end]
			leftedge = int(start+(end-start)*edge_tol)
			rightedge = int(end-(end-start)*edge_tol)
			leftmin = 0
			rightmin = 0
			nineeight = np.percentile(smoothed[start:end],5)
			if(leftedge!=start and rightedge!=end):
				leftmin=np.min(smoothed[start:leftedge])
				rightmin=np.min(smoothed[rightedge:end])
			if((leftmin<nineeight or rightmin<nineeight) and not noburst):
				if(leftmin<nineeight and rightmin<nineeight):
					events[square] += 1
				else:
					events[burst] += 1
			else:
				posbinary = 0
				leftcut = smoothed[start:peakind]
				rightcut = smoothed[peakind:end]
				lefttimecut = 0
				righttimecut = end-peakind
				for m in range(peakind-start-1):
					timedifference = HJD[start+m+1]-HJD[start+m]
					if(timedifference>.1):
						lefttimecut = m+1
				for m in range(end-peakind-1):
					timedifference = HJD[end-m-1]-HJD[end-m-2]
					if(timedifference>.1):
						righttimecut = end-m-1-peakind

				base = np.percentile(smoothed[start+lefttimecut:peakind+righttimecut],80)
				leftcut = smoothed[start+lefttimecut:peakind]
				rightcut = smoothed[peakind:peakind+righttimecut]
				leftsub = abs(leftcut-base)
				rightsub = abs(rightcut-base)
				peakheight = np.percentile(smoothed[start+lefttimecut:peakind+righttimecut],80)-np.percentile(smoothed[start+lefttimecut:peakind+righttimecut],0)
				if(len(leftsub)<50 or len(rightsub)<50 or (base-np.percentile(leftcut,60))/peakheight>.8 or (base-np.percentile(rightcut,60))/peakheight>.8):
					nobinary = True
				elif(np.min(leftsub)>stdev/2 or np.min(rightsub)>stdev/2 or (eventline-base)/peakheight<.08):
					base = np.percentile(smoothed[start+lefttimecut:peakind+righttimecut],60)
					leftsub = abs(leftcut-base)
					rightsub = abs(rightcut-base)
					peakheight = np.percentile(smoothed[start+lefttimecut:peakind+righttimecut],60)-np.percentile(smoothed[start+lefttimecut:peakind+righttimecut],0)
					if(np.min(leftsub)>stdev/2 or np.min(rightsub)>stdev/2 or (eventline-base)/peakheight<.08):
						base = np.percentile(smoothed[start+lefttimecut:peakind+righttimecut],40)
						leftsub = abs(leftcut-base)
						rightsub = abs(rightcut-base)
						peakheight = np.percentile(smoothed[start+lefttimecut:peakind+righttimecut],40)-np.percentile(smoothed[start+lefttimecut:peakind+righttimecut],0)
						if(np.min(leftsub)>stdev/2 or np.min(rightsub)>stdev/2 or (eventline-base)/peakheight<.08):
							nobinary = True
				if(not nobinary and len(leftcut)>50 and len(rightcut)>50):
					leftstart = np.argmin(leftsub)
					rightend = np.argmin(rightsub)
					leftcut = smoothed[start+lefttimecut+leftstart:peakind]
					rightcut = smoothed[peakind:peakind+rightend+1]
					maxbinary = 0

					for m in range(21):
						leftdifference = (np.percentile(leftcut,m*5)-np.percentile(leftcut,max(m-1,0)*5))/peakheight
						rightdifference = (np.percentile(rightcut,m*5)-np.percentile(rightcut,max(m-1,0)*5))/peakheight
						if(abs((np.percentile(leftcut,m*5)-np.percentile(rightcut,m*5))/max(peakheight,50*stdev))>.05 and (leftdifference<.1 or rightdifference<.1)):
							posbinary += 1
							maxbinary = posbinary
					if(maxbinary>1):
						events[pbinary] += 1
					else:
						events[pspl] += 1
				else:
					events[pspl] += 1

	#interpret lightcurve based on events
	if(events[pspl]+events[binary]+events[burst]+events[pbinary]+events[square]+events[longpspl]>1):
		data = [str(filenumber),'Periodic',str(events[pspl]),str(events[binary]),str(events[burst]),str(events[pbinary]),str(events[square])]
		periodics.append(str(filenumber))
	elif(events[pbinary]>0):
		data = [str(filenumber),'Pos Binary',str(events[pspl]),str(events[binary]),str(events[burst]),str(events[pbinary]),str(events[square])]
		pbinaries.append(str(filenumber))
	elif(events[pspl]==1):
		data = [str(filenumber),'PSPL',str(events[pspl]),str(events[binary]),str(events[burst]),str(events[pbinary]),str(events[square])]
		pspls.append(str(filenumber))
	elif(events[binary]==1):
		data = [str(filenumber),'Binary',str(events[pspl]),str(events[binary]),str(events[burst]),str(events[pbinary]),str(events[square])]
		binaries.append(str(filenumber))
	elif(events[burst]==1):
		data = [str(filenumber),'Burst',str(events[pspl]),str(events[binary]),str(events[burst]),str(events[pbinary]),str(events[square])]
		bursts.append(str(filenumber))
	elif(events[square]==1):
		data = [str(filenumber),'Square',str(events[pspl]),str(events[binary]),str(events[burst]),str(events[pbinary]),str(events[square])]
		squares.append(str(filenumber))
	elif(events[longpspl]==1):
		data = [str(filenumber),'Long PSPL',str(events[pspl]),str(events[binary]),str(events[burst]),str(events[pbinary]),str(events[square])]
		squares.append(str(filenumber))
	else:
		data = [str(filenumber),'Flat',str(events[pspl]),str(events[binary]),str(events[burst]),str(events[pbinary]),str(events[square])]
		flats.append(str(filenumber))
	data_line = ''
	for value in data:
		data_line += ('{:<15}'.format(value))
	data_line += '\r\n'
	outfile.write(data_line)
data = ['PSPL',str(len(pspls)),'Binary',str(len(binaries)),'Burst',str(len(bursts)),'Pbinary',str(len(pbinaries)),'Square',str(len(squares)),'Periodic',str(len(periodics)),'Long PSPL',str(len(longpspl)),'Flat',str(len(flats)),'No Base',str(len(nobases))]
data_line = ''
for value in data:
	data_line += ('{:<10}'.format(value))
data_line += '\r\n'
outfile.write(data_line)

data_line = 'PSPLs\r\n'
outfile.write(data_line)
data_line = ''
counter = 0
for value in pspls:
	data_line += ('{:<4}'.format(value))
	counter += 1
	if(counter == 10):
		data_line += '\r\n'
		counter = 0
data_line += '\r\n'
outfile.write(data_line)

data_line = 'Binaries\r\n'
outfile.write(data_line)
data_line = ''
counter = 0
for value in binaries:
	data_line += ('{:<4}'.format(value))
	counter += 1
	if(counter == 10):
		data_line += '\r\n'
		counter = 0
data_line += '\r\n'
outfile.write(data_line)

data_line = 'Bursts\r\n'
outfile.write(data_line)
data_line = ''
counter = 0
for value in bursts:
	data_line += ('{:<4}'.format(value))
	counter += 1
	if(counter == 10):
		data_line += '\r\n'
		counter = 0
data_line += '\r\n'
outfile.write(data_line)

data_line = 'Possible Binaries\r\n'
outfile.write(data_line)
data_line = ''
counter = 0
for value in pbinaries:
	data_line += ('{:<4}'.format(value))
	counter += 1
	if(counter == 10):
		data_line += '\r\n'
		counter = 0
data_line += '\r\n'
outfile.write(data_line)

data_line = 'Squares\r\n'
outfile.write(data_line)
data_line = ''
counter = 0
for value in squares:
	data_line += ('{:<4}'.format(value))
	counter += 1
	if(counter == 10):
		data_line += '\r\n'
		counter = 0
data_line += '\r\n'
outfile.write(data_line)

data_line = 'Periodics\r\n'
outfile.write(data_line)
data_line = ''
counter = 0
for value in periodics:
	data_line += ('{:<4}'.format(value))
	counter += 1
	if(counter == 10):
		data_line += '\r\n'
		counter = 0
data_line += '\r\n'
outfile.write(data_line)

data_line = 'Long PSPL\r\n'
outfile.write(data_line)
data_line = ''
counter = 0
for value in longpspl:
	data_line += ('{:<4}'.format(value))
	counter += 1
	if(counter == 10):
		data_line += '\r\n'
		counter = 0
data_line += '\r\n'
outfile.write(data_line)

data_line = 'Flats\r\n'
outfile.write(data_line)
data_line = ''
counter = 0
for value in flats:
	data_line += ('{:<4}'.format(value))
	counter += 1
	if(counter == 10):
		data_line += '\r\n'
		counter = 0
data_line += '\r\n'
outfile.write(data_line)

data_line = 'No Bases\r\n'
outfile.write(data_line)
data_line = ''
counter = 0
for value in nobases:
	data_line += ('{:<4}'.format(value))
	counter += 1
	if(counter == 10):
		data_line += '\r\n'
		counter = 0
data_line += '\r\n'
outfile.write(data_line)