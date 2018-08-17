'''
Filter script to separate non-events from variables (usually periodic), candidate microlensing (one magnification), 
and binary+ microlensing events (one main magnification with second peak during event)
'''

import numpy as np

import matplotlib.pyplot as plt

import os, sys, argparse

#Change depending on amount of variation during non-variable phase
mag_tol = .002
#Determine how many points are smoothed together
smoothlength = 10
#Take inputs
parser = argparse.ArgumentParser(description='Input file')
parser.add_argument('infile',type=str)

args = parser.parse_args()

data = np.loadtxt(args.infile)
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
	span = np.zeros(100)
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
		print 'Cant find baseline'
		print startind
		print counts[startind]
		quit()
#Take average of lower index to upper index
baseline = (percents[startind] + percents[endind])/2
print 'Baseline = ',baseline

print counts[startind]
print startind

#To calculate standard deviation of background
background = []
for i in range(0,length-smoothlength):
	if(baseline-smoothed[i]<0):
		background.append(smoothed[i])
stdev = np.std(background)
print 'stdev = ',stdev
#Correction factor for noisy lightcurves
corr = max(1,3*stdev/mag_tol)
#used to identify burst events which would either have a fast incline or fast dropoff
edge_tol = .1
peak_tol = 3*mag_tol*corr

event_tol = mag_tol * corr
event_tol2 = event_tol * .3

area_thresh = .05 * max(1,corr/2)

#Useful when calculating cuts later on
eventline = baseline - event_tol

print 'eventline = ',eventline

#Find events
startev = []
endev = []
#few strikes
endnumber = 0
for i in range(0,length-smoothlength):
	if(len(startev)==len(endev)):
		if(baseline-smoothed[i]>event_tol):
			startev.append(i)
	elif(len(startev)>len(endev)):
		smoothlength2 = min((i-startev[-1])/10,length-smoothlength-i-1)
		if(smoothlength2>20):
			smoothed2 = 0
			for j in range(smoothlength2):
				smoothed2 += smoothed[i+j]
			smoothed2 = smoothed2/smoothlength2
		else:
			smoothed2 = smoothed[i]
		if(baseline-smoothed2<event_tol2 and baseline-smoothed[i]<event_tol2):
			if(endnumber<5):
				endnumber += 1
			else:
				endev.append(i)
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
	if(evlength>1000):
		for i in range(21):
			print i*5,' percent = ',np.percentile(smoothed[start:end],i*5)
		for i in range(20):
			print (i+1)*5,' difference = ',np.percentile(smoothed[start:end],(i+1)*5) - np.percentile(smoothed[start:end],i*5)
#	Start cutting from the top
	if(evmin-eventline<-.02):
		ncuts = min(int(ncuts*(evmin-eventline)/-.02),3000)
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
	elif(end < length - 2*smoothlength):
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
				endgoal = min(1+int(np.log(k-startcut[-1])),8)
				cutbase2 = eventline + (evmin-eventline)*0.01*(100-j*percentcut)+base_tol
				if(cuttime > 2):
#					So that long runs can be more easily kept
					cutbase2 = eventline + (evmin-eventline)*0.01*(100-j*percentcut)+base_tol*min(cuttime/2,3)
				if(cutbase2-smoothed[k]<0):
					if(endnumber<endgoal and k != end-1):
						endnumber += 1
					else:
						endnumber = 0
						endcut.append(k-endgoal)
						#If long enough, use increased tolerance to step backwards			
						if(cuttime > 2):
							endnumber = 0
							backwards = startcut[-1]
							endgoal = min(1+int(np.log(k-startcut[-1])),8)
							tempstart = startcut[-1]
							tempend = endcut[-1]
							done = False
							cutbase2 = eventline + (evmin-eventline)*0.01*(100-j*percentcut)+base_tol*min(cuttime/2,3)
							while(backwards>start and not done):
								backwards -= 1
								if(cutbase2-smoothed[backwards]<0):
									if(endnumber<endgoal):
										endnumber += 1
									else:
										done = True
										endnumber = 0
										newstart = backwards + endnumber +1
										neglength = 0
										if(newstart != tempstart):
											for i in range(len(startcut)):
												if(startcut[-i-1]>=newstart):
													neglength = -i-1
												else:
													if(newstart <= endcut[-i-1]):
														neglength = -i-1
														newstart = startcut[neglength]
													break
											startcut = startcut[0:neglength]
											endcut = endcut[0:neglength]
											startcut.append(newstart)
											endcut.append(tempend)
								else:
									endnumber = 0
				else:
					endnumber = 0
			else:
				print'endcut>startcut'
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
				if(HJD[endcut[l]]-HJD[startcut[l]]>1):
					print 'l = ',1,HJD[startcut[l]],HJD[endcut[l]]
					print 'area = ',area,' cutline = ',cutbase
				if(np.max(cutbase-smoothed[cutstart:cutend])>peak_tol):
					peaks.append(np.argmin(smoothed[cutstart:cutend])+cutstart)
					noburst = True
					print 'high peak'
					print 'starttime = ',HJD[cutstart],', endtime = ',HJD[cutend]
				elif(area>area_thresh):
					peaks.append(np.argmin(smoothed[cutstart:cutend])+cutstart)
					print 'area = ',area
					print 'cutstart = ',cutstart,'cutend = ',cutend
					print 'starttime = ',HJD[cutstart],', endtime = ',HJD[cutend]
					print 'cutbase = ',cutbase
#	Interpret the event based on the peaks
	npeak = len(peaks)
	if(npeak >= 1):
		print 'starttime = ',HJD[start],', endtime = ',HJD[end]
		print 'last area = ',area
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
		print leftmin
		print nineeight
		if((leftmin<nineeight) and not noburst):
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

			if(peakind+righttimecut-(start+lefttimecut)>100):
				base = np.percentile(smoothed[start+lefttimecut:peakind+righttimecut],80)
				leftcut = smoothed[start+lefttimecut:peakind]
				rightcut = smoothed[peakind:peakind+righttimecut]
				leftsub = abs(leftcut-base)
				rightsub = abs(rightcut-base)
				peakheight = np.percentile(smoothed[start+lefttimecut:peakind+righttimecut],80)-np.percentile(smoothed[start+lefttimecut:peakind+righttimecut],0)
				print 'eventline: ',eventline
				print '80: ',base,' leftmin: ',np.min(leftsub)/stdev,' rightmin: ',np.min(rightsub)/stdev
				print 'leftmin value: ',leftcut[np.argmin(leftsub)],' '
#				if(len(leftsub)<50 or len(rightsub)<50 or (base-np.percentile(leftcut,60))/peakheight>.8 or (base-np.percentile(rightcut,60))/peakheight>.8):
				if(len(leftsub)<50 or len(rightsub)<50):
					nobinary = True
				elif(np.min(leftsub)>stdev/2 or np.min(rightsub)>stdev/2 or (eventline-base)/peakheight<.08):
					base = np.percentile(smoothed[start+lefttimecut:peakind+righttimecut],60)
					leftsub = abs(leftcut-base)
					rightsub = abs(rightcut-base)
					print '60: ',base,' leftmin: ',np.min(leftsub)/stdev,' rightmin: ',np.min(rightsub)/stdev
					peakheight = np.percentile(smoothed[start+lefttimecut:peakind+righttimecut],60)-np.percentile(smoothed[start+lefttimecut:peakind+righttimecut],0)	
					if(np.min(leftsub)>stdev/2 or np.min(rightsub)>stdev/2 or (eventline-base)/peakheight<.08):
						base = np.percentile(smoothed[start+lefttimecut:peakind+righttimecut],40)
						leftsub = abs(leftcut-base)
						rightsub = abs(rightcut-base)
						print '40: ',base,' leftmin: ',np.min(leftsub)/stdev,' rightmin: ',np.min(rightsub)/stdev
						peakheight = np.percentile(smoothed[start+lefttimecut:peakind+righttimecut],40)-np.percentile(smoothed[start+lefttimecut:peakind+righttimecut],0)	
						if(np.min(leftsub)>stdev/2 or np.min(rightsub)>stdev/2 or (eventline-base)/peakheight<.08):
							nobinary = True
				if(not nobinary and len(leftcut)>50 and len(rightcut)>50):
					leftstart = np.argmin(leftsub)
					rightend = np.argmin(rightsub)
					leftcut = smoothed[start+lefttimecut+leftstart:peakind]
					rightcut = smoothed[peakind:peakind+rightend+1]
					print 'leftstart = ',leftcut[0],'rightend = ',rightcut[rightend]

					maxbinary = 0
	#				Cut off early so don't get fluctuation statistics at the bottom
					for m in range(21):
						print m*5,' left percent = ',np.percentile(leftcut,m*5),' right percent = ',np.percentile(rightcut,m*5)
						print m*5,' difference = ',(np.percentile(leftcut,m*5)-np.percentile(rightcut,m*5))/max(peakheight,50*stdev)
						leftdifference = (np.percentile(leftcut,m*5)-np.percentile(leftcut,max(m-1,0)*5))/peakheight
						rightdifference = (np.percentile(rightcut,m*5)-np.percentile(rightcut,max(m-1,0)*5))/peakheight
						if(abs((np.percentile(leftcut,m*5)-np.percentile(rightcut,m*5))/max(peakheight,50*stdev))>.05 and (leftdifference<.1 or rightdifference<.1)):
							posbinary += 1
							maxbinary = posbinary
						else:
							posbinary = 0
					if(maxbinary>1):
						events[pbinary] += 1
					else:
						if(evlength > 2000):
							events[longpspl] += 1
						else:
							events[pspl] += 1
				else:
					if(evlength > 2000):
						events[longpspl] += 1
					else:
						events[pspl] += 1
			else:
				events[pspl] += 1
	if(npeak >= 1):
		print 'event = ',i,', npeak = ',npeak
		print [HJD[peak] for peak in peaks]

#interpret lightcurve based on events
if(events[pspl]+events[binary]+events[burst]+events[pbinary]+events[square]+events[longpspl]>1):
	print 'periodic'
elif(events[pbinary]>0):
	print 'possible binary'
elif(events[pspl]==1):
	print 'pspl'
elif(events[longpspl]==1):
	print 'longpspl'
elif(events[binary]==1):
	print 'binary'
elif(events[burst]==1):
	print 'burst'
elif(events[square]==1):
	print 'square'
else:
	print 'flat'