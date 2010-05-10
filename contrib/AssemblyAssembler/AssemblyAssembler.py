#!/usr/bin/python
#
#       AssemblyAssembler.pl
#
#		Conduct velvet assemblies across a range of kmer values, find parameter region 
#		with good assembly stats, conduct additional assemblies in that region and then
#		assemble contigs from all previous assemblies in one final assembly.  Potentially 
#		useful for de novo transcriptome assembly.  
#
#       Copyright 2010 Jacob Crawford <jc598@cornell.edu>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.
#
# ==========================================================================================

## import modules
import sys
import os

## Retrieve user input
if (len(sys.argv) == 9):
	print '\nRetrieving your input values'
	for i in  range(len(sys.argv)):
		if (sys.argv[i] == '-s'):
			if(int(sys.argv[i+1])%2 == 1):
				hashs = sys.argv[i+1]
			else:
				sys.exit('Invalid hash value.  Please enter odd number')
		elif (sys.argv[i] == '-e'):
			if(int(sys.argv[i+1])%2 == 1):
				hashe = sys.argv[i+1]
			else:
				sys.exit('Invalid hash value.  Please enter odd number')
		elif (sys.argv[i] == '-f'):
			## Velvethcall must be entered as a string at the command line with ' ' or " "
			velvethcall = sys.argv[i+1]
		elif (sys.argv[i] == '-v'):
			velvDir = sys.argv[i+1]
	print str('Okay, will conduct assemblies across range of kmer values from '+hashs+' to '+hashe+'\n')
else: 
	sys.exit('Incorrect command line entry.  Please enter parameters values for -s, -e, -f -v and -d')

## Check and set hash value range
print 'Checking your hash value range'
if ((int(hashe)-int(hashs))>=16):
	print 'Looks good\n'
	hashrange = range(int(hashs),(int(hashe)+1),2)
	hashlist = [hashrange[0],hashrange[int(round(len(hashrange)*0.12))],hashrange[int(round(len(hashrange)*0.25))],hashrange[int(round(len(hashrange)*0.37))],hashrange[int(round(len(hashrange)*0.5))],hashrange[int(round(len(hashrange)*0.6))],hashrange[int(round(len(hashrange)*0.7))],hashrange[int(round(len(hashrange)*0.82))],hashrange[len(hashrange)-1]]
else:
	sys.exit('Insufficient hash range. This routine is only useful if kmer range >= 16. \n\n')
	
## Make directories, one for each hash value
parDir = os.getcwd()
try:
	os.system(str('rm -rf SceneOfTheCrime'))
	os.system(str('rm -rf FinalDir'))
	os.mkdir('SceneOfTheCrime')
	os.mkdir('FinalDir')
	os.chdir(str(os.getcwd()+'/SceneOfTheCrime'))
except:
	os.mkdir('SceneOfTheCrime')
	os.mkdir('FinalDir')
	os.chdir(str(os.getcwd()+'/SceneOfTheCrime'))	

hashes = [hashlist[0],hashlist[2],hashlist[4],hashlist[6],hashlist[8]]
for i in range(5):
	os.mkdir(str('Dir'+str(hashes[i])))

## Call velveth and velvetg for hash values across coarse sampling of range
print 'Running velvet with intitial set of hash values \n'
contigs = []
for i in range(5):
	os.chdir(str(os.getcwd()+'/Dir'+str(hashes[i])))
	print 'Running velveth with k= '+str(hashes[i])
	os.system(str(velvDir+'/velveth ./ '+str(hashes[i])+' '+velvethcall+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	print 'Running velvetg on k= '+str(hashes[i])+' velveth files'
	os.system(str(velvDir+'/velvetg ./'+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	os.system(str('mv contigs.fa '+parDir+'/FinalDir/contigs'+str(hashes[i])+'.fa'))
	os.system(str('mv Log '+parDir+'/FinalDir/Log'+str(hashes[i])))
	contigs.append(str('./contigs'+str(hashes[i])+'.fa '))
	os.chdir('..')
	os.system(str('rm -rf Dir'+str(hashes[i])))
	
## Read log files for first 5 assemblies and determine which is best
MAXCONT = {}
for j in range(5):
	input = open(str(parDir+'/FinalDir/Log'+str(hashes[j])),'r')
	for i in input:
		i = i.split(' ')
		if (i[0] == 'Final'):
			if not MAXCONT.has_key(hashes[j]):
				MAXCONT[hashes[j]] = 0
			MAXCONT[hashes[j]] = int(i[10][:-1])
	input.close()

BEST = max(MAXCONT, key = MAXCONT.get)
print str('Of the first 5 assemblies k='+str(BEST)+' had the longest contig')
print str('Sampling deeper around '+str(BEST)+'\n')
## Make directories and do additional velvet runs
if (BEST == hashlist[0]):
	hashes.append(hashlist[1])
	os.mkdir(str('Dir'+str(hashlist[1])))
	os.chdir(str(os.getcwd()+'/Dir'+str(hashlist[1])))
	os.system(str(velvDir+'/velveth ./ '+str(hashlist[1])+' '+velvethcall+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	os.system(str(velvDir+'/velvetg ./'+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	os.system(str('mv contigs.fa '+parDir+'/FinalDir/contigs'+str(hashlist[1])+'.fa'))
	os.system(str('mv Log '+parDir+'/FinalDir/Log'+str(hashlist[1])))
	contigs.append(str('./contigs'+str(hashlist[1])+'.fa '))
	os.chdir('..')
	os.system(str('rm -rf Dir'+str(hashlist[1])))
	
elif (BEST == hashlist[8]):
	hashes.append(hashlist[7])
	os.mkdir(str('Dir'+str(hashlist[7])))
	os.chdir(str(os.getcwd()+'/Dir'+str(hashlist[7])))
	os.system(str(velvDir+'/velveth ./ '+str(hashlist[7])+' '+velvethcall+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	os.system(str(velvDir+'/velvetg ./'+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	os.system(str('mv contigs.fa '+parDir+'/FinalDir/contigs'+str(hashlist[7])+'.fa'))
	os.system(str('mv Log '+parDir+'/FinalDir/Log'+str(hashlist[7])))
	contigs.append(str('./contigs'+str(hashlist[7])+'.fa '))
	os.chdir('..')
	os.system(str('rm -rf Dir'+str(hashlist[7])))
	
else:
	for i in range(9):
		if (hashlist[i] == BEST):
			HASHUP = hashlist[i+1]
			HASHDWN = hashlist[i-1]
	hashes.append(HASHUP)
	hashes.append(HASHDWN)
	os.mkdir(str('Dir'+str(HASHUP)))
	os.chdir(str(os.getcwd()+'/Dir'+str(HASHUP)))
	os.system(str(velvDir+'/velveth ./ '+str(HASHUP)+' '+velvethcall+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	os.system(str(velvDir+'/velvetg ./'+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	os.system(str('mv contigs.fa '+parDir+'/FinalDir/contigs'+str(HASHUP)+'.fa'))
	os.system(str('mv Log '+parDir+'/FinalDir/Log'+str(HASHUP)))
	contigs.append(str('./contigs'+str(HASHUP)+'.fa '))
	os.chdir('..')
	os.system(str('rm -rf Dir'+str(HASHUP)))
	os.mkdir(str('Dir'+str(HASHDWN)))
	os.chdir(str(os.getcwd()+'/Dir'+str(HASHDWN)))
	os.system(str(velvDir+'/velveth ./ '+str(HASHDWN)+' '+velvethcall+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	os.system(str(velvDir+'/velvetg ./'+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	os.system(str('mv contigs.fa '+parDir+'/FinalDir/contigs'+str(HASHDWN)+'.fa'))
	os.system(str('mv Log '+parDir+'/FinalDir/Log'+str(HASHDWN)))
	contigs.append(str('./contigs'+str(HASHDWN)+'.fa '))
	os.chdir('..')
	os.system(str('rm -rf Dir'+str(HASHDWN)))
	
os.chdir(str(parDir+'/FinalDir'))
os.system(str(velvDir+'/velveth ./ 35 '+'-fasta -long '+''.join(contigs)+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
os.system(str(velvDir+'/velvetg ./'+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
os.system('mv contigs.fa Finalcontigs.fa')
	
## Remove interim directories and their contents
os.chdir('..')
os.system('rm -rf SceneOfTheCrime')
	
print 'Assembly Assembler Job Complete \n' 

input = open(str(parDir+'/FinalDir/Log'),'r')
for i in input:
	i = i.split(' ')
	if (i[0] == 'Final'):
		print ' '.join(i)
input.close()
























