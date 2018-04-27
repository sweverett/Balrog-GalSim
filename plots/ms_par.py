"""Finds the FOF groups changed after Balrog-injection and unchanged after Balrog-injection.
Author: Brian Yanny"""

#!/usr/bin/env python

import sys
import csv

if(len(sys.argv) != 3):
  sys.exit("Usage:  python par.py matchlist outprefixname")

fn=sys.argv[1]
outprefix=sys.argv[2]

gidlist={} 
gidsizelist={}
fofidlist={}
dumplist={}

othergidlist={}

otherlist={}
othersizelist={}

outok=open(outprefix+'.ok','w')
outrerunmof=open(outprefix+'.rerun','w')

with open(fn) as csvfile:
	reader=csv.DictReader(csvfile)
        for row in reader:
		num=row['number_1_1']
		gid=row['GroupID_1']
		gsize=row['GroupSize_1']
		if (gid==""):
		  gid=0
		  gsize=0
		othergid=row['GroupID_2']
		othersize=row['GroupSize_2']
		if (othergid==""):
		  othergid=0
	          othersize=0
		fofid=row['fofid_1']
		#print "num:",num," gid:",gid,gsize,othergid,othersize,fofid
                if (gid!=0):
		  if gid not in gidlist:
		    if(gsize != othersize):
		      outrerunmof.write(str(num)+' '+str(fofid)+'\n')
		      dumplist[int(num)] = 1
		      #print "group size mismatch"
		    else:
		      gidlist[gid]=[int(num)]
		      gidsizelist[gid]=int(gsize)
		      fofidlist[gid]=int(fofid)
		      othergidlist[gid]=othergid
		      otherlist[othergid]=[int(num)]
		      othersizelist[othergid]=int(othersize)
	          else:	
                      if (othergid != othergidlist[gid]):
		        b=gidlist[gid]
                        for a in range(0,len(b)):
			  tnum=int(b[a])
			  if tnum not in dumplist:
		            outrerunmof.write(str(b[a])+' '+str(fofid)+'\n')
		            dumplist[int(b[a])] = 1
			if num not in dumplist:
		          outrerunmof.write(str(num)+' '+str(fofid)+'\n')
		          dumplist[int(num)] = 1
		        #print "otherlist id mismatch"
		      else:
                        gidlist[gid].append(int(num))
			otherlist[othergid].append(int(num))
			if len(gidlist[gid]) == gidsizelist[gid]:
			  b=gidlist[gid]
                          for a in range(0,len(b)):
			    tnum=int(b[a])
			    if tnum not in dumplist:
		  	      outok.write(str(b[a])+' '+str(fofid)+'\n')
                              dumplist[int(b[a])] = 1
		else: 
		  if (othergid==0):
		    outok.write(str(num)+' '+str(fofid)+'\n')
                    dumplist[int(num)]=1
 	          else:
		    if num not in dumplist:
		      outrerunmof.write(str(num)+' '+str(fofid)+'\n')
	              dumplist[int(num)]=1
		    #print "othergid not 0"

for q in gidlist.keys():
	if(len(gidlist[q]) != gidsizelist[q]):
	    b=gidlist[q]
	    #print q,b
	    for a in range(0,len(b)):
	      num=int(b[a])
              if num not in dumplist:
	        outrerunmof.write(str(num)+' '+str(fofidlist[q])+'\n')
                dumplist[int(num)]=1
		
outrerunmof.close()
outok.close()
