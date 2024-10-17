#!/usr/bin/pyhton3
# Extracellular Vesicle Image ANalysis
# Version 2.6
# Copyright © 2000 David Fabijan <david.fabijan@ijs.si> & Mario Kurtjak <mario.kurtjak@ijs.si> &
# This work is free. You can redistribute it and/or modify it under the
# terms of the WTFPL, Version 2.  See the COPYING file for more details.
# See http://www.wtfpl.net/ for more details.


#This program is designed to help with the manual classification of vesicles (or other notable features) imaged by AFM methods
#Internally the vesicles are called grains
#To use this program you must first prepare the input files.
#The program requires an AFM image in the .gwy or .xyz format
#It also requires that locations of the grains of interest in the .dat format as can be achieved with the Gwydion software
#The two files MUST have the same name (excluding the type) and be in the same folder as the program itself
#See the provided Example.dat and Example.gwy
#There can be more than one pair of files in the folder with EVIAN.py, in that case the pairs will be handled in order

#Built-in python libraries used
from dataclasses import dataclass
import json
import os

#Common pyhton libraries used
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import numpy as np
from scipy.interpolate import make_interp_spline


skipFirst=0 #Change this number if you exited before analysing all grains and would like to continue where you stopped. If you leave it 0 in such case, you will get duplicates in the output files.

#It is recomended that you install the gwyfile library if not this is the workaround
#Checks if it can read Gwyddion files directly. If not, another .xyz file with all coordinates needs to be exported from Gwyddion.
gwySupport=False
try:
    import gwyfile
    gwySupport=True
except:
    gwySupport = False

if gwySupport:
    print("GWYs are supported")
else:
    print("GWYs are NOT supported")

#Here we define the types of grains we want to classify
#typesOfGrains=["Round","Single lobed", "Concave","Multi lobed", "Flat","Neglect"] was the default selection used in our research
#but you can use any other array of descriptiors that would fit your particular problem
typesOfGrains=["Round","Single lobed", "Concave","Multi lobed", "Flat","Neglect"]
numberOfGrains=[0]*len(typesOfGrains)

#Setting start here
directOutput=False #Direct output turns off buttons and sorting. Change to True if you want to export all data and images without shape classification

zFilterValue=0.0 #Skip grains with max height less than this value in nanometers

#Some settings regaridng the plotting
spline=True
smooth=False
equalAspectRatio=False

eps=0.00001
epsDiag=0.1



def closeTo(a,b):
    return (abs(a-b)/a)<eps

def onPositiveDiagonal(ax,ay,bx,by):
    if ax==bx:
        if ay==by:
            return True
        return False
    d=(by-ay)/(bx-ax)
    return abs(d+1)<epsDiag

def onNegativeDiagonal(ax,ay,bx,by):
    if ax==bx:
        if ay==by:
            return True
        return False
    d=(by-ay)/(bx-ax)
    return abs(d-1)<epsDiag

def smoothing(arr):    
    l=len(arr)
    outputArr=np.zeros(l)
    if l<=2:
        return arr
    for i in range(l):
        if i==0:
            outputArr[i]=(arr[i]+arr[i+1]/2)*(2/3)
        elif i==l-1:
             outputArr[i]=(arr[i]+arr[i-1]/2)*(2/3)
        else:
            outputArr[i]=(arr[i]+(arr[i-1]+arr[i+1])/2)*(2/4)
    return outputArr



def hideButtons():
    for a in axs:
        a.set_visible(False)
    axexit.set_visible(False)
#    axneglect.set_visible(False)



def writeOutput(folderDesignation,i):
    jsonOutputData[int(i)]={"type":folderDesignation, "matrix":pointData.toImshowFormat(pointsForGrain)}


    with open(directoryname +"_"+ folderDesignation + "/" + str(i) + "_D2.tsv", "w") as output:
        for p in pointsForGrain:
            output.write(str(p))
    with open(directoryname +"_"+ folderDesignation + "/" + str(i) + "_sections.tsv", "w") as output:
        for k in range(4):
            output.write(str(labels[k]))
            for v in plotPoints[k]:
                output.write("\t" + str(v))
            output.write("\n")

    pat=directoryname+ "_"+folderDesignation+"_sizes.txt"
    check_file = os.path.isfile(pat)
    with open(pat, "a") as output:
        if not check_file:
            output.write("#Index\tMax z[nm]\tFeret[nm]\n")
        output.write(f"{i}\t{g.z*1e9}\t{g.diameter*1e9}\n") #Creates files with maximum heights and Feret sizes for each shape
    if "Neglect" not in folderDesignation:
     pat = directoryname + "_all_sizes.txt"
     nonNeglectDiam.append(g.diameter*1e9)
     nonNeglectHeight.append(g.z*1e9)
     nonNeglectRatio.append(g.z/g.diameter)
     check_file = os.path.isfile(pat)
     with open(pat, "a") as output:
         if not check_file:
             output.write("#Index\tMax z[nm]\tFeret[nm]\n")
         output.write(f"{i}\t{g.z*1e9}\t{g.diameter*1e9}\n") #Creates a file with maximum heights and Feret sizes for all grains that were not neglected

def generalClicked(event, typ,i):
    # print("General clicked ",typ)
    hideButtons()
    plt.savefig(directoryname + "_"+typ + "/" + str(i) + ".png")
    plt.close()
    idx=typesOfGrains.index(typ)
    numberOfGrains[idx]+=1
    writeOutput(typ,i)



def exitClicked(event):
    if len(jsonOutputData)>0:
        writeOutJsonAndCSV(mainInputeName)
    for idx, val in enumerate(typesOfGrains):
        print(val,": ",numberOfGrains[idx])
    nonNeglectDiam.sort()
    nonNeglectHeight.sort()
    nonNeglectRatio.sort()
    print("Max Martin diameter range: ",nonNeglectDiam[0]," - ",nonNeglectDiam[-1]," nm")
    print("Height range: ",nonNeglectHeight[0]," - ",nonNeglectHeight[-1]," nm")
    print("Aspect ratio (height/diameter) range: ",nonNeglectRatio[0]," - ",nonNeglectRatio[-1])
    exit()

def writeOutJsonAndCSV(fname):
    with open(fname+"dataDump.json", "w") as myf:
        json.dump(jsonOutputData, myf, indent = 1)

    with open(fname + "dataDump.tsv", "w") as myf:
        for (key,value) in jsonOutputData.items():
            outputStr = ""
            outputStr+=str(key)+"\t"
            outputStr += value["type"] + "\t"
            outputStr += str(len(value["matrix"])) + "\t"
            outputStr += str(len(value["matrix"][0]) )+'\t"'
            for a in value["matrix"]:
                for b in a:
                    outputStr += " "+str(b)
            outputStr +='"\n'
            myf.write(outputStr)



@dataclass
class pointData:
    x: float
    y: float
    z: float 

    def __str__(self):
        return "%g\t%g\t%g\n"%(self.x,self.y,self.z)

    #Important to allow us to sort the points onto a grid for plotting
    def __lt__(self, other):
        if abs(self.y-other.y)<eps:
            return self.x<other.x
        else:
            return self.y<other.y

    @staticmethod
    def toImshowFormat(points):
        a=points[0]
        arr=[]
        inerarr=[]
        for p in points:
            if a.y==p.y:
                inerarr.append(p.z)
            else:
                arr.append(inerarr)
                a=p
                inerarr=[p.z]
        if inerarr:#If not empty
            arr.append(inerarr)
        return arr
    
    @staticmethod
    def toSurfaceFormat(points):
        X=[]
        Y=[]
        Z=[]
        for p in points:
            X.append(p.x*1e9)
            Y.append(p.y*1e9)
            Z.append(p.z*1e9)
        return (X,Y,Z)


@dataclass
class grainData:
    x: float
    y: float
    z: float
    diameter: float 

    def minx(self):
        return self.x-self.diameter/2
    def maxx(self):
        return self.x+self.diameter/2
    def miny(self):
        return self.y-self.diameter/2
    def maxy(self):
        return self.y+self.diameter/2

    def __str__(self):
        return "[%g,%g,%g]"%(self.x,self.y,self.diameter)


#Here we start the main function
if __name__ == "__main__":
    thisDir=os.listdir()
    jsonOutputData={}
    gwys=[]
    nonNeglectDiam=[]
    nonNeglectHeight=[]
    nonNeglectRatio=[]

    for file in thisDir:
        if ".gwy" in file:
            if not(file[0]=="."): #This was added because Mac OS X adds hidden files starting with . that also contain .gwy in their names
                gwys.append(file[0:-4])
    print(gwys)

    #For each .gwy (or .xyz file)
    for mainInputeName in gwys:
        jsonOutputData = {}
        grainFile=mainInputeName+".dat"

        directoryname=grainFile+"_outputs"


        try:
            if directOutput:
                os.mkdir(directoryname+"_directOutput")
        except:
            pass

        # Create output directories
        for d in typesOfGrains:
            try:
                os.mkdir(directoryname+"_"+d)
            except:
                pass

        grains=[]
        points=[]

        #Read in the positions of the grains from the .dat file
        with open(grainFile) as fi:
            for line in fi:
                if "#" in line:
                    continue
                sp=line.split()
                grains.append(grainData(float(sp[0]),float(sp[1]),float(sp[2]),float(sp[4])))

        #Read in the full AFM image if in .xyz format
        if not gwySupport:
            xyzFile=mainInputeName+".xyz"
            with open(xyzFile) as fi:
                for line in fi:
                    if "#" in line:
                        continue
                    sp=line.split()
                    points.append(pointData(float(sp[0]),float(sp[1]),float(sp[2])))

        # Read in the full AFM image if in .gwy format
        else:
            gwyObj = gwyfile.load(mainInputeName+".gwy")
            channels = gwyfile.util.get_datafields(gwyObj)
            try:
                ZChannel = channels["ZSensor"] #Check the name of your z Channel in Gwyddion. If it is something else, modify the word in [] brackets
            except:
                print("Could not read the ZSensor chanel from the .gwy file. Check the name of your z Channel in Gwyddion. If it is not ZSensor, change it. Or modify the code at around line 312 in EVIAN.py")
                exit(1)

            xstep = ZChannel.xreal / ZChannel['xres']
            xstart = xstep / 2

            ystep = ZChannel.yreal / ZChannel['yres']
            ystart = ystep / 2

            data = ZChannel.data
            y = ystart
            for line in data:
                x = xstart
                for value in line:
                    points.append(pointData(float(x),float(y),float(value)))
                    x = x + xstep
                y = y + ystep

        #Grain counter
        i=0
        for g in grains:
            i+=1
            print(i)
            if i<skipFirst:
                continue
            pointsForGrain=[]

            plotPoints=[ [] for _ in range(4) ]
            labels=["y","x","+d","-d"]
            plotOffset=[0]*4

            #Limiting coordinates for grain
            mx=g.minx()
            Mx=g.maxx()
            my=g.miny()
            My=g.maxy()

            #Closest points to center of grain
            cx=0
            dist_cx=1e40
            cy=0
            dist_cy=1e40

            maxZ=0
            minZ=1e40

            for p in points:
                if (p.x>=mx and p.x<=Mx and p.y>=my and p.y<=My):
                    pointsForGrain.append(p)
                    if p.z>maxZ:
                        maxZ=p.z
                    if p.z<minZ:
                        minZ=p.z
                    if abs(p.x-g.x)<dist_cx:
                        cx=p.x
                        dist_cx=abs(p.x-g.x)
                    if abs(p.y-g.y)<dist_cy:
                        cy=p.y
                        dist_cy=abs(p.y-g.y)

            #Skip if below z filter
            if maxZ*1e9<zFilterValue:
                continue
            #Skip if only one point in grain
            if len(pointsForGrain)<2:
                continue

            pointsForGrain.sort()

            distanceFactor=abs(pointsForGrain[0].x-pointsForGrain[1].x)
            #Prepeare the cross sections (horizontal, vertical and both diagonals)
            for p in pointsForGrain:
                if closeTo(p.x,cx):
                    plotPoints[0].append(p.z)
                if closeTo(p.y,cy):
                    plotPoints[1].append(p.z)
                if onPositiveDiagonal(p.x,p.y,cx,cy):
                    plotPoints[2].append(p.z)
                if onNegativeDiagonal(p.x,p.y,cx,cy):
                    plotPoints[3].append(p.z)
                #Found the center
                if p.x==cx and p.y==cy:
                    for j in range(4):
                        plotOffset[j]=len(plotPoints[j])-1

            #Since we're going from top left to bottom right, we need to reverse the positive diagonal to fit
            plotPoints[2].reverse()
            plotOffset[2]=len(plotPoints[2])- plotOffset[2]-1

            #Prepare the figure and axes
            fig = plt.figure()
            ax0 = fig.add_subplot(2,2,1)
            ax3d=fig.add_subplot(2,2,2, projection='3d')
            ax1=fig.add_subplot(2,2,(3,4))
            whereWeAre= str(i) + '/' + str(len(grains))
            hdratio=(maxZ-minZ)/g.diameter
            fig.suptitle(whereWeAre + "; " + "height/diameter=" + str(hdratio))


            #Plot contours
            imshowData=pointData.toImshowFormat(pointsForGrain)
            ax0.imshow(imshowData,extent=(1e9*mx,1e9*Mx,1e9*my,1e9*My))

            #Plot 3D
            X,Y,Z=pointData.toSurfaceFormat(pointsForGrain)
            ax3d.plot_trisurf(X,Y,Z,cmap='viridis')


            ax1.set(xlabel="l [nm]", ylabel='z [nm]')

            #Add buttons
            numOfStates=len(typesOfGrains)
            axs=[]
            btns=[]
            for jj in range(numOfStates):
                axs.append(fig.add_axes([jj*0.9/numOfStates, 0.01,0.85/numOfStates, 0.03]))
                btns.append(Button(axs[jj], typesOfGrains[jj]))
                btns[jj].on_clicked(lambda ev,t=typesOfGrains[jj],iii=i: generalClicked(ev, t, iii))


            axexit = fig.add_axes([0.91, 0.01, 0.08, 0.03])
            bexit = Button(axexit, 'Exit')
            bexit.on_clicked(exitClicked)


            #Plot the cross sections
            for j in range(4):
                xx=np.array(range(-plotOffset[j],len(plotPoints[j])-plotOffset[j],1))
                if j<2:
                    #For x and y intersections
                    xx=xx*1e9*distanceFactor
                else:
                    #For diagonals
                    xx=xx*1e9*distanceFactor*1.41421356
                yy=np.array(plotPoints[j])*1e9
                if smooth:
                    xx=smoothing(xx)

                if spline and len(xx)>3: #Trying to do splines on too small a sample breaks things
                    X_Y_Spline = make_interp_spline(xx, yy)
                    X_ = np.linspace(xx.min(), xx.max(), 500)
                    Y_ = X_Y_Spline(X_)
                    ax1.plot(X_,Y_,label=labels[j])
                else:
                    ax1.plot(xx,yy,label=labels[j])

            if equalAspectRatio:
                ax1.set_aspect('equal')
            else:
                ax1.set_aspect('auto')
            ax1.legend()

            if not directOutput:
                aaa=plt.show()
            else:
                hideButtons()
                plt.savefig(directoryname+"_directOutput"+"/"+str(i)+".png")
                plt.close()
                writeOutput("directOutput",i)

        #At the end of gwy file handling
        writeOutJsonAndCSV(mainInputeName)

    exitClicked(0)