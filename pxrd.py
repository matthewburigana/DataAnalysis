import csv
import os
from itertools import islice
import matplotlib.pyplot as plt

filepath = 'C:\\Users\\Matt\\OneDrive - McGill University\\Work\\Analysis\\PXRD\\'
 
# Lists to store the x and y data to be fit
xData = []
yData = []
yFit= []
numFiles = 0

# Loads a new file and clears any previously stored data if there is any.
def newFile():
    desiredFile = input('Enter the CSV file name: ')
    if (desiredFile == ''):
        for filename in os.listdir(filepath):
            processor(filename)
            plotXRD(numFiles)
    else:
        desiredFile = desiredFile + '.csv'
        processor(desiredFile)
        plotXRD(numFiles)
        
# Used by newFile() to process the file name passed to it.
def processor(name):
    global numFiles
    
    xData.append([])
    yData.append([])
    yFit.append([])
             
    with open(filepath + name, 'r') as inputData:
        reader = csv.reader(inputData, dialect = 'excel')
        for row in islice(reader, 1, None):
            xData[numFiles].append(float(row[0]))
            yData[numFiles].append(float(row[1]))
            if(len(row) > 2):
                yFit[numFiles].append(float(row[2]))
        
        maxInt = max(yData[numFiles])
        if(len(yFit[numFiles]) != 0):
            maxFit = max(yFit[numFiles]) 
        for i in range(len(xData[numFiles])):
            yData[numFiles][i] = (yData[numFiles][i]/maxInt) + (numFiles)
            if(len(yFit[numFiles]) != 0):
                yFit[numFiles][i] = (yFit[numFiles][i]/maxFit) + (numFiles)             
    numFiles += 1
    print('File %d: %s' % (numFiles, name))
        
# Plots the PXRD data of the input file number. If there are 10 or less files, the plot
# will be the color corresponding to the file number and if there are more than 10 loaded
# files, all spectra will be plotted as a single color. 
def plotXRD(fileNum = 1):
    fileNum -= 1
    plt.figure(1)
    if(len(yFit[fileNum]) != 0):
        plt.plot(xData[fileNum], yData[fileNum], 'bo', ms = 1.5)
        plt.plot(xData[fileNum], yFit[fileNum], 'r', linewidth = 1.5)
    else:
        if(numFiles <= 10):
            plt.plot(xData[fileNum], yData[fileNum], 'C%d' % fileNum, linewidth = 1.5)
        else:
            plt.plot(xData[fileNum], yData[fileNum], linewidth = 1.5)
    plt.xlabel('Scattering Angle (deg.)', fontsize = 'x-large')
    plt.ylabel('Intensity',  fontsize = 'x-large')
    plt.tick_params('y', labelleft = False)
    plt.tick_params(direction='in', labelsize = 'large', length = 5.0)
    plt.xlim(min(xData[fileNum]), max(xData[fileNum]))
    plt.show()
    
# Plots PXRD data with the fit and difference of the experimental data and Rietica fit.
def plotXRDError(fileNum = 1):
    fileNum -= 1
    error = []
    plt.figure(1)
    for i in range(len(yData[fileNum])):
        error.append(yData[fileNum][i] - yFit[fileNum][i] + fileNum)
    plt.plot(xData[fileNum], yData[fileNum], 'bo', ms = 1.5)
    plt.plot(xData[fileNum], yFit[fileNum], 'r', linewidth = 1.5)
    plt.plot(xData[fileNum], error, 'g', linewidth = 1.5)
    plt.xlabel('Scattering Angle (deg.)', fontsize = 'x-large')
    plt.ylabel('Intensity',  fontsize = 'x-large')
    plt.tick_params('y', labelleft = False)
    plt.tick_params(direction='in', labelsize = 'large', length = 5.0)
    plt.xlim(min(xData[fileNum]), max(xData[fileNum]))
    plt.show()
    
# Plots all loaded PXRD files as a stacked vertical plot. 
def plotAllXRD():
    plt.figure(1)
    for i in range(len(xData)):
        plt.plot(xData[i], yData[i], linewidth = 1.5, label = 'File %d' % (i+1))
    plt.xlabel('Scattering Angle (deg.)', fontsize = 'x-large')
    plt.ylabel('Intensity',  fontsize = 'x-large')
    plt.tick_params('y', labelleft = False)
    plt.tick_params(direction='in', labelsize = 'large', length = 5.0)
    plt.xlim(min(xData[i]), max(xData[i]))
    plt.show()
    
# Clears all stored data from the loaded files.
def clear():
    global numFiles
    numFiles = 0
    xData.clear()
    yData.clear()
    yFit.clear()
    
print('\nAvailable Functions:')
print('newFile()')
print('clear()')
print('plotXRD(fileNum = 1)')
print('plotXRDError(fileNum = 1')
print('plotAllXRD()\n')

newFile()