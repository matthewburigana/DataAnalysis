import csv
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
j = 1j
 
filepath = 'C:\\Users\\Matt\\OneDrive - McGill University\\Work\\Temp\\'
 
# Lists to store the x and y data to be fit
xData = []
yData = []
xName = ''
yName = ''
 
# Loads a new file and clears any previously stored data if there is any. regression is then
# called on the desired function to be fit to the data. 
def newData(function):
    desiredFile = input('Enter the CSV file name: ')
    if (desiredFile == ''):
        for filename in os.listdir(filepath):
            processor(filename, function)
    else:
        
        desiredFile = desiredFile + '.csv'
        processor(desiredFile, function)
             
# Used by newFile() to process the file name passed to it and do a regression on the data
# with the input function
def processor(name, function):
    print('File: ' + name)
    if (xData != []):
        xData.clear()
        yData.clear()
                 
    global xName
    global yName
             
    with open(filepath + name, 'r') as inputData:
        reader = csv.reader(inputData, dialect = 'excel')
        firstRow = True
        for row in reader:
            if(firstRow):
                xName = row[0]
                yName = row[1]
                firstRow = False
            else:
                xData.append(float(row[0]))
                yData.append(float(row[1]))
             
    regression(function)
     
def regression(function):
    xTemp = np.ndarray(shape = (len(xData),), dtype = float, buffer = np.array(xData))
    popt, cov1 = leastsq(residual, (1, 1, 1, 1, 1, 1), args = (xTemp, yData), maxfev = 100000)
    #popt2, cov2 = leastsq(residual2, (popt1[0], popt1[1], 1, 1), args = (xTemp, yData), maxfev = 100000)
    #popt3, cov3 = leastsq(residual3, (popt2[0], popt2[1], popt2[2], 1, popt2[3], 1), args = (xTemp, yData), maxfev = 100000)
    fit = 'Fit: '
    for i in range(len(popt)):
        if (i < len(popt)-1):
            fit = fit + chr(ord('a') + i) + ' = %5.3f, ' % popt[i]
        else:
            fit = fit + chr(ord('a') + i) + ' = %5.3f' % popt[i]
    print(fit)
    plt.plot(xData, yData, 'b-', label = 'Data')
    plt.plot(xData, function(xTemp, *popt), 'r-', label = fit)
    plt.xlabel(xName)
    plt.ylabel(yName)
    plt.legend()
    plt.show()
     
def residual(params, x, Z):
    a, b, c, d, e, f = params
    diff = Z - func(x, a, b, c, d, e, f)
    return diff.imag**2 + diff.real**2

def residual1(params, x, Z):
    a, b = params
    diff = Z - func(x, a, b, 1, 1, 1, 1)
    return diff.astype(np.complex128).view(np.float64)

def residual2(params, x, Z):
    a, b, c, e = params
    diff = Z - func(x, a, b, c, 1, e, 1)
    return diff.astype(np.complex128).view(np.float64)

def residual3(params, x, Z):
    a, b, c, d, e, f = params
    diff = Z - func(x, a, b, c, d, e, f)
    return diff.astype(np.complex128).view(np.float64)
 
def func(x, a, b, c, d, e, f):
    Z = a*x*2.0*np.pi*j + b/(b*c*(x*2.0*np.pi*j)**d + 1.0) + 1.0/(e*(x*2.0*np.pi*j)**f)
    return Z
 
print('\nAvailable Function:')
print('newData(function)\n')
 
newData(func)