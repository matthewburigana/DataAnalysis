import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

filepath = 'C:\\Users\\Matt\\OneDrive - McGill University\\Work\\Temp\\'

xData = []
yData = []
xName = ''
yName = ''

def newFile():
    if (xData != []):
        xData.clear()
        yData.clear()
        
    filename = input('Enter the CSV file name: ') + '.csv'
    global xName
    global yName
    
    with open(filepath + filename, 'r') as inputData:
        reader = csv.reader(inputData, dialect = 'excel')
        rowNum = 1
        for row in reader:
            if(rowNum == 1):
                xName = row[0]
                yName = row[1]
                rowNum = 0
            else:
                xData.append(float(row[0]))
                yData.append(float(row[1]))

def regression(function):
    popt, pcov = curve_fit(function, xData, yData)
    print(popt)
    plt.plot(xData, yData, 'b-', label = 'Data')
    plt.plot(xData, function(xData, *popt), 'r-', label = 'Fit: a = %5.3f, b = %5.3f' % tuple(popt))
    plt.xlabel(xName)
    plt.ylabel(yName)
    plt.legend()
    plt.show()

newFile()