import csv
import ternary
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

filepath = 'C:\\Users\\McCalla Lab\\Documents\\Data\\Matthew\\'

# Empty global lists to store the fractions of each of the three elements and their respective capacities.
# The element names are then stored in the global xName, yName, and zName strings
xVal = []
yVal = []
zVal = []
capacities = []
xName = ''
yName = ''
zName = ''

# Clears data from a previous file if it is present then loads the new data
def newFile():
    filename = input('Enter the CSV file name: ') + '.csv'
    if (xVal != []):
        xVal.clear()
        yVal.clear()
        zVal.clear()
        capacities.clear()
        
    global xName
    global yName
    global zName
    
    with open(filepath + filename, 'r') as inputData:
        reader = csv.reader(inputData, dialect = 'excel')
        rowNum = 1
        for row in reader:
            if (rowNum == 1):
                xName = row[0]
                yName = row[1]
                zName = row[2]
                rowNum = 0
            else:
                if (len(row) == 3):
                    xVal.append(float(row[0]))
                    yVal.append(float(row[1]))
                    zVal.append(float(row[2]))
                    capacities.append(0)
                else:
                    xVal.append(float(row[0]))
                    yVal.append(float(row[1]))
                    zVal.append(float(row[2]))
                    capacities.append(float(row[3]))
    plotDiagram()            
    
# Plots the ternary data using python-ternary
def plotDiagram():
    points = []
    colors = []
    
    # Creates tuples that will be used as data points
    for i in range(len(xVal)):
        points.append((xVal[i], yVal[i], zVal[i]))
        
    # Assigns colors to each composition based on the capacity found. The highest capacity is
    # given the color red (1, 0, 0) then the color goes down to black (0, 0, 0) for the lowest
    # capacity. 
    for cap in capacities:
        if(min(capacities) != max(capacities)):
            colors.append((((cap-min(capacities))/(max(capacities)-min(capacities))), 0.0, 0.0))
        else:
            colors.append((1,0,0))

    # Creates a ternary set of axes to plot the diagram
    figure, tax = ternary.figure(scale = 1.0)
    tax.boundary()
    tax.gridlines(multiple = 0.1, color = 'blue')
    if(min(capacities) != max(capacities)):
        figure.set_size_inches(6.5,7.06)
        cm = LinearSegmentedColormap.from_list('Capacities', ['black', 'red'], N = 1024)
        cb_args = {'shrink':0.5, 'aspect':20, 'pad':0.02, 'orientation':'horizontal'}
        tax.scatter(points, marker = 'o', c = colors, s = 20, colorbar = True, colormap = cm, vmin = min(capacities), vmax = max(capacities), cbarlabel = 'Capacity (mAh/g)', cb_kwargs = cb_args)
        tax.left_axis_label(zName, fontsize = 14, offset = 0.11)
        tax.right_axis_label(yName, fontsize = 14, offset = 0.105)
        tax.bottom_axis_label(xName, fontsize = 14, offset = 0.007)
    else:
        figure.set_size_inches(6.5,6)
        tax.left_axis_label(zName, fontsize = 14, offset = 0.11)
        tax.right_axis_label(yName, fontsize = 14, offset = 0.105)
        tax.bottom_axis_label(xName, fontsize = 14, offset = 0.005)
        tax.scatter(points, marker = 'o', c = colors, s = 20, label = '%.4f mAh/g' % max(capacities))
        if (max(capacities) != 0):
            tax.legend()
    tax.ticks(axis = 'blr', linewidth = 1.0, multiple = 0.1, tick_formats = '%.1f', offset = 0.02)
    plt.axis('off')
    print('Ternary plots generated using the python-ternary library\nMarc Harper et al. (2015). python-ternary: Ternary Plots in Python. Zenodo. 10.5281/zenodo.34938')
    tax.show()
        
print('\nAvailable Functions:')
print('newFile()')
print('plotDiagram()\n')

newFile()