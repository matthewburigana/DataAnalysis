import csv
from itertools import islice
from matplotlib import rc, rcParams
import matplotlib.pyplot as plt
import numpy as np

filepath = 'C:\\Users\\Matt\\OneDrive - McGill University\\Work\\Cells\\'

# Empty global lists to store the cycle numbers, times, and voltages in chronological order and 
# charge and discharge capacities and resistances
cycles = []
times = []
voltages = []
chgCapacities = []
dchgCapacities = []
resistance = []
meanResistances = []
stdevResistances = []

# A mass dictionary for the mass of each channel
# mass = {x: 0.0 for x in range(1, 65)}
mass = []

# A dictionary to store the currents from each channel where the channel
# number is the key and the value is a list of currents
# channels = {x: [] for x in range(1, 65)}
channels = []

numFiles = 0

# Initializes a new file by clearing all previously stored data then processing the data stored
# by the new input file
def newFile():
    global numFiles
    filename = input('Enter the CSV file name: ') + '.csv'
    massFile = input('Enter the mass file name: ') + '.csv'

    cycles.append([])
    times.append([])
    voltages.append([])
    chgCapacities.append([])
    dchgCapacities.append([])
    resistance.append([])
    meanResistances.append([])
    stdevResistances.append([])
    mass.append({x: 0.0 for x in range(1, 65)})
    channels.append({x: [] for x in range(1, 65)})  
            
    with open(filepath + filename, 'r') as inputData:
        reader = csv.reader(inputData, dialect = 'excel')
        
        # Appends the cycle numbers, times, and average voltages to their respective
        # list for each row in the CSV file
        for row in islice(reader, 8, None):
            if(row[0] == ''):
                break
            cycles[numFiles].append(int(row[0]))
            times[numFiles].append(float(row[1]))
            voltages[numFiles].append(float(row[5]))
        
            # Appends the current for each channel in the current list at each channel
            # in the channels dictionary
            for j in range(1, 65):
                channels[numFiles][j].append(float(row[j+5]))

    # Appends empty lists to store resistances and capacities when they are calculated
    for _ in range(int(max(cycles[numFiles])/2)):
        resistance[numFiles].append([])
        chgCapacities[numFiles].append([])
        dchgCapacities[numFiles].append([])
        
    # Reads masses stored in the same file format as an exported file
    if(massFile != ''):
        with open(filepath + massFile, 'r') as masses:
            massRead = csv.reader(masses, dialect = 'excel')
            for row in islice(massRead, 1, None):
                for i in range(1,65):
                    mass[numFiles][i] = float(row[i+1])
                break
    numFiles += 1
            
def clear():
    global numFiles
    cycles.clear()
    times.clear()
    voltages.clear()
    for file in range(numFiles):
        for key in channels[file]:
            mass[file][key] = 0.0
            channels[file][key].clear()
    chgCapacities.clear()
    dchgCapacities.clear()
    resistance.clear()
    meanResistances.clear()
    stdevResistances.clear()
    numFiles = 0
    
# Set the mass of the sample on a channel
def setMass(channel, fileNum = 1):
    fileNum -= 1
    newMass = float(input('Enter the new sample mass for channel %d in grams: ' % channel))
    if(chgCapacities[fileNum][0] != []):
        for i in range(min(cycles[fileNum]), int(max(cycles[fileNum])/2) + 1):
            chgCapacities[fileNum][i-1][channel-1] = chgCapacities[fileNum][i-1][channel-1]*mass[fileNum][channel]/newMass
            dchgCapacities[fileNum][i-1][channel-1] = dchgCapacities[fileNum][i-1][channel-1]*mass[fileNum][channel]/newMass
    mass[fileNum][channel] = newMass
    
# Exports the masses, charge and discharge capacities and resistances to a CSV file
def export(fileNum = 1):
    fileNum -= 1
    outputFile = input('Enter the export file name: ') + '.csv'
    allCapacities()
    allResistances()
    with open(filepath + outputFile, 'w', newline = '') as output:
        writer = csv.writer(output, dialect = 'excel')
        channel = []
        channelMass = []
        for i in range(1,65):
            channel.append('Channel %d' % i)
            channelMass.append(mass[fileNum][i])
        writer.writerow(['Cycle', 'Name'] + channel)
        writer.writerow(['','Mass (g)'] + channelMass)
        for cycle in range(min(cycles[fileNum]), int(max(cycles[fileNum])/2) + 1):
            writer.writerow([cycle, 'Charge (mAh/g)'] + chgCapacities[fileNum][cycle-1])
            writer.writerow(['', 'Discharge (mAh/g)'] + dchgCapacities[fileNum][cycle-1])
            writer.writerow(['', 'Resistance'] + resistance[fileNum][cycle-1])
            writer.writerow(['', 'Mean Resistance'] + [meanResistances[fileNum][cycle-1]])
            writer.writerow(['', 'Standard Deviation'] + [stdevResistances[fileNum][cycle-1]])
        
# Calculates the charge and discharge capacity of a channel during the specified cycle
# between the minimum and maximum voltages. Integrates dQ = Idt
def capacity(channel, cycle = 1, fileNum = 1):
    fileNum -= 1
    qCharge = 0.0
    qDischarge = 0.0
    if(mass[fileNum][channel] == 0.0):
        mass[fileNum][channel] = float(input('Enter the sample mass for channel %d in grams: ' % channel))
    for i in range(cycles[fileNum].index(2*cycle - 1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle)):
        
        # Only adds to charge capacity during a charging cycle and only if the voltages are 
        # between those specified. The average of the capacity and the next capacity are taken
        # and multiplied by the difference in time of those measurements
        if(2*cycle - 1 == cycles[fileNum][i] and i != len(cycles[fileNum]) - 1):
            qCharge += 0.5*(channels[fileNum][channel][i] + channels[fileNum][channel][i+1])*(times[fileNum][i+1] - times[fileNum][i])
            
        # Adds to discharge capacity during a discharge cycle if the voltages are within the 
        # specified voltages and the iterator is not at the last value of the cycles
        elif(2*cycle == cycles[fileNum][i] and i != len(cycles[fileNum]) - 1):
            qDischarge += 0.5*(channels[fileNum][channel][i] + channels[fileNum][channel][i+1])*(times[fileNum][i+1] - times[fileNum][i])
                
    # Converts the capacities to Ah/g then returns a dictionary listing the channel, cycle,
    # charge and discharge capacities
    qCharge = round(qCharge/mass[fileNum][channel]/1000.0, 8)
    qDischarge = round(qDischarge/mass[fileNum][channel]/1000.0, 8)
    if(len(chgCapacities[fileNum][cycle-1]) < 64):
        chgCapacities[fileNum][cycle-1].append(qCharge)
        dchgCapacities[fileNum][cycle-1].append(qDischarge)
    return({'Channel': channel, 'Cycle': cycle, 'Charge capacity (mAh/g)': qCharge, 'Discharge capacity (mAh/g)': qDischarge})

# Prints the capacity of all 64 channels in the input cycle
def capacity64(cycle = 1, fileNum = 1):
    for i in range(1,65):
        print(capacity(i, cycle, fileNum))
      
# Prints the capacity of each cycle for the input channel
def capacityAllCycles(channel, fileNum = 1):
    for i in range(min(cycles[fileNum-1]), int(max(cycles[fileNum-1])/2) + 1):
        print(capacity(channel, i, fileNum))
    
# Prints the capacities of all the channels for all the cycles    
def allCapacities(fileNum = 1):
    for i in range(1,65):
        for j in range(min(cycles[fileNum-1]), int(max(cycles[fileNum-1])/2) + 1):
            print(capacity(i, j, fileNum))

# Calculates then prints the average voltages of each channel using Vavg = VdQ/Q
def averageVoltages(startTime, endTime, fileNum = 1):
    fileNum -= 1
    capacities = {}
    for i in range(1, 65):
        capacities[i] = []
        voltage = 0.0
        
        # Appends the specific capacity at each time to the capacities dictionary then calculates the 
        # differential voltage and adds that to the sum of VdQ then divides that sum by the total capacity
        # to give the average voltage
        for j in range(len(times[fileNum])):
            if(times[fileNum][j] >= startTime and times[fileNum][j] <= endTime and j < len(times[fileNum]) - 1):
                capacities[i].append(0.5*(channels[fileNum][i][j] + channels[fileNum][i][j+1])*(times[fileNum][j+1] - times[fileNum][j])/1000000.0)
        totalCapacity = sum(capacities[i])
        for k in range(len(times[fileNum])):
            if(times[fileNum][k] >= startTime and times[fileNum][k] <= endTime and k < len(capacities[i]) - 1):
                voltage += 0.5*(voltages[fileNum][k] + voltages[fileNum][k+1])*(capacities[i][k+1] - capacities[i][k])/totalCapacity
        print('Channel %d average voltage: %f V' % (i, voltage))
    
# Prints the resistances of all 64 channels based on the linear regression of 
# voltage vs current for the charging cycle
def resistances(cycle = 1, fileNum = 1):
    fileNum -= 1
    for i in range(1, 65):
        cycleCurrents = []
        cycleVoltages = []
        for j in range(cycles[fileNum].index(2*cycle-1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle-1)):
            cycleCurrents.append(channels[fileNum][i][j]/1000000.0)
            cycleVoltages.append(voltages[fileNum][j])
        param = np.polyfit(cycleCurrents, cycleVoltages, 1)
        print('Channel %d resistance: %f Ohms' % (i, param[0]))
        if(len(resistance[fileNum][cycle-1]) < 64):
            resistance[fileNum][cycle-1].append(param[0])
        
    # Calculates and prints the mean, standard deviation and variance of the 64 resistance
    # values
    stdev = np.std(resistance[fileNum][cycle-1], ddof = 1)
    variance = np.var(resistance[fileNum][cycle-1], ddof = 1)
    mean = np.mean(resistance[fileNum][cycle-1])
    if(len(meanResistances[fileNum]) < max(cycles[fileNum])/2):
        meanResistances[fileNum].append(mean)
        stdevResistances[fileNum].append(stdev)
    print('Average: %f' % mean)
    print('Standard Deviation: %f' % stdev)
    print('Variance: %f' % variance)

# Prints resistances for ever channel during each cycle
def allResistances(fileNum = 1):
    for i in range(min(cycles[fileNum-1]), int(max(cycles[fileNum-1])/2) + 1):
        print('\nCycle %d' % i)
        resistances(i, fileNum)
    
# Plots the current vs voltage for one channel during the input cycle
def plotCurrentVsVolts(channel, xMin = 3.1, xMax = 4.5, cycle = 1, fileNum = 1):
    fileNum -= 1
    chargeVoltages = []
    chargeCurrents = []
    dischargeVoltages = []
    dischargeCurrents = []
    rc('font', weight = 'bold')
    plt.figure(1)
    ax = plt.subplot(111)
    
#   Loops through all the voltages for the specified cycle and stores them in the 
#   chargeVoltages list for the charge cycle and dischargeVoltages for the discharge 
#   cycle. Loops through the currents for the specified cycle stored in the channels 
#   dictionary at the specified channel and stores the charge currents in the 
#   chargeCurrents list and discharge currents in the dischargeCurrents list.
    if cycle*2 not in cycles[fileNum]:
        for i in range(cycles[fileNum].index(2*cycle - 1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle-1)):
            chargeVoltages.append(voltages[fileNum][i])
            chargeCurrents.append(channels[fileNum][channel][i])
        chargePlot, = plt.plot(chargeVoltages, chargeCurrents, 'b', linewidth = 2.0, label = 'Charge')
        dischargePlot, = plt.plot(dischargeVoltages, dischargeCurrents, 'r', linewidth = 2.0, label = 'Discharge')
        
        # Positions the legend to the top left corner outside the plot
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(handles = [chargePlot, dischargePlot], loc = 2, bbox_to_anchor = (1,1), fontsize = 'small')
        plt.title('Current on Channel %d During Cycle %d' % (channel, cycle))
        plt.xlabel('Voltage (V)')
        plt.ylabel('Current (\u03BCA)')
        plt.tick_params(direction='in')
        plt.xlim((xMin, xMax))
    else:
        for i in range(cycles[fileNum].index(2*cycle - 1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle)):
            if(cycles[fileNum][i] == 2*cycle - 1):
                chargeVoltages.append(voltages[fileNum][i])
                chargeCurrents.append(channels[fileNum][channel][i])
            else:
                dischargeVoltages.append(voltages[fileNum][i])
                dischargeCurrents.append(channels[fileNum][channel][i])
         
    # Plots the charge currents vs charge voltages in blue then the discharge currents vs 
    # discharge voltages in red
        chargePlot, = plt.plot(chargeVoltages, chargeCurrents, 'b', linewidth = 2.0, label = 'Charge')
        dischargePlot, = plt.plot(dischargeVoltages, dischargeCurrents, 'r', linewidth = 2.0, label = 'Discharge')
        
        # Positions the legend to the top left corner outside the plot
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(handles = [chargePlot, dischargePlot], loc = 2, bbox_to_anchor = (1,1), fontsize = 'small')
        plt.title('Current on Channel %d During Cycle %d' % (channel, cycle), fontsize = 'xx-large', weight = 'bold')
        plt.xlabel('Voltage (V)', fontsize = 'x-large', weight = 'bold')
        plt.ylabel('Current (\u03BCA)', fontsize = 'x-large', weight = 'bold')
        plt.tick_params(direction='in', labelsize = 'large', length = 5.0, width = 1.5, top = True, right = True)
        plt.xlim((xMin, xMax))
    plt.show()
    plt.close()

# Plots the current vs voltage for all 64 channels in one window
def plotCurrentVsVolts64(cycle = 1, xMin = 3.1, xMax = 4.5, fileNum = 1):
    fileNum -= 1
    rc('font', weight = 'bold')
    plt.figure(1)   

    # Loops through the 64 channels and adds each new plot as a subplot to figure 1,
    # using the same plotting logic as the plot function    
    for i in range(1,65):
        loc = 0
        if(i <= 8):
            loc = i + (7*(i-1))
        elif(i <= 16):
            loc = i + (7*(i-10))
        elif(i <= 24):
            loc = i + (7*(i-19))
        elif(i <= 32):
            loc = i + (7*(i-28))
        elif(i <= 40):
            loc = i + (7*(i-37))
        elif(i <= 48):
            loc = i + (7*(i-46))
        elif(i <= 56):
            loc = i + (7*(i-55))
        else:
            loc = i + (7*(i-64))
        ax = plt.subplot(8, 8, loc, title = 'Channel %d' % i)
        ax.title.set_visible(False)
        chargeVoltages = []
        chargeCurrents = []
        dischargeVoltages = []
        dischargeCurrents = []
        if cycle*2 not in cycles[fileNum]:
            for j in range(cycles[fileNum].index(2*cycle-1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle-1)):
                chargeVoltages.append(voltages[fileNum][j])
                chargeCurrents.append(channels[fileNum][i][j])
            rc('font', weight = 'bold')
            plt.tick_params('x', labelbottom = False)
            plt.tick_params(direction='in', labelsize = 'large', length = 5.0, width = 1.5, top = True, right = True)
            ax.spines['top'].set_linewidth(1.5)
            ax.spines['right'].set_linewidth(1.5)
            ax.spines['bottom'].set_linewidth(1.5)
            ax.spines['left'].set_linewidth(1.5)
            plt.xlim((xMin, xMax))
            if(i%8 == 0):
                plt.tick_params('x', labelbottom = True)
        else:
            for j in range(cycles[fileNum].index(2*cycle - 1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle)):
                if(cycles[fileNum][j] == 2*cycle - 1):
                    chargeVoltages.append(voltages[fileNum][j])
                    chargeCurrents.append(channels[fileNum][i][j])
                else:
                    dischargeVoltages.append(voltages[fileNum][j])
                    dischargeCurrents.append(channels[fileNum][i][j])
            plt.plot(chargeVoltages, chargeCurrents, 'b', linewidth = 1.75, label = 'Charge')
            plt.plot(dischargeVoltages, dischargeCurrents, 'r', linewidth = 1.75, label = 'Discharge')
            plt.tick_params('x', labelbottom = False)
            plt.tick_params(direction='in', labelsize = 'large', length = 5.0, width = 1.5, top = True, right = True)
            ax.spines['top'].set_linewidth(1.5)
            ax.spines['right'].set_linewidth(1.5)
            ax.spines['bottom'].set_linewidth(1.5)
            ax.spines['left'].set_linewidth(1.5)
            plt.xlim((xMin, xMax))
            if(i%8 == 0):
                plt.tick_params('x', labelbottom = True)
    plt.show()
    plt.close()

# Plots the current vs voltage for every cycle on the input channel
def plotAllCyclesCurrentVsVolts(channel, xMin = 3.1, xMax = 4.5, fileNum = 1):
    fileNum -= 1
    rc('font', weight = 'bold')
    plt.plot(voltages[fileNum], channels[fileNum][channel], 'b', linewidth = 2.0)
    plt.title('Current on Channel %d During each Cycle' % channel, fontsize = 'xx-large', weight = 'bold')
    plt.xlabel('Voltage (V)', fontsize = 'x-large', weight = 'bold')
    plt.ylabel('Current (\u03BCA)', fontsize = 'x-large', weight = 'bold')
    plt.tick_params(direction='in', labelsize = 'large', length = 5.0, width = 1.5, top = True, right = True)
    plt.xlim((xMin, xMax))
    plt.show()
    plt.close()

# Plots current vs voltage for every cycle on all 64 channels in one window
def plotAllCyclesCurrentVsVolts64(xMin = 3.1, xMax = 4.5, fileNum = 1):
    fileNum -= 1
    rc('font', weight = 'bold')
    plt.figure(1)
    for i in range(1,65):
        if(i <= 8):
            loc = i + (7*(i-1))
        elif(i <= 16):
            loc = i + (7*(i-10))
        elif(i <= 24):
            loc = i + (7*(i-19))
        elif(i <= 32):
            loc = i + (7*(i-28))
        elif(i <= 40):
            loc = i + (7*(i-37))
        elif(i <= 48):
            loc = i + (7*(i-46))
        elif(i <= 56):
            loc = i + (7*(i-55))
        else:
            loc = i + (7*(i-64))
        ax = plt.subplot(8, 8, loc, title = 'Channel %d' % i)
        ax.title.set_visible(False)
        plt.plot(voltages[fileNum], channels[fileNum][i], 'b', linewidth = 1.75)
        plt.tick_params('x', labelbottom = False)
        plt.tick_params(direction='in', labelsize = 'large', length = 5.0, width = 1.5, top = True, right = True)
        ax.spines['top'].set_linewidth(1.5)
        ax.spines['right'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        plt.xlim((xMin, xMax))
        if(i%8 == 0):
            plt.tick_params('x', labelbottom = True)
    plt.show()
    plt.close()
    
def plotNormalizedCurrentVsVolts(channel, cycle = 1, xMin = 3.1, xMax = 4.5, fileNum = 1):
    fileNum -= 1
    chargeVoltages = []
    chargeCurrents = []
    dischargeVoltages = []
    dischargeCurrents = []
    rc('font', weight = 'bold')
    plt.figure(1)
    ax = plt.subplot(111)
    if(mass[fileNum][channel] == 0.0):
        mass[fileNum][channel] = float(input('Enter the sample mass for channel %d in grams: ' % channel))
    
#   Loops through all the voltages for the specified cycle and stores them in the 
#   chargeVoltages list for the charge cycle and dischargeVoltages for the discharge 
#   cycle. Loops through the currents for the specified cycle stored in the channels 
#   dictionary at the specified channel and stores the charge currents in the 
#   chargeCurrents list and discharge currents in the dischargeCurrents list.
    if cycle*2 not in cycles[fileNum]:
        for i in range(cycles[fileNum].index(2*cycle - 1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle-1)):
            chargeVoltages.append(voltages[fileNum][i])
            chargeCurrents.append(channels[fileNum][channel][i]/mass[fileNum][channel])
        chargePlot, = plt.plot(chargeVoltages, chargeCurrents, 'b', linewidth = 2.0, label = 'Charge')
        dischargePlot, = plt.plot(dischargeVoltages, dischargeCurrents, 'r', linewidth = 2.0, label = 'Discharge')
        
        # Positions the legend to the top left corner outside the plot
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(handles = [chargePlot, dischargePlot], loc = 2, bbox_to_anchor = (1,1), fontsize = 'small')
        plt.title('Current on Channel %d During Cycle %d' % (channel, cycle))
        plt.xlabel('Voltage (V)')
        plt.ylabel('Current (\u03BCA/g)')
        plt.tick_params(direction='in')
        plt.xlim((xMin, xMax))
    else:
        for i in range(cycles[fileNum].index(2*cycle - 1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle)):
            if(cycles[fileNum][i] == 2*cycle - 1):
                chargeVoltages.append(voltages[fileNum][i])
                chargeCurrents.append(channels[fileNum][channel][i]/mass[fileNum][channel])
            else:
                dischargeVoltages.append(voltages[fileNum][i])
                dischargeCurrents.append(channels[fileNum][channel][i]/mass[fileNum][channel])
         
    # Plots the charge currents vs charge voltages in blue then the discharge currents vs 
    # discharge voltages in red
        chargePlot, = plt.plot(chargeVoltages, chargeCurrents, 'b', linewidth = 2.0, label = 'Charge')
        dischargePlot, = plt.plot(dischargeVoltages, dischargeCurrents, 'r', linewidth = 2.0, label = 'Discharge')
        
        # Positions the legend to the top left corner outside the plot
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(handles = [chargePlot, dischargePlot], loc = 2, bbox_to_anchor = (1,1), fontsize = 'small')
        plt.title('Current on Channel %d During Cycle %d' % (channel, cycle), fontsize = 'xx-large', weight = 'bold')
        plt.xlabel('Voltage (V)', fontsize = 'x-large', weight = 'bold')
        plt.ylabel('Current (\u03BCA/g)', fontsize = 'x-large', weight = 'bold')
        plt.tick_params(direction='in', labelsize = 'large', length = 5.0, width = 1.5, top = True, right = True)
        plt.xlim((xMin, xMax))
    plt.show()
    plt.close()    

def plotNormalizedCurrentVsVolts64(cycle = 1, xMin = 3.1, xMax = 4.5, fileNum = 1):
    fileNum -= 1
    rc('font', weight = 'bold')
    plt.figure(1)   

    # Loops through the 64 channels and adds each new plot as a subplot to figure 1,
    # using the same plotting logic as the plot function    
    for i in range(1,65):
        loc = 0
        if(i <= 8):
            loc = i + (7*(i-1))
        elif(i <= 16):
            loc = i + (7*(i-10))
        elif(i <= 24):
            loc = i + (7*(i-19))
        elif(i <= 32):
            loc = i + (7*(i-28))
        elif(i <= 40):
            loc = i + (7*(i-37))
        elif(i <= 48):
            loc = i + (7*(i-46))
        elif(i <= 56):
            loc = i + (7*(i-55))
        else:
            loc = i + (7*(i-64))
        ax = plt.subplot(8, 8, loc, title = 'Channel %d' % i)
        ax.title.set_visible(False)
        chargeVoltages = []
        chargeCurrents = []
        dischargeVoltages = []
        dischargeCurrents = []
        if cycle*2 not in cycles[fileNum]:
            for j in range(cycles[fileNum].index(2*cycle-1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle-1)):
                chargeVoltages.append(voltages[fileNum][j])
                chargeCurrents.append(channels[fileNum][i][j]/mass[fileNum][i])
            rc('font', weight = 'bold')
            plt.tick_params('x', labelbottom = False)
            plt.tick_params(direction='in', labelsize = 'large', length = 5.0, width = 1.5, top = True, right = True)
            ax.spines['top'].set_linewidth(1.5)
            ax.spines['right'].set_linewidth(1.5)
            ax.spines['bottom'].set_linewidth(1.5)
            ax.spines['left'].set_linewidth(1.5)
            plt.xlim((xMin, xMax))
            if(i%8 == 0):
                plt.tick_params('x', labelbottom = True)
        else:
            for j in range(cycles[fileNum].index(2*cycle - 1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle)):
                if(cycles[fileNum][j] == 2*cycle - 1):
                    chargeVoltages.append(voltages[fileNum][j])
                    chargeCurrents.append(channels[fileNum][i][j]/mass[fileNum][i])
                else:
                    dischargeVoltages.append(voltages[fileNum][j])
                    dischargeCurrents.append(channels[fileNum][i][j]/mass[fileNum][i])
            plt.plot(chargeVoltages, chargeCurrents, 'b', linewidth = 1.75, label = 'Charge')
            plt.plot(dischargeVoltages, dischargeCurrents, 'r', linewidth = 1.75, label = 'Discharge')
            plt.tick_params('x', labelbottom = False)
            plt.tick_params(direction='in', labelsize = 'large', length = 5.0, width = 1.5, top = True, right = True)
            ax.spines['top'].set_linewidth(1.5)
            ax.spines['right'].set_linewidth(1.5)
            ax.spines['bottom'].set_linewidth(1.5)
            ax.spines['left'].set_linewidth(1.5)
            plt.xlim((xMin, xMax))
            if(i%8 == 0):
                plt.tick_params('x', labelbottom = True)
    plt.show()
    plt.close()

def plotAllCyclesNormalizedCurrentVsVolts(channel, xMin = 3.1, xMax = 4.5, fileNum = 1):
    fileNum -= 1
    rc('font', weight = 'bold')
    channelCurrents = []
    for current in channels[fileNum][channel]:
        channelCurrents.append(current/mass[fileNum][channel])
    plt.plot(voltages[fileNum], channelCurrents, 'b', linewidth = 2.0)
    plt.title('Current on Channel %d During each Cycle' % channel, fontsize = 'xx-large', weight = 'bold')
    plt.xlabel('Voltage (V)', fontsize = 'x-large', weight = 'bold')
    plt.ylabel('Current (\u03BCA/g)', fontsize = 'x-large', weight = 'bold')
    plt.tick_params(direction='in', labelsize = 'large', length = 5.0, width = 1.5, top = True, right = True)
    plt.xlim((xMin, xMax))
    plt.show()
    plt.close()
    
def plotAllCyclesNormalizedCurrentVsVolts64(xMin = 3.1, xMax = 4.5, fileNum = 1):
    fileNum -= 1
    rc('font', weight = 'bold')
    plt.figure(1)
    for i in range(1,65):
        if(i <= 8):
            loc = i + (7*(i-1))
        elif(i <= 16):
            loc = i + (7*(i-10))
        elif(i <= 24):
            loc = i + (7*(i-19))
        elif(i <= 32):
            loc = i + (7*(i-28))
        elif(i <= 40):
            loc = i + (7*(i-37))
        elif(i <= 48):
            loc = i + (7*(i-46))
        elif(i <= 56):
            loc = i + (7*(i-55))
        else:
            loc = i + (7*(i-64))
        channelCurrents = []
        for current in channels[fileNum][i]:
            channelCurrents.append(current/mass[fileNum][i])
        ax = plt.subplot(8, 8, loc, title = 'Channel %d' % i)
        ax.title.set_visible(False)
        plt.plot(voltages[fileNum], channelCurrents, 'b', linewidth = 1.75)
        plt.tick_params('x', labelbottom = False)
        plt.tick_params(direction='in', labelsize = 'large', length = 5.0, width = 1.5, top = True, right = True)
        ax.spines['top'].set_linewidth(1.5)
        ax.spines['right'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        plt.xlim((xMin, xMax))
        if(i%8 == 0):
            plt.tick_params('x', labelbottom = True)
    plt.show()
    plt.close()
    
# Plots voltage as a function of current for the input channel during the input cycle
def plotVoltsVsCurrent(channel, cycle = 1, fileNum = 1):
    fileNum -= 1
    chargeCurrents = []
    chargeVoltages = []
    dischargeCurrents = []
    dischargeVoltages = []
    plt.figure(1)
    ax = plt.subplot(111)
    if 2*cycle in cycles[fileNum]:
        for i in range(cycles[fileNum].index(2*cycle - 1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle)):
            if(cycles[fileNum][i] == 2*cycle - 1): 
                chargeCurrents.append(channels[fileNum][channel][i])
                chargeVoltages.append(voltages[fileNum][i])
            else: 
                dischargeCurrents.append(channels[fileNum][channel][i])
                dischargeVoltages.append(voltages[fileNum][i])
        chargePlot, = plt.plot(chargeCurrents, chargeVoltages, 'b', linewidth = 2, label = 'Charge')
        dischargePlot, = plt.plot(dischargeCurrents, dischargeVoltages, 'r', linewidth = 2, label = 'Discharge')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(handles = [chargePlot, dischargePlot], loc = 2, bbox_to_anchor = (1,1), fontsize = 'small')
        plt.title('Voltage on Channel %d During Cycle %d' % (channel, cycle))
        plt.xlabel('Current (\u03BCA)')
        plt.ylabel('Voltage (V)')
    else:
        for i in range(cycles[fileNum].index(2*cycle - 1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle-1)):
            chargeCurrents.append(channels[fileNum][channel][i])
            chargeVoltages.append(voltages[fileNum][i])
        chargePlot, = plt.plot(chargeCurrents, chargeVoltages, 'b', linewidth = 2, label = 'Charge')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(handles = [chargePlot, dischargePlot], loc = 2, bbox_to_anchor = (1,1), fontsize = 'small')
        plt.title('Voltage on Channel %d During Cycle %d' % (channel, cycle))
        plt.xlabel('Current (\u03BCA)')
        plt.ylabel('Voltage (V)')
    plt.tick_params(direction='in', labelsize = 'large', length = 5.0)
    plt.show()
    plt.close()

# Plots voltage vs current for all 64 channels in one window for the imput cycle
def plotVoltsVsCurrent64(cycle = 1, fileNum = 1):
    fileNum -= 1
    plt.figure(1)
    for i in range(1,65):
        if(i <= 8):
            loc = i + (7*(i-1))
        elif(i <= 16):
            loc = i + (7*(i-10))
        elif(i <= 24):
            loc = i + (7*(i-19))
        elif(i <= 32):
            loc = i + (7*(i-28))
        elif(i <= 40):
            loc = i + (7*(i-37))
        elif(i <= 48):
            loc = i + (7*(i-46))
        elif(i <= 56):
            loc = i + (7*(i-55))
        else:
            loc = i + (7*(i-64))
        chargeCurrents = []
        chargeVoltages = []
        dischargeCurrents = []
        dischargeVoltages = []
        if 2*cycle in cycles[fileNum]:
            for j in range(cycles[fileNum].index(2*cycle - 1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle)):
                if(cycles[fileNum][j] == 2*cycle - 1): 
                    chargeCurrents.append(channels[fileNum][i][j])
                    chargeVoltages.append(voltages[fileNum][j])
                else: 
                    dischargeCurrents.append(channels[fileNum][i][j])
                    dischargeVoltages.append(voltages[fileNum][j])
            
            ax = plt.subplot(8, 8, loc, title = 'Channel %d' % i)
            ax.title.set_visible(False)
            plt.plot(chargeCurrents, chargeVoltages, 'b', linewidth = 1.75, label = 'Charge')
            plt.plot(dischargeCurrents, dischargeVoltages, 'r', linewidth = 1.75, label = 'Discharge')
            if(i%8 == 0):
                plt.xlabel('Current (\u03BCA)')
            if(i <= 8):
                plt.ylabel('Voltage (V)')
        else:
            for j in range(cycles[fileNum].index(2*cycle - 1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle-1)):
                chargeCurrents.append(channels[fileNum][i][j])
                chargeVoltages.append(voltages[fileNum][j])
            ax = plt.subplot(8, 8, loc, title = 'Channel %d' % i)
            ax.title.set_visible(False)
            plt.plot(chargeCurrents, chargeVoltages, 'b', linewidth = 1.75, label = 'Charge')
            if(i%8 == 0):
                plt.xlabel('Current (\u03BCA)')
            if(i <= 8):
                plt.ylabel('Voltage (V)')
    plt.tick_params(direction='in', labelsize = 'large', length = 5.0)
    plt.show()
    plt.close()

# Plots the voltage vs current for the input channel
def plotAllCyclesVoltsVsCurrent(channel, fileNum = 1):
    fileNum -= 1
    plt.plot(channels[fileNum][channel], voltages[fileNum], 'b', linewidth = 2.0)
    plt.title('Voltage on Channel %d During every Cycle' % channel)
    plt.xlabel('Current (\u03BCA)')
    plt.ylabel('Voltage (V)')
    plt.tick_params(direction='in', labelsize = 'large', length = 5.0)
    plt.show()
    plt.close()
    
# Plots voltage vs current for all cycles on all 64 channels in one window
def plotAllCyclesVoltsVsCurrent64(fileNum = 1):
    fileNum -= 1
    plt.figure(1)
    for i in range(1,65):
        if(i <= 8):
            loc = i + (7*(i-1))
        elif(i <= 16):
            loc = i + (7*(i-10))
        elif(i <= 24):
            loc = i + (7*(i-19))
        elif(i <= 32):
            loc = i + (7*(i-28))
        elif(i <= 40):
            loc = i + (7*(i-37))
        elif(i <= 48):
            loc = i + (7*(i-46))
        elif(i <= 56):
            loc = i + (7*(i-55))
        else:
            loc = i + (7*(i-64))
        ax = plt.subplot(8, 8, loc, title = 'Channel %d' % i)
        ax.title.set_visible(False)
        plt.plot(channels[fileNum][i], voltages[fileNum], 'b', linewidth = 1.75)
        plt.tick_params(direction='in', labelsize = 'large', length = 5.0)
        if(i%8 == 0):
            plt.xlabel('Current (\u03BCA)')
        if(i <= 8):
            plt.ylabel('Voltage (V)')
    plt.show()
    plt.close()

# Plots voltage vs time during the input cycle
def plotVoltsVsTime(cycle = 1, fileNum = 1):
    fileNum -= 1
    cycleVoltages = []
    cycleTimes = []
    
    # Stores the voltages and times during the input cycle in the cycleVoltages and
    # cycleTimes lists respectively
    if cycle*2 in cycles[fileNum]:
        for i in range(cycles[fileNum].index(2*cycle-1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle)):
            cycleVoltages.append(voltages[fileNum][i])
            cycleTimes.append(times[fileNum][i])
        plt.plot(cycleTimes, cycleVoltages, linewidth = 2.0)
        plt.title('Change in Voltage During Cycle %d' % cycle)
        plt.xlabel('Time (hr)')
        plt.ylabel('Voltage (V)')
    else:
        for i in range(cycles[fileNum].index(2*cycle-1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle-1)):
            cycleVoltages.append(voltages[fileNum][i])
            cycleTimes.append(times[fileNum][i])
        plt.plot(cycleTimes, cycleVoltages, linewidth = 2.0)
        plt.title('Change in Voltage During Cycle %d' % cycle)
        plt.xlabel('Time (hr)')
        plt.ylabel('Voltage (V)')
    plt.tick_params(direction='in', labelsize = 'large', length = 5.0)
    plt.show()
    plt.close()
    
# Plots voltage vs time for every cycle
def plotAllCyclesVoltsVsTime(fileNum = 1):
    fileNum -= 1
    plt.plot(times[fileNum], voltages[fileNum], linewidth = 2.0)
    plt.title('Change in Voltage During each Cycle')
    plt.xlabel('Time (hr)')
    plt.ylabel('Voltage (V)')
    plt.show()
    plt.close()
    
# Plots current vs time for the input channel and cycle
def plotCurrentVsTime(channel, cycle = 1, fileNum = 1):
    fileNum -= 1
    cycleCurrents = []
    cycleTimes = []
    
    # Stores the currents and times for the input cycle in cycleCurrents and cycleTimes
    if cycle*2 in cycles[fileNum]:
        for i in range(cycles[fileNum].index(2*cycle-1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle)):
            cycleCurrents.append(channels[fileNum][channel][i])
            cycleTimes.append(times[fileNum][i])
        plt.plot(cycleTimes, cycleCurrents, linewidth = 2.0)
        plt.title('Change in Current on Channel %d During Cycle %d' % (channel, cycle))
        plt.xlabel('Time (hr)')
        plt.ylabel('Current (\u03BCA)')
    else:
        for i in range(cycles[fileNum].index(2*cycle-1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle-1)):
            cycleCurrents.append(channels[fileNum][channel][i])
            cycleTimes.append(times[fileNum][i])
        plt.plot(cycleTimes, cycleCurrents, linewidth = 2.0)
        plt.title('Change in Current on Channel %d During Cycle %d' % (channel, cycle))
        plt.xlabel('Time (hr)')
        plt.ylabel('Current (\u03BCA)')
    plt.show()
    plt.close()

# Plots current vs time for all the cycles on the input channel
def plotAllCyclesCurrentVsTime(channel,fileNum = 1):
    fileNum -= 1
    plt.plot(times[fileNum], channels[fileNum][channel], linewidth = 2.0)
    plt.title('Change in Current on Channel %d During each Cycle' % channel)
    plt.xlabel('Time (hr)')
    plt.ylabel('Current (\u03BCA)')
    plt.show()
    plt.close()

# Plots the voltage vs capacity for the input channel during the input cycle
def plotVoltsVsCapacity(channel, cycle = 1, fileNum = 1):
    fileNum -= 1
    chargeCapacities = []
    chargeVoltages = []
    dischargeCapacities = []
    dischargeVoltages = []
    plt.figure(1)
    ax = plt.subplot(111)
    
    if(mass[fileNum][channel] == 0.0):
        mass[fileNum][channel] = float(input('Enter the sample mass for channel %d in grams: ' % channel))
    
    # Calculates specific capacity for the charge and discharge cycles and appends those
    # values to the cycleCapacities list and appends voltages during those cycles to the
    # cycleVoltages list
    if 2*cycle in cycles[fileNum]:
        for i in range(cycles[fileNum].index(2*cycle-1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle)):
            if(cycles[fileNum][i] == 2*cycle - 1):
                if(i != len(cycles[fileNum]) - 1):
                    chargeCapacities.append(((0.5*(channels[fileNum][channel][i] + channels[fileNum][channel][i+1])*(times[fileNum][i+1] - times[fileNum][i]))/mass[fileNum][channel])/1000.0)
                chargeVoltages.append(voltages[fileNum][i])
            else:
                if(i != len(cycles[fileNum]) - 1):
                    dischargeCapacities.append(((0.5*(channels[fileNum][channel][i] + channels[fileNum][channel][i+1])*(times[fileNum][i+1] - times[fileNum][i]))/mass[fileNum][channel])/1000.0)
                dischargeVoltages.append(voltages[fileNum][i])
        chargePlot, = ax.plot(chargeCapacities, chargeVoltages, 'b', linewidth = 2.0, label = 'Charge')
        dischargePlot, = ax.plot(dischargeCapacities, dischargeVoltages, 'r', linewidth = 2.0, label = 'Discharge')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(handles = [chargePlot, dischargePlot], loc = 2, bbox_to_anchor = (1,1), fontsize = 'small')
        plt.title('Voltage as a Function of Specific Capacity on Channel %d' % channel)
        plt.xlabel('Specific Capacity (mAh/g)')
        plt.ylabel('Voltage (V)')
    else:
        for i in range(cycles[fileNum].index(2*cycle-1), len(cycles[fileNum]) - cycles[fileNum][::-1].index(2*cycle-1)):
            if(i != len(cycles[fileNum]) - 1):
                chargeCapacities.append(((0.5*(channels[fileNum][channel][i] + channels[fileNum][channel][i+1])*(times[fileNum][i+1] - times[fileNum][i]))/mass[fileNum][channel])/1000.0)
            chargeVoltages.append(voltages[fileNum][i])
        chargePlot, = ax.plot(chargeCapacities, chargeVoltages, 'b', linewidth = 2.0, label = 'Charge')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(handles = [chargePlot, dischargePlot], loc = 2, bbox_to_anchor = (1,1), fontsize = 'small')
        plt.title('Voltage as a Function of Specific Capacity on Channel %d' % channel)
        plt.xlabel('Specific Capacity (mAh/g)')
        plt.ylabel('Voltage (V)')
    plt.show()
    plt.close()

# Plots the charge and discharge capacities of the input channel for each cycle
def plotCapacityVsCycle(channel, fileNum = 1):
    fileNum -= 1
    chgCapacity = []
    dchgCapacity = []
    cycle = []
    plt.figure(1)
    ax = plt.subplot(111)
    
    # Adds each cycle number to the cycle list
    for i in range(1, int(max(cycles[fileNum])/2) + 1):
        cycle.append(i)
        
    # Adds the capacity for each cycle to the cycleCapacity list
    for j in cycle:
        charge = capacity(channel, j, fileNum+1)['Charge capacity (mAh/g)']
        discharge = capacity(channel, j, fileNum+1)['Discharge capacity (mAh/g)']
        chgCapacity.append(charge)
        dchgCapacity.append(-1.0*discharge)
    chgplot, = ax.plot(cycle, chgCapacity, 'bo', label = 'Charge', ms = 4.0)
    dchgplot, = ax.plot(cycle, dchgCapacity, 'ro', label = 'Discharge', ms = 4.0)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    ax.legend(handles = [chgplot, dchgplot], loc = 2, bbox_to_anchor = (1,1), fontsize = 'small')
    plt.title('Capacity at each Cycle for Channel %d' % channel)
    if(len(cycle) > 10):
        ticks = []
        for x in cycle[4::5]:
            ticks.append(x)
        plt.xticks(ticks)
    else:
        plt.xticks([x for x in cycle])
    plt.xlabel('Cycle Number')
    plt.ylabel('Specific Capacity (mAh/g)')
    plt.show() 
    plt.close()   

# Plots the charge capacity of the input channel for each cycle
def plotChargeCapacityVsCycle(channel, fileNum = 1):
    fileNum -= 1
    cycleCapacity = []
    cycle = []
    
    # Adds each cycle number to the cycle list
    for i in range(1, int(max(cycles[fileNum])/2) + 1):
        cycle.append(i)
        
    # Adds the capacity for each cycle to the cycleCapacity list
    for j in cycle:
        charge = capacity(channel, j, fileNum+1)['Charge capacity (mAh/g)']
        cycleCapacity.append(charge)
    plt.plot(cycle, cycleCapacity, 'bo', ms = 4.0)
#     for value in cycle:
#         plt.annotate('%f' % cycleCapacity[value-1], xy = (value, cycleCapacity[value-1]), xytext = (value, cycleCapacity[value-1] + 0.007*cycleCapacity[value-1]))
    plt.title('Charge Capacity at each Cycle for Channel %d' % channel)
    if(len(cycle) > 10):
        ticks = []
        for x in cycle[4::5]:
            ticks.append(x)
        plt.xticks(ticks)
    else:
        plt.xticks([x for x in cycle])
    plt.xlabel('Cycle Number')
    plt.ylabel('Charge Capacity (mAh/g)')
    plt.show()
    plt.close()

# Plots the discharge capacity of the input channel for each cycle
def plotDischargeCapacityVsCycle(channel, fileNum = 1):
    fileNum -= 1
    cycleCapacity = []
    cycle = []
    
    # Adds each cycle number to the cycle list
    for i in range(1, int(max(cycles[fileNum])/2) + 1):
        cycle.append(i)
        
    # Adds the capacity for each cycle to the cycleCapacity list
    for j in cycle:
        discharge = capacity(channel, j, fileNum+1)['Discharge capacity (mAh/g)']
        cycleCapacity.append(-1.0*discharge)
    plt.plot(cycle, cycleCapacity, 'ro', ms = 4.0)
#     for value in cycle:
#         plt.annotate('%f' % cycleCapacity[value-1], xy = (value, cycleCapacity[value-1]), xytext = (value, cycleCapacity[value-1] + 0.007*cycleCapacity[value-1]))
    plt.title('Discharge Capacity at each Cycle for Channel %d' % channel)
    if(len(cycle) > 10):
        ticks = []
        for x in cycle[4::5]:
            ticks.append(x)
        plt.xticks(ticks)
    else:
        plt.xticks([x for x in cycle])
    plt.xlabel('Cycle Number')
    plt.ylabel('Discharge Capacity (mAh/g)')
    plt.show()
    plt.close()
    
print('\nAvailable functions:')
print('newFile()')
print('setMass(channel)')
print('export()')
print('capacity(channel, cycle = 1, fileNum = 1)')
print('capacity64(cycle = 1, fileNum = 1)')
print('capacityAllCycles(channel)')
print('allCapacities()')
print('averageVoltages(startTime, endTime)')
print('resistances(cycle = 1, fileNum = 1)')
print('allResistances()')
print('plotCurrentVsVolts(channel, cycle = 1, fileNum = 1)')
print('plotCurrentVsVolts64(cycle)')
print('plotAllCyclesCurrentVsVolts(channel)')
print('plotAllCyclesCurrentVsVolts64()')
print('plotVoltsVsCurrent(channel, cycle = 1, fileNum = 1)')
print('plotVoltsVsCurrent64(cycle = 1, fileNum = 1)')
print('plotAllCyclesVoltsVsCurrent(channel)')
print('plotAllCyclesVoltsVsCurrent64(channel)')
print('plotVoltsVsTime(cycle = 1, fileNum = 1)')
print('plotAllCyclesVoltsVsTime()')
print('plotCurrentVsTime(channel, cycle = 1, fileNum = 1)')
print('plotAllCyclesCurrentVsTime(channel)')
print('plotVoltsVsCapacity(channel, cycle = 1, fileNum = 1)')
print('plotCapacityVsCycle(channel)')
print('plotChargeCapacityVsCycle(channel)')
print('plotDischargeCapacityVsCycle(channel)\n')

newFile()
# allCapacities()
# allResistances()
# setMass(1)
# export()
# capacity(1,1)
# capacityAllCycles(1)
# setMass(1)
# capacity64(1)
# averageVoltages(0, 13.01)
# resistances(1)
# plotVoltsVsCurrent(59, 1)
# plotVoltsVsCurrent64(1)
# plotCurrentVsVolts(1)
# plotCurrentVsVolts64()
# plotAllCyclesCurrentVsVolts(1)
# plotAllCyclesCurrentVsVolts64()
# plotNormalizedCurrentVsVolts(1)
# plotNormalizedCurrentVsVolts64()
# plotAllCyclesNormalizedCurrentVsVolts(1)
# plotAllCyclesNormalizedCurrentVsVolts64()
# plotAllCyclesVoltsVsCurrent64()
# plotVoltsVsTime(1)
# plotAllCyclesVoltsVsTime()
# plotChargeCapacityVsCycle(1)
# plotDischargeCapacityVsCycle(1)
# plotAllCyclesCurrentVsTime(1)
# plotCurrentVsTime(1, 1)
# plotVoltsVsCapacity(59, 1)
# plotCapacityVsCycle(1)