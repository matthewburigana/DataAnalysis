import csv
from itertools import islice
import matplotlib.pyplot as plt
import numpy as np

filepath = 'C:\\Users\\McCalla Lab\\Documents\\Data\\Matthew\\Cells\\'

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
mass = {x: 0.0 for x in range(1, 65)}

# A dictionary to store the currents from each channel where the channel
# number is the key and the value is a list of currents
channels = {x: [] for x in range(1, 65)}

# Initializes a new file by clearing all previously stored data then processing the data stored
# by the new input file
def newFile():
    filename = input('Enter the CSV file name: ') + '.csv'
    massFile = input('Enter the mass file name: ') + '.csv'
    if (cycles != []):
        cycles.clear()
        times.clear()
        voltages.clear()
        for key in channels:
            mass[key] = 0.0
            channels[key].clear()
        chgCapacities.clear()
        dchgCapacities.clear()
        resistance.clear()
        meanResistances.clear()
        stdevResistances.clear()
        
    with open(filepath + filename, 'r') as inputData:
        reader = csv.reader(inputData, dialect = 'excel')
        
        # Appends the cycle numbers, times, and average voltages to their respective
        # list for each row in the CSV file
        for row in islice(reader, 8, None):
            if(row[0] == ''):
                break
            cycles.append(int(row[0]))
            times.append(float(row[1]))
            voltages.append(float(row[5]))
        
            # Appends the current for each channel in the current list at each channel
            # in the channels dictionary
            for j in range(1, 65):
                channels[j].append(float(row[j+5]))

    # Appends empty lists to store resistances and capacities when they are calculated
    for k in range(int(max(cycles)/2)):
        resistance.append([])
        chgCapacities.append([])
        dchgCapacities.append([])
        
    # Reads masses stored in the same file format as an exported file
    if(massFile != ''):
        with open(filepath + massFile, 'r') as masses:
            massRead = csv.reader(masses, dialect = 'excel')
            for row in islice(massRead, 1, None):
                for i in range(1,65):
                    mass[i] = float(row[i+1])
                break
    
# Set the mass of the sample on a channel
def setMass(channel):
    newMass = float(input('Enter the new sample mass for channel %d in grams: ' % channel))
    for i in range(min(cycles), int(max(cycles)/2) + 1):
        chgCapacities[i-1][channel-1] = chgCapacities[i-1][channel-1]*mass[channel]/newMass
        dchgCapacities[i-1][channel-1] = dchgCapacities[i-1][channel-1]*mass[channel]/newMass
    mass[channel] = newMass
    
# Exports the masses, charge and discharge capacities and resistances to a CSV file
def export():
    outputFile = input('Enter the export file name: ') + '.csv'
    allCapacities()
    allResistances()
    with open(filepath + outputFile, 'w', newline = '') as output:
        writer = csv.writer(output, dialect = 'excel')
        channel = []
        channelMass = []
        for i in range(1,65):
            channel.append('Channel %d' % i)
            channelMass.append(mass[i])
        writer.writerow(['Cycle', 'Name'] + channel)
        writer.writerow(['','Mass (g)'] + channelMass)
        for cycle in range(min(cycles), int(max(cycles)/2) + 1):
            writer.writerow([cycle, 'Charge (mAh/g)'] + chgCapacities[cycle-1])
            writer.writerow(['', 'Discharge (mAh/g)'] + dchgCapacities[cycle-1])
            writer.writerow(['', 'Resistance'] + resistance[cycle-1])
            writer.writerow(['', 'Mean Resistance'] + [meanResistances[cycle-1]])
            writer.writerow(['', 'Standard Deviation'] + [stdevResistances[cycle-1]])
        
# Calculates the charge and discharge capacity of a channel during the specified cycle
# between the minimum and maximum voltages. Integrates dQ = Idt
def capacity(channel, cycle):
    qCharge = 0.0
    qDischarge = 0.0
    if(mass[channel] == 0.0):
        mass[channel] = float(input('Enter the sample mass for channel %d in grams: ' % channel))
    for i in range(cycles.index(2*cycle - 1), len(cycles) - cycles[::-1].index(2*cycle)):
        
        # Only adds to charge capacity during a charging cycle and only if the voltages are 
        # between those specified. The average of the capacity and the next capacity are taken
        # and multiplied by the difference in time of those measurements
        if(2*cycle - 1 == cycles[i] and i != len(cycles) - 1):
            qCharge += 0.5*(channels[channel][i] + channels[channel][i+1])*(times[i+1] - times[i])
            
        # Adds to discharge capacity during a discharge cycle if the voltages are within the 
        # specified voltages and the iterator is not at the last value of the cycles
        elif(2*cycle == cycles[i] and i != len(cycles) - 1):
            qDischarge += 0.5*(channels[channel][i] + channels[channel][i+1])*(times[i+1] - times[i])
                
    # Converts the capacities to Ah/g then returns a dictionary listing the channel, cycle,
    # charge and discharge capacities
    qCharge = round(qCharge/mass[channel]/1000.0, 8)
    qDischarge = round(qDischarge/mass[channel]/1000.0, 8)
    if(len(chgCapacities[cycle-1]) < 64):
        chgCapacities[cycle-1].append(qCharge)
        dchgCapacities[cycle-1].append(qDischarge)
    return({'Channel': channel, 'Cycle': cycle, 'Charge capacity (mAh/g)': qCharge, 'Discharge capacity (mAh/g)': qDischarge})

# Prints the capacity of all 64 channels in the input cycle
def capacity64(cycle):
    for i in range(1,65):
        print(capacity(i, cycle))
      
# Prints the capacity of each cycle for the input channel
def capacityAllCycles(channel):
    for i in range(min(cycles), int(max(cycles)/2) + 1):
        print(capacity(channel, i))
    
# Prints the capacities of all the channels for all the cycles    
def allCapacities():
    for i in range(1,65):
        for j in range(min(cycles), int(max(cycles)/2) + 1):
            print(capacity(i, j))

# Calculates then prints the average voltages of each channel using Vavg = VdQ/Q
def averageVoltages(startTime, endTime):
    capacities = {}
    for i in range(1, 65):
        capacities[i] = []
        voltage = 0.0
        
        # Appends the specific capacity at each time to the capacities dictionary then calculates the 
        # differential voltage and adds that to the sum of VdQ then divides that sum by the total capacity
        # to give the average voltage
        for j in range(len(times)):
            if(times[j] >= startTime and times[j] <= endTime and j < len(times) - 1):
                capacities[i].append(0.5*(channels[i][j] + channels[i][j+1])*(times[j+1] - times[j])/1000000.0)
        totalCapacity = sum(capacities[i])
        for k in range(len(times)):
            if(times[k] >= startTime and times[k] <= endTime and k < len(capacities[i]) - 1):
                voltage += 0.5*(voltages[k] + voltages[k+1])*(capacities[i][k+1] - capacities[i][k])/totalCapacity
        print('Channel %d average voltage: %f V' % (i, voltage))
    
# Prints the resistances of all 64 channels based on the linear regression of 
# voltage vs current for the charging cycle
def resistances(cycle):
    for i in range(1, 65):
        cycleCurrents = []
        cycleVoltages = []
        for j in range(cycles.index(2*cycle-1), len(cycles) - cycles[::-1].index(2*cycle-1)):
            cycleCurrents.append(channels[i][j]/1000000.0)
            cycleVoltages.append(voltages[j])
        param = np.polyfit(cycleCurrents, cycleVoltages, 1)
        print('Channel %d resistance: %f Ohms' % (i, param[0]))
        if(len(resistance[cycle-1]) < 64):
            resistance[cycle-1].append(param[0])
        
    # Calculates and prints the mean, standard deviation and variance of the 64 resistance
    # values
    stdev = np.std(resistance[cycle-1], ddof = 1)
    variance = np.var(resistance[cycle-1], ddof = 1)
    mean = np.mean(resistance[cycle-1])
    if(len(meanResistances) < max(cycles)/2):
        meanResistances.append(mean)
        stdevResistances.append(stdev)
    print('Average: %f' % mean)
    print('Standard Deviation: %f' % stdev)
    print('Variance: %f' % variance)

# Prints resistances for ever channel during each cycle
def allResistances():
    for i in range(min(cycles), int(max(cycles)/2) + 1):
        print('\nCycle %d' % i)
        resistances(i)
    
# Plots the current vs voltage for one channel during the input cycle
def plotCurrentVsVolts(channel, cycle):
    chargeVoltages = []
    chargeCurrents = []
    dischargeVoltages = []
    dischargeCurrents = []
    plt.figure(1)
    ax = plt.subplot(111)
    
#   Loops through all the voltages for the specified cycle and stores them in the 
#   chargeVoltages list for the charge cycle and dischargeVoltages for the discharge 
#   cycle. Loops through the currents for the specified cycle stored in the channels 
#   dictionary at the specified channel and stores the charge currents in the 
#   chargeCurrents list and discharge currents in the dischargeCurrents list.
    if cycle*2 not in cycles:
        for i in range(cycles.index(2*cycle - 1), len(cycles) - cycles[::-1].index(2*cycle-1)):
            chargeVoltages.append(voltages[i])
            chargeCurrents.append(channels[channel][i])
        chargePlot, = plt.plot(chargeVoltages, chargeCurrents, 'b', linewidth = 0.75, label = 'Charge')
        dischargePlot, = plt.plot(dischargeVoltages, dischargeCurrents, 'r', linewidth = 0.75, label = 'Discharge')
        
        # Positions the legend to the top left corner outside the plot
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(handles = [chargePlot, dischargePlot], loc = 2, bbox_to_anchor = (1,1), fontsize = 'small')
        plt.title('Current on Channel %d During Cycle %d' % (channel, cycle))
        plt.xlabel('Voltage (V)')
        plt.ylabel('Current (\u03BCA)')
    else:
        for i in range(cycles.index(2*cycle - 1), len(cycles) - cycles[::-1].index(2*cycle)):
            if(cycles[i] == 2*cycle - 1):
                chargeVoltages.append(voltages[i])
                chargeCurrents.append(channels[channel][i])
            else:
                dischargeVoltages.append(voltages[i])
                dischargeCurrents.append(channels[channel][i])
         
    # Plots the charge currents vs charge voltages in blue then the discharge currents vs 
    # discharge voltages in red
        chargePlot, = plt.plot(chargeVoltages, chargeCurrents, 'b', linewidth = 0.75, label = 'Charge')
        dischargePlot, = plt.plot(dischargeVoltages, dischargeCurrents, 'r', linewidth = 0.75, label = 'Discharge')
        
        # Positions the legend to the top left corner outside the plot
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(handles = [chargePlot, dischargePlot], loc = 2, bbox_to_anchor = (1,1), fontsize = 'small')
        plt.title('Current on Channel %d During Cycle %d' % (channel, cycle))
        plt.xlabel('Voltage (V)')
        plt.ylabel('Current (\u03BCA)')
    plt.tick_params(direction='in', labelsize = 'large', length = 5.0)
    plt.show()
    plt.close()

# Plots the current vs voltage for all 64 channels in one window
def plotCurrentVsVolts64(cycle):
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
        if cycle*2 not in cycles:
            for j in range(cycles.index(2*cycle-1), len(cycles) - cycles[::-1].index(2*cycle-1)):
                chargeVoltages.append(voltages[j])
                chargeCurrents.append(channels[i][j])
            plt.plot(chargeVoltages, chargeCurrents, 'b', linewidth = 0.75, label = 'Charge')
            if(i%8 == 0):
                plt.xlabel('Voltage (V)')
            if(i <= 8):
                plt.ylabel('Current (\u03BCA)')
        else:
            for j in range(cycles.index(2*cycle - 1), len(cycles) - cycles[::-1].index(2*cycle)):
                if(cycles[j] == 2*cycle - 1):
                    chargeVoltages.append(voltages[j])
                    chargeCurrents.append(channels[i][j])
                else:
                    dischargeVoltages.append(voltages[j])
                    dischargeCurrents.append(channels[i][j])
            plt.plot(chargeVoltages, chargeCurrents, 'b', linewidth = 0.75, label = 'Charge')
            plt.plot(dischargeVoltages, dischargeCurrents, 'r', linewidth = 0.75, label = 'Discharge')
            if(i%8 == 0):
                plt.xlabel('Voltage (V)')
            if(i <= 8):
                plt.ylabel('Current (\u03BCA)')
    plt.tick_params(direction='in', labelsize = 'large', length = 5.0)
    plt.show()
    plt.close()

# Plots the current vs voltage for every cycle on the input channel
def plotAllCyclesCurrentVsVolts(channel):
    plt.plot(voltages, channels[channel], linewidth = 1.0)
    plt.title('Current on Channel %d During each Cycle' % channel)
    plt.xlabel('Voltage (V)')
    plt.ylabel('Current (\u03BCA)')
    plt.show()
    plt.close()

# Plots current vs voltage for every cycle on all 64 channels in one window
def plotAllCyclesCurrentVsVolts64():
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
        plt.plot(voltages, channels[i], linewidth = 0.75)
        if(i%8 == 0):
            plt.xlabel('Voltage (V)')
        if(i <= 8):
            plt.ylabel('Current (\u03BCA)')
    plt.tick_params(direction='in', labelsize = 'large', length = 5.0)
    plt.show()
    plt.close()
    
# Plots voltage as a function of current for the input channel during the input cycle
def plotVoltsVsCurrent(channel, cycle):
    chargeCurrents = []
    chargeVoltages = []
    dischargeCurrents = []
    dischargeVoltages = []
    plt.figure(1)
    ax = plt.subplot(111)
    if 2*cycle in cycles:
        for i in range(cycles.index(2*cycle - 1), len(cycles) - cycles[::-1].index(2*cycle)):
            if(cycles[i] == 2*cycle - 1): 
                chargeCurrents.append(channels[channel][i])
                chargeVoltages.append(voltages[i])
            else: 
                dischargeCurrents.append(channels[channel][i])
                dischargeVoltages.append(voltages[i])
        chargePlot, = plt.plot(chargeCurrents, chargeVoltages, 'b', linewidth = 0.75, label = 'Charge')
        dischargePlot, = plt.plot(dischargeCurrents, dischargeVoltages, 'r', linewidth = 0.75, label = 'Discharge')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(handles = [chargePlot, dischargePlot], loc = 2, bbox_to_anchor = (1,1), fontsize = 'small')
        plt.title('Voltage on Channel %d During Cycle %d' % (channel, cycle))
        plt.xlabel('Current (\u03BCA)')
        plt.ylabel('Voltage (V)')
    else:
        for i in range(cycles.index(2*cycle - 1), len(cycles) - cycles[::-1].index(2*cycle-1)):
            chargeCurrents.append(channels[channel][i])
            chargeVoltages.append(voltages[i])
        chargePlot, = plt.plot(chargeCurrents, chargeVoltages, 'b', linewidth = 0.75, label = 'Charge')
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
def plotVoltsVsCurrent64(cycle):
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
        if 2*cycle in cycles:
            for j in range(cycles.index(2*cycle - 1), len(cycles) - cycles[::-1].index(2*cycle)):
                if(cycles[j] == 2*cycle - 1): 
                    chargeCurrents.append(channels[i][j])
                    chargeVoltages.append(voltages[j])
                else: 
                    dischargeCurrents.append(channels[i][j])
                    dischargeVoltages.append(voltages[j])
            
            ax = plt.subplot(8, 8, loc, title = 'Channel %d' % i)
            ax.title.set_visible(False)
            plt.plot(chargeCurrents, chargeVoltages, 'b', linewidth = 0.75, label = 'Charge')
            plt.plot(dischargeCurrents, dischargeVoltages, 'r', linewidth = 0.75, label = 'Discharge')
            if(i%8 == 0):
                plt.xlabel('Current (\u03BCA)')
            if(i <= 8):
                plt.ylabel('Voltage (V)')
        else:
            for j in range(cycles.index(2*cycle - 1), len(cycles) - cycles[::-1].index(2*cycle-1)):
                chargeCurrents.append(channels[i][j])
                chargeVoltages.append(voltages[j])
            ax = plt.subplot(8, 8, loc, title = 'Channel %d' % i)
            ax.title.set_visible(False)
            plt.plot(chargeCurrents, chargeVoltages, 'b', linewidth = 0.75, label = 'Charge')
            if(i%8 == 0):
                plt.xlabel('Current (\u03BCA)')
            if(i <= 8):
                plt.ylabel('Voltage (V)')
    plt.tick_params(direction='in', labelsize = 'large', length = 5.0)
    plt.show()
    plt.close()

# Plots the voltage vs current for the input channel
def plotAllCyclesVoltsVsCurrent(channel):
    plt.plot(channels[channel], voltages, linewidth = 0.75)
    plt.title('Voltage on Channel %d During every Cycle' % channel)
    plt.xlabel('Current (\u03BCA)')
    plt.ylabel('Voltage (V)')
    plt.tick_params(direction='in', labelsize = 'large', length = 5.0)
    plt.show()
    plt.close()
    
# Plots voltage vs current for all cycles on all 64 channels in one window
def plotAllCyclesVoltsVsCurrent64():
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
        plt.plot(channels[i], voltages, linewidth = 0.75)
        if(i%8 == 0):
            plt.xlabel('Current (\u03BCA)')
        if(i <= 8):
            plt.ylabel('Voltage (V)')
    plt.tick_params(direction='in', labelsize = 'large', length = 5.0)
    plt.show()
    plt.close()

# Plots voltage vs time during the input cycle
def plotVoltsVsTime(cycle):
    cycleVoltages = []
    cycleTimes = []
    
    # Stores the voltages and times during the input cycle in the cycleVoltages and
    # cycleTimes lists respectively
    if cycle*2 in cycles:
        for i in range(cycles.index(2*cycle-1), len(cycles) - cycles[::-1].index(2*cycle)):
            cycleVoltages.append(voltages[i])
            cycleTimes.append(times[i])
        plt.plot(cycleTimes, cycleVoltages, linewidth = 0.75)
        plt.title('Change in Voltage During Cycle %d' % cycle)
        plt.xlabel('Time (hr)')
        plt.ylabel('Voltage (V)')
    else:
        for i in range(cycles.index(2*cycle-1), len(cycles) - cycles[::-1].index(2*cycle-1)):
            cycleVoltages.append(voltages[i])
            cycleTimes.append(times[i])
        plt.plot(cycleTimes, cycleVoltages, linewidth = 0.75)
        plt.title('Change in Voltage During Cycle %d' % cycle)
        plt.xlabel('Time (hr)')
        plt.ylabel('Voltage (V)')
    plt.tick_params(direction='in', labelsize = 'large', length = 5.0)
    plt.show()
    plt.close()
    
# Plots voltage vs time for every cycle
def plotAllCyclesVoltsVsTime():
    plt.plot(times, voltages, linewidth = 0.75)
    plt.title('Change in Voltage During each Cycle')
    plt.xlabel('Time (hr)')
    plt.ylabel('Voltage (V)')
    plt.show()
    plt.close()
    
# Plots current vs time for the input channel and cycle
def plotCurrentVsTime(channel, cycle):
    cycleCurrents = []
    cycleTimes = []
    
    # Stores the currents and times for the input cycle in cycleCurrents and cycleTimes
    if cycle*2 in cycles:
        for i in range(cycles.index(2*cycle-1), len(cycles) - cycles[::-1].index(2*cycle)):
            cycleCurrents.append(channels[channel][i])
            cycleTimes.append(times[i])
        plt.plot(cycleTimes, cycleCurrents, linewidth = 0.75)
        plt.title('Change in Current on Channel %d During Cycle %d' % (channel, cycle))
        plt.xlabel('Time (hr)')
        plt.ylabel('Current (\u03BCA)')
    else:
        for i in range(cycles.index(2*cycle-1), len(cycles) - cycles[::-1].index(2*cycle-1)):
            cycleCurrents.append(channels[channel][i])
            cycleTimes.append(times[i])
        plt.plot(cycleTimes, cycleCurrents, linewidth = 0.75)
        plt.title('Change in Current on Channel %d During Cycle %d' % (channel, cycle))
        plt.xlabel('Time (hr)')
        plt.ylabel('Current (\u03BCA)')
    plt.show()
    plt.close()

# Plots current vs time for all the cycles on the input channel
def plotAllCyclesCurrentVsTime(channel):
    plt.plot(times, channels[channel], linewidth = 0.75)
    plt.title('Change in Current on Channel %d During each Cycle' % channel)
    plt.xlabel('Time (hr)')
    plt.ylabel('Current (\u03BCA)')
    plt.show()
    plt.close()

# Plots the voltage vs capacity for the input channel during the input cycle
def plotVoltsVsCapacity(channel, cycle):
    chargeCapacities = []
    chargeVoltages = []
    dischargeCapacities = []
    dischargeVoltages = []
    plt.figure(1)
    ax = plt.subplot(111)
    
    if(mass[channel] == 0.0):
        mass[channel] = float(input('Enter the sample mass for channel %d in grams: ' % channel))
    
    # Calculates specific capacity for the charge and discharge cycles and appends those
    # values to the cycleCapacities list and appends voltages during those cycles to the
    # cycleVoltages list
    if 2*cycle in cycles:
        for i in range(cycles.index(2*cycle-1), len(cycles) - cycles[::-1].index(2*cycle)):
            if(cycles[i] == 2*cycle - 1):
                if(i != len(cycles) - 1):
                    chargeCapacities.append(((0.5*(channels[channel][i] + channels[channel][i+1])*(times[i+1] - times[i]))/mass[channel])/1000.0)
                chargeVoltages.append(voltages[i])
            else:
                if(i != len(cycles) - 1):
                    dischargeCapacities.append(((0.5*(channels[channel][i] + channels[channel][i+1])*(times[i+1] - times[i]))/mass[channel])/1000.0)
                dischargeVoltages.append(voltages[i])
        chargePlot, = ax.plot(chargeCapacities, chargeVoltages, 'b', linewidth = 0.75, label = 'Charge')
        dischargePlot, = ax.plot(dischargeCapacities, dischargeVoltages, 'r', linewidth = 0.75, label = 'Discharge')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(handles = [chargePlot, dischargePlot], loc = 2, bbox_to_anchor = (1,1), fontsize = 'small')
        plt.title('Voltage as a Function of Specific Capacity on Channel %d' % channel)
        plt.xlabel('Specific Capacity (mAh/g)')
        plt.ylabel('Voltage (V)')
    else:
        for i in range(cycles.index(2*cycle-1), len(cycles) - cycles[::-1].index(2*cycle-1)):
            if(i != len(cycles) - 1):
                chargeCapacities.append(((0.5*(channels[channel][i] + channels[channel][i+1])*(times[i+1] - times[i]))/mass[channel])/1000.0)
            chargeVoltages.append(voltages[i])
        chargePlot, = ax.plot(chargeCapacities, chargeVoltages, 'b', linewidth = 0.75, label = 'Charge')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(handles = [chargePlot, dischargePlot], loc = 2, bbox_to_anchor = (1,1), fontsize = 'small')
        plt.title('Voltage as a Function of Specific Capacity on Channel %d' % channel)
        plt.xlabel('Specific Capacity (mAh/g)')
        plt.ylabel('Voltage (V)')
    plt.show()
    plt.close()

def plotCapacityVsCycle(channel):
    chgCapacity = []
    dchgCapacity = []
    cycle = []
    plt.figure(1)
    ax = plt.subplot(111)
    
    # Adds each cycle number to the cycle list
    for i in range(1, int(max(cycles)/2) + 1):
        cycle.append(i)
        
    # Adds the capacity for each cycle to the cycleCapacity list
    for j in cycle:
        charge = capacity(channel, j)['Charge capacity (mAh/g)']
        discharge = capacity(channel, j)['Discharge capacity (mAh/g)']
        chgCapacity.append(charge)
        dchgCapacity.append(-1.0*discharge)
    chgplot, = ax.plot(cycle, chgCapacity, 'bo', label = 'Charge')
    dchgplot, = ax.plot(cycle, dchgCapacity, 'rD', label = 'Discharge')
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
def plotChargeCapacityVsCycle(channel):
    cycleCapacity = []
    cycle = []
    
    # Adds each cycle number to the cycle list
    for i in range(1, int(max(cycles)/2) + 1):
        cycle.append(i)
        
    # Adds the capacity for each cycle to the cycleCapacity list
    for j in cycle:
        charge = capacity(channel, j)['Charge capacity (mAh/g)']
        cycleCapacity.append(charge)
    plt.plot(cycle, cycleCapacity, 'bo')
    for value in cycle:
        plt.annotate('%f' % cycleCapacity[value-1], xy = (value, cycleCapacity[value-1]), xytext = (value, cycleCapacity[value-1] + 0.007*cycleCapacity[value-1]))
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

def plotDischargeCapacityVsCycle(channel):
    cycleCapacity = []
    cycle = []
    
    # Adds each cycle number to the cycle list
    for i in range(1, int(max(cycles)/2) + 1):
        cycle.append(i)
        
    # Adds the capacity for each cycle to the cycleCapacity list
    for j in cycle:
        discharge = capacity(channel, j)['Discharge capacity (mAh/g)']
        cycleCapacity.append(-1.0*discharge)
    plt.plot(cycle, cycleCapacity, 'rD')
    for value in cycle:
        plt.annotate('%f' % cycleCapacity[value-1], xy = (value, cycleCapacity[value-1]), xytext = (value, cycleCapacity[value-1] + 0.007*cycleCapacity[value-1]))
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
print('capacity(channel, cycle)')
print('capacity64(cycle)')
print('capacityAllCycles(channel)')
print('allCapacities()')
print('averageVoltages(startTime, endTime)')
print('resistances(cycle)')
print('allResistances()')
print('plotCurrentVsVolts(channel, cycle)')
print('plotCurrentVsVolts64(cycle)')
print('plotAllCyclesCurrentVsVolts(channel)')
print('plotAllCyclesCurrentVsVolts64()')
print('plotVoltsVsCurrent(channel, cycle)')
print('plotVoltsVsCurrent64(cycle)')
print('plotAllCyclesVoltsVsCurrent(channel)')
print('plotAllCyclesVoltsVsCurrent64(channel)')
print('plotVoltsVsTime(cycle)')
print('plotAllCyclesVoltsVsTime()')
print('plotCurrentVsTime(channel, cycle)')
print('plotAllCyclesCurrentVsTime(channel)')
print('plotVoltsVsCapacity(channel, cycle)')
print('plotCapacityVsCycle(channel)')
print('plotChargeCapacityVsCycle(channel)')
print('plotDischargeCapacityVsCycle(channel)\n')

newFile()
# allCapacities()
# allResistances()
# export()
# capacity(1,1)
# capacityAllCycles(1)
# setMass(1)
# capacity64(1)
# averageVoltages(0, 13.01)
# resistances(1)
# plotVoltsVsCurrent(59, 1)
# plotVoltsVsCurrent64(1)
# plotCurrentVsVolts(59,1)
# plotCurrentVsVolts64(2)
# plotAllCyclesCurrentVsVolts(1)
# plotAllCyclesCurrentVsVolts64()
# plotAllCyclesVoltsVsCurrent64()
# plotVoltsVsTime(1)
# plotAllCyclesVoltsVsTime()
# plotChargeCapacityVsCycle(1)
# plotDischargeCapacityVsCycle(1)
# plotAllCyclesCurrentVsTime(1)
# plotCurrentVsTime(1, 1)
# plotVoltsVsCapacity(59, 1)
# plotCapacityVsCycle(1)