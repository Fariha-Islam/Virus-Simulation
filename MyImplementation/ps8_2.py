from ps8_1 import *

'''
Problem 2 Implementation
Simulation with administration of 1 drug: Gutaggonol after 150 Timesteps
'''

def virusCollection(numViruses, maxBirthProb, clearProb, drugList):
    viruses = []
    resistances = {}
    for drug in drugList:
        resistances[drug] = False
    for virusNum in range(numViruses):
        viruses.append(ResistantVirus(maxBirthProb, clearProb, resistances, 0.005))
    return viruses

def simulationWithDrug(numTrials = 100, numTimeSteps = 300):

    """
    Run the simulation and plot the graph for problem 2 (Guttagonol drug is used,
    viruses do not have any drug resistance initially). Mutation probability = 0.005
    Instantiates a Patient, runs a simulation for 300 timesteps, Guttagonol is 
    administered at the 150th timestep. Graph is plotted with the total virus population
    as a function of time.
    """
    random.seed()

    # Virus Characteristics.
    maxPop = 1000
    numViruses = 100
    maxBirthProb = 0.1
    clearProb = 0.05
    
    gutResistVirusMatrix = numpy.zeros(shape = (numTrials, numTimeSteps))
    dataMatrix = numpy.zeros(shape = (numTrials, numTimeSteps))    
    for trial in range(numTrials):        

        # Model a random patient with the given virus charateristics.        
        viruses = virusCollection(numViruses, maxBirthProb, clearProb, ['guttagonol'])
        randPatientX = Patient(viruses, maxPop)

        # Simulate the time-steps.
        dataMatrix[trial][0] = numViruses
        for time in range(1, numTimeSteps):
            if time == 150:
                randPatientX.addPrescription('guttagonol')
            dataMatrix[trial][time] = randPatientX.update()
            gutResistVirusMatrix[trial][time] = randPatientX.getResistPop(['guttagonol'])            
            
    # Statistical Analysis.
    meanData = dataMatrix.mean(0)
    time = numpy.arange(numTimeSteps) 
    stdData95_CI = dataMatrix.std(0) * 2
    selectedTime = numpy.arange(0, numTimeSteps, 10)

    meanResistVirus = gutResistVirusMatrix.mean(0)

    #f = pylab.figure(figsize=(15, 7))

    # Plotting.
    #pylab.subplot(121)
    pylab.plot(time, meanData, label='Mean Virus Population')
    pylab.errorbar(time[selectedTime], meanData[selectedTime], stdData95_CI[selectedTime], fmt = 'o', color = 'blue')
    pylab.grid()    
    pylab.xlabel('Time Steps')
    pylab.ylabel('Total Virus Population')
    pylab.title('Effect of Guttagonol on Virus Population being administered\nafter {} Timesteps over a total period of {} Timesteps'.format('150', '300'), fontsize='medium')

    stdDevGutVirusPop = gutResistVirusMatrix.std(0) * 2

    # Plotting 2nd graph
    #pylab.subplot(122)
    pylab.plot(time, meanResistVirus, label='Mean Guttagonol-resistant Virus Population', color = 'red')
    pylab.errorbar(time[selectedTime], meanResistVirus[selectedTime], stdDevGutVirusPop[selectedTime], fmt = 'o', color = 'red')
    pylab.legend(fontsize='x-small', loc='best')
    #pylab.grid()
    #pylab.xlabel('Time Steps')
    #pylab.ylabel('Total Guttagonol-Resistant Virus Population')
    #pylab.title('Total Number of Guttagonol-Resistant Virus Population after {} Timesteps\nDrug administered after {} Timesteps'.format('300', '150'), fontsize='medium')
    pylab.show()
    

if __name__ == '__main__':
    simulationWithDrug()