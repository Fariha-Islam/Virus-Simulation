from ps8_1 import *

def virusCollection(numViruses, maxBirthProb, clearProb):
    viruses = []
    resistances = {}
    resistances['guttagonol'] = False
    for virusNum in range(numViruses):
        viruses.append(ResistantVirus(maxBirthProb, clearProb, resistances, 0.005))
    return viruses 

def simulationDelayedTreatment(numTrials = 30):
"""Runs simulations and make histograms for problem 5.
Runs multiple simulations to show the relationship between delayed
treatment
and patient outcome.
Histograms of final total virus populations are displayed for delays of
300, 150, 75, 0 timesteps (followed by an additional 150 timesteps of
simulation).
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
        viruses = virusCollection(numViruses, maxBirthProb, clearProb)
        randPatientX = Patient(viruses, maxPop)  
        
        # The simulation runs for an additional 150 steps after administering guttagonol
        # Therefore, simulate the time-steps for:
        # (450,300), (300,150), (225,75), (150,0) 
        dataMatrix[trial][0] = numViruses
        for time in range(1, 300):
            if time == 150:
                randPatientX.addPrescription('guttagonol') # guttagonol added after 150 steps
            dataMatrix[trial][time] = randPatientX.update()
            gutResistVirusMatrix[trial][time] = randPatientX.getResistPop(['guttagonol'])

        for time in range(1, 450):
            if time == 300:
                randPatientX.addPrescription('guttagonol') # guttagonol added after 300 steps
            dataMatrix[trial][time] = randPatientX.update()
            gutResistVirusMatrix[trial][time] = randPatientX.getResistPop(['guttagonol'])

        for time in range(1, 225):
            if time == 75:
                randPatientX.addPrescription('guttagonol') # guttagonol added after 75 steps
            dataMatrix[trial][time] = randPatientX.update()
            gutResistVirusMatrix[trial][time] = randPatientX.getResistPop(['guttagonol'])
        
        for time in range(1, 150):
            if time == 0:
                randPatientX.addPrescription('guttagonol') # guttagonol added after 0 steps (meaning at the beginning)
            dataMatrix[trial][time] = randPatientX.update()
            gutResistVirusMatrix[trial][time] = randPatientX.getResistPop(['guttagonol'])


        

#if __name__ == '__main__':
    #simulationDelayedTreatment()
