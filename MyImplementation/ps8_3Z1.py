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

    ExptOneNumTimeSteps = 450
    ExptTwoNumTimeSteps = 300
    ExptThreeNumTimeSteps = 225
    ExptFourNumTimeSteps = 150
 
    ExptOneDataMatrix = numpy.zeros(shape = (numTrials, ExptOneNumTimeSteps))
    ExptTwoDataMatrix = numpy.zeros(shape = (numTrials, ExptTwoNumTimeSteps))
    ExptThreeDataMatrix = numpy.zeros(shape = (numTrials, ExptThreeNumTimeSteps)) 
    ExptFourDataMatrix = numpy.zeros(shape = (numTrials, ExptFourNumTimeSteps))
     
    for trial in range(numTrials): 
        # Model a random patient with the given virus charateristics.        
        virusesW = virusCollection(numViruses, maxBirthProb, clearProb)
        virusesX = virusCollection(numViruses, maxBirthProb, clearProb)
        virusesY = virusCollection(numViruses, maxBirthProb, clearProb)
        virusesZ = virusCollection(numViruses, maxBirthProb, clearProb)

        randPatientW = Patient(virusesW, maxPop)  
        randPatientX = Patient(virusesX, maxPop)  
        randPatientY = Patient(virusesY, maxPop)  
        randPatientZ = Patient(virusesZ, maxPop)  
        
        # The simulation runs for an additional 150 steps after administering guttagonol
        # Therefore, simulate the time-steps for:
        # (450,300), (300,150), (225,75), (150,0) 
        ExptOneDataMatrix[trial][0] = numViruses
        for time in range(1, ExptOneNumTimeSteps):
            if time == 300:
                randPatientW.addPrescription('guttagonol') # guttagonol added after 150 steps
            ExptOneDataMatrix[trial][time] = randPatientW.update()

        ExptTwoDataMatrix[trial][0] = numViruses
        for time in range(1, ExptTwoNumTimeSteps):
            if time == 150:
                randPatientX.addPrescription('guttagonol') # guttagonol added after 300 steps
            ExptTwoDataMatrix[trial][time] = randPatientX.update()
        
        ExptThreeDataMatrix[trial][0] = numViruses
        for time in range(1, ExptThreeNumTimeSteps):
            if time == 75:
                randPatientY.addPrescription('guttagonol') # guttagonol added after 75 steps
            ExptThreeDataMatrix[trial][time] = randPatientY.update()
        
        ExptFourDataMatrix[trial][0] = numViruses
        for time in range(0, ExptFourNumTimeSteps):
            if time == 0:
                randPatientZ.addPrescription('guttagonol') # guttagonol added after 0 steps (meaning at the beginning)
            ExptFourDataMatrix[trial][time] = randPatientZ.update()

    # Statistical Analysis.
    ExptOneMeanData = ExptOneDataMatrix.mean(0)
    ExptTwoMeanData = ExptTwoDataMatrix.mean(0)
    ExptThreeMeanData = ExptThreeDataMatrix.mean(0)
    ExptFourMeanData = ExptFourDataMatrix.mean(0)

    
    ExptOneTime = numpy.arange(ExptOneNumTimeSteps)
    ExptTwoTime = numpy.arange(ExptTwoNumTimeSteps)
    ExptThreeTime = numpy.arange(ExptThreeNumTimeSteps)
    ExptFourTime = numpy.arange(ExptFourNumTimeSteps)

    ExptOneStdData95_CI = ExptOneDataMatrix.std(0) * 2
    ExptTwoStdData95_CI = ExptTwoDataMatrix.std(0) * 2
    ExptThreeStdData95_CI = ExptThreeDataMatrix.std(0) * 2
    ExptFourStdData95_CI = ExptFourDataMatrix.std(0) * 2

    ExptOneSelectedTime = numpy.arange(0, ExptOneNumTimeSteps, 10)
    ExptTwoSelectedTime = numpy.arange(0, ExptTwoNumTimeSteps, 10)
    ExptThreeSelectedTime = numpy.arange(0, ExptThreeNumTimeSteps, 10)
    ExptFourSelectedTime = numpy.arange(0, ExptFourNumTimeSteps, 10)

    # Plotting 1st graph
    f = pylab.figure(figsize=(15, 15))

    pylab.subplot(221)
    pylab.plot(ExptOneTime, ExptOneMeanData)
    pylab.errorbar(ExptOneTime[ExptOneSelectedTime], ExptOneMeanData[ExptOneSelectedTime], ExptOneStdData95_CI[ExptOneSelectedTime], fmt = 'o')
    pylab.grid()    
    pylab.xlabel('Time Steps')
    pylab.ylabel('Total Virus Population')
    pylab.title('300 -> 150 = 450', fontsize='medium')

    # Plotting 2nd graph
    pylab.subplot(222)
    pylab.plot(ExptTwoTime, ExptTwoMeanData)
    pylab.errorbar(ExptTwoTime[ExptTwoSelectedTime], ExptTwoMeanData[ExptTwoSelectedTime], ExptTwoStdData95_CI[ExptTwoSelectedTime], fmt = 'o')
    pylab.grid()
    pylab.xlabel('Time Steps')
    pylab.ylabel('Total Guttagonol-Resistant Virus Population')
    pylab.title('150 -> 150 = 300', fontsize='medium')    

    # Plotting 3rd graph
    pylab.subplot(223)
    pylab.plot(ExptThreeTime, ExptThreeMeanData)
    pylab.errorbar(ExptThreeTime[ExptThreeSelectedTime], ExptThreeMeanData[ExptThreeSelectedTime], ExptThreeStdData95_CI[ExptThreeSelectedTime], fmt = 'o')
    pylab.grid()    
    pylab.xlabel('Time Steps')
    pylab.ylabel('Total Virus Population')
    pylab.title('75 -> 150 = 225', fontsize='medium')

    # Plotting 4th graph
    pylab.subplot(224)
    pylab.plot(ExptFourTime, ExptFourMeanData)
    pylab.errorbar(ExptFourTime[ExptFourSelectedTime], ExptFourMeanData[ExptFourSelectedTime], ExptFourStdData95_CI[ExptFourSelectedTime], fmt = 'o')
    pylab.grid()
    pylab.xlabel('Time Steps')
    pylab.ylabel('Total Guttagonol-Resistant Virus Population')
    pylab.title('0 -> 150 = 150', fontsize='medium')

    pylab.show()
    


if __name__ == '__main__':
    simulationDelayedTreatment()
