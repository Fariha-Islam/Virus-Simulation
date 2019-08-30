from ps8_1 import *
from ps8_2 import virusCollection
from ps8_4 import histData

def simulationTwoDrugsVirusPopulations(numTrials = 1000):
    """
    Run simulations and plot graphs examining the relationship between
    administration of multiple drugs and patient outcome.
    Plots of total and drug-resistant viruses vs. time are made for asimulation with a 300 time step delay between administering the 2 drugs
    and
    a simulations for which drugs are administered simultaneously.
    """
    random.seed()

    # Virus Characteristics.
    maxPop = 1000
    numViruses = 100
    maxBirthProb = 0.1
    clearProb = 0.05

    numTimeStepsSimOne = 600
    firstSimMatrix = numpy.zeros(shape = (numTrials, numTimeStepsSimOne))
    gutResistVirusMatrixSimOne = numpy.zeros(shape = (numTrials, numTimeStepsSimOne))
    grimpResistVirusMatrixSimOne = numpy.zeros(shape = (numTrials, numTimeStepsSimOne))
    bothResistVirusMatrixSimOne = numpy.zeros(shape = (numTrials, numTimeStepsSimOne))

    numTimeStepsSimTwo = 300
    secondSimMatrix = numpy.zeros(shape = (numTrials, numTimeStepsSimTwo))
    gutResistVirusMatrixSimTwo = numpy.zeros(shape = (numTrials, numTimeStepsSimTwo))
    grimpResistVirusMatrixSimTwo = numpy.zeros(shape = (numTrials, numTimeStepsSimTwo))
    bothResistVirusMatrixSimTwo = numpy.zeros(shape = (numTrials, numTimeStepsSimTwo))

    for trial in range(numTrials):        

        # Model two random patients with the given virus charateristics.        
        virusesOne = virusCollection(numViruses, maxBirthProb, clearProb, ['guttagonol', 'grimpex'])
        virusesTwo = virusCollection(numViruses, maxBirthProb, clearProb, ['guttagonol', 'grimpex'])
        randPatientSimOne = Patient(virusesOne, maxPop)
        randPatientSimTwo = Patient(virusesTwo, maxPop)

        # Simulation One
        firstSimMatrix[trial][0] = numViruses
        for time in range(1, numTimeStepsSimOne):

            if time == 150:
                randPatientSimOne.addPrescription('guttagonol')
            if time == 450:
                randPatientSimOne.addPrescription('grimpex')

            firstSimMatrix[trial][time] = randPatientSimOne.update()
            gutResistVirusMatrixSimOne[trial][time] = randPatientSimOne.getResistPop(['guttagonol'])
            grimpResistVirusMatrixSimOne[trial][time] = randPatientSimOne.getResistPop(['grimpex'])
            bothResistVirusMatrixSimOne[trial][time] = randPatientSimOne.getResistPop(['guttagonol', 'grimpex'])
        
        # Simulation Two
        secondSimMatrix[trial][0] = numViruses
        for time in range(1, numTimeStepsSimTwo):

            if time == 150:
                randPatientSimTwo.addPrescription('guttagonol')
                randPatientSimTwo.addPrescription('grimpex')

            secondSimMatrix[trial][time] = randPatientSimTwo.update()
            gutResistVirusMatrixSimTwo[trial][time] = randPatientSimTwo.getResistPop(['guttagonol'])
            grimpResistVirusMatrixSimTwo[trial][time] = randPatientSimTwo.getResistPop(['grimpex'])
            bothResistVirusMatrixSimTwo[trial][time] = randPatientSimTwo.getResistPop(['guttagonol', 'grimpex'])
            
    # Statistical Analysis.

    # Simulation One
    totalMeanSimOne = firstSimMatrix.mean(0)
    gutResistMeanSimOne = gutResistVirusMatrixSimOne.mean(0)
    grimpResistMeanSimOne = grimpResistVirusMatrixSimOne.mean(0)
    bothResistMeanSimOne = bothResistVirusMatrixSimOne.mean(0)
    timeOne = numpy.arange(numTimeStepsSimOne) 

    # Simulation Two
    totalMeanSimTwo = secondSimMatrix.mean(0)
    gutResistMeanSimTwo = gutResistVirusMatrixSimTwo.mean(0)
    grimpResistMeanSimTwo = grimpResistVirusMatrixSimTwo.mean(0)
    bothResistMeanSimTwo = bothResistVirusMatrixSimTwo.mean(0)
    timeTwo = numpy.arange(numTimeStepsSimTwo) 

    # Plotting first Simulation
    pylab.figure(figsize=(14,12))
    pylab.subplot(211)
    pylab.plot(timeOne, totalMeanSimOne, label='Mean Virus Population', color='crimson', marker='o', linewidth=2.0, markevery=10)
    pylab.plot(timeOne, gutResistMeanSimOne, label='Guttagonol-only Resistant Virus Population', color='deepskyblue', marker='o', linewidth=2.0, markevery=10)
    pylab.plot(timeOne, grimpResistMeanSimOne, label='Grimpex-only Resistant Population', color='forestgreen', marker='o', linewidth=2.0, markevery=10)
    pylab.plot(timeOne, bothResistMeanSimOne, label='Both drugs Resistant Population', color='orchid', marker='o', linewidth=2.0, markevery=10)
    pylab.grid()    
    pylab.xlabel('Time Steps')
    pylab.ylabel('Mean Total Virus Population')
    pylab.title('Guttagonol administered after 150 timesteps. Grimpex administered after\n300 more timesteps. Total Simulation = 600 timesteps', fontsize='medium')
    pylab.legend(fontsize='x-small', loc='best')

    # Plotting second Simulation
    pylab.subplot(212)
    pylab.plot(timeTwo, totalMeanSimTwo, label='Mean Virus Population', color='crimson', marker='o', linewidth=2.0, markevery=10)
    pylab.plot(timeTwo, gutResistMeanSimTwo, label='Guttagonol-only Resistant Virus Population', color='deepskyblue', marker='o', linewidth=2.0, markevery=10)
    pylab.plot(timeTwo, grimpResistMeanSimTwo, label='Grimpex-only Resistant Population', color='forestgreen', marker='o', linewidth=2.0, markevery=10)
    pylab.plot(timeTwo, bothResistMeanSimTwo, label='Both drugs Resistant Population', color='orchid', marker='o', linewidth=2.0, markevery=10)
    pylab.grid()    
    pylab.xlabel('Time Steps')
    pylab.ylabel('Mean Total Virus Population')
    pylab.title('Guttagonol and Grimpex both simultaneously administered after 150 timesteps.\nTotal Simulation = 300 timesteps', fontsize='medium')
    pylab.legend(fontsize='x-small', loc='best')

    pylab.tight_layout()
    pylab.show()

if __name__ == '__main__':
    simulationTwoDrugsVirusPopulations()

