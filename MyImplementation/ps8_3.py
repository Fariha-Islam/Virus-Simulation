from ps8_1 import *
from ps8_2 import virusCollection

def simulationDelayedTreatment(numTrials = 10):
    """Runs simulations and make histograms for problem 5.
    Runs multiple simulations to show the relationship between delayed
    treatment
    and patient outcome.
    Histograms of final total virus populations are displayed for delays of
    300, 150, 75, 0 timesteps (followed by an additional 150 timesteps of
    simulation).
    """
    random.seed()

    bins = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550]

    finalPopDrugAt_300, curedAfter_450 = histData(450, 300, numTrials)
    finalPopDrugAt_150, curedAfter_300 = histData(300, 150, numTrials)
    finalPopDrugAt_75, curedAfter_225 = histData(225, 75, numTrials)
    finalPopDrugAt_0, curedAfter_150 = histData(150, 0, numTrials)

    pylab.figure(figsize=(10, 10))

    #Plotting 1st Histogram
    pylab.subplot(221)
    curedPercentage = curedAfter_450/float(numTrials)*100
    pylab.hist(finalPopDrugAt_300, bins)
    pylab.title('Drug administered at\n300 timesteps, followed by additional 150 timesteps (Total = 450 timesteps)', fontsize='small')
    pylab.xlabel('Final Total Virus population ({}% cured)'.format(curedPercentage))
    pylab.ylabel('Number of Patients (Trials)')
    pylab.grid()

    #Plotting 2nd Histogram
    pylab.subplot(222)
    curedPercentage = curedAfter_300/float(numTrials)*100
    pylab.hist(finalPopDrugAt_150, bins, color='red')
    pylab.title('Drug administered at\n150 timesteps, followed by additional 150 timesteps (Total = 300 timesteps)', fontsize='small')
    pylab.xlabel('Final Total Virus population ({}% cured)'.format(curedPercentage))
    pylab.ylabel('Number of Patients (Trials)')
    pylab.grid()

    #Plotting 3rd Histogram
    pylab.subplot(223)
    curedPercentage = curedAfter_225/float(numTrials)*100
    pylab.hist(finalPopDrugAt_75, bins, color='magenta')
    pylab.title('Drug administered at\n75 timesteps, followed by additional 150 timesteps (Total = 225 timesteps)', fontsize='small')
    pylab.xlabel('Final Total Virus population ({}% cured)'.format(curedPercentage))
    pylab.ylabel('Number of Patients (Trials)')
    pylab.grid()

    #Plotting 4th Histogram
    pylab.subplot(224)
    curedPercentage = curedAfter_150/float(numTrials)*100
    pylab.hist(finalPopDrugAt_0, bins, color='green')
    pylab.title('Drug administered at\n0 timesteps, followed by additional 150 timesteps (Total = 150 timesteps)', fontsize='small')
    pylab.xlabel('Final Total Virus population ({}% cured)'.format(curedPercentage))
    pylab.ylabel('Number of Patients (Trials)')
    pylab.grid()
    
    pylab.tight_layout()
    pylab.show()


def histData(totalSteps, drugStep, numTrials):
    '''
    Helper method which returns a list of final virus population for each trial (patient)
    At drugStep, it administers the drug guttagonol
    If final virus population at each trial is less than or equal to 50, increments cured
    counter by 1. Cured is set to 0 by default at the start of each function call.
    '''
    # Virus Characteristics.
    maxPop = 1000
    maxBirthProb = 0.1
    clearProb = 0.05
    numViruses = 100

    finalVirusList = []
    cured = 0

    for i in range(numTrials):
        viruses = virusCollection(numViruses, maxBirthProb, clearProb, ['guttagonol'])
        patient = Patient(viruses, 1000)

        for step in range(totalSteps):
            if step == drugStep:
                patient.addPrescription('guttagonol')
            virusPop = patient.update()
        if virusPop <= 50:
            cured += 1       
        finalVirusList.append(virusPop)

    return finalVirusList, cured

if __name__ == '__main__':
    simulationDelayedTreatment()
