from ps8_1 import *
from ps8_2 import virusCollection

def simulationTwoDrugsDelayedTreatment(numTrials = 1000):
    """
    Runs simulations and make histograms for problem 6.
    Runs multiple simulations to show the relationship between administration
    of multiple drugs and patient outcome.
    Histograms of final total virus populations are displayed for lag times
    of
    150, 75, 0 timesteps between adding drugs (followed by an additional 150
    timesteps of simulation).
    """

    random.seed()

    # Virus Characteristics.
    maxPop = 1000
    numViruses = 100
    maxBirthProb = 0.1
    clearProb = 0.05

    viruses = virusCollection(numViruses, maxBirthProb, clearProb, ['guttagonol', 'grimpex'])
    bins = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550]

    finalPopSecondDrugAt_450, curedAfter_600 = histData(600, 150, 450, numTrials)
    finalPopSecondDrugAt_300, curedAfter_450 = histData(450, 150, 300, numTrials)
    finalPopSecondDrugAt_225, curedAfter_375 = histData(375, 150, 225, numTrials)
    finalPopSecondDrugAt_150, curedAfter_300 = histData(300, 150, 150, numTrials)

    pylab.figure(figsize=(15, 12))

    #Plotting 1st Histogram
    pylab.subplot(221)
    curedPercentage = curedAfter_600/float(numTrials)*100
    pylab.hist(finalPopSecondDrugAt_450, bins)
    pylab.title('Guttagonol Drug administered at 150 timesteps, Grimpex administered after\n300 more timesteps (Total = 600 timesteps)', fontsize='small')
    pylab.xlabel('Final Total Virus population ({}% cured)'.format(curedPercentage))
    pylab.ylabel('Number of Patients (Trials)')
    pylab.grid()

    #Plotting 2nd Histogram
    pylab.subplot(222)
    curedPercentage = curedAfter_450/float(numTrials)*100
    pylab.hist(finalPopSecondDrugAt_300, bins, color='red')
    pylab.title('Guttagonol Drug administered at 150 timesteps, Grimpex administered after\n150 more timesteps (Total = 450 timesteps)', fontsize='small')
    pylab.xlabel('Final Total Virus population ({}% cured)'.format(curedPercentage))
    pylab.ylabel('Number of Patients (Trials)')
    pylab.grid()

    #Plotting 3rd Histogram
    pylab.subplot(223)
    curedPercentage = curedAfter_375/float(numTrials)*100
    pylab.hist(finalPopSecondDrugAt_225, bins, color='magenta')
    pylab.title('Guttagonol Drug administered at 150 timesteps, Grimpex administered after\n75 more timesteps (Total = 375 timesteps)', fontsize='small')
    pylab.xlabel('Final Total Virus population ({}% cured)'.format(curedPercentage))
    pylab.ylabel('Number of Patients (Trials)')
    pylab.grid()

    #Plotting 4th Histogram
    pylab.subplot(224)
    curedPercentage = curedAfter_300/float(numTrials)*100
    pylab.hist(finalPopSecondDrugAt_150, bins, color='green')
    pylab.title('Guttagonol Drug administered at 150 timesteps, Grimpex administered after\n0 more timesteps (Total = 300 timesteps)', fontsize='small')
    pylab.xlabel('Final Total Virus population ({}% cured)'.format(curedPercentage))
    pylab.ylabel('Number of Patients (Trials)')
    pylab.grid()
    
    pylab.tight_layout()
    pylab.show()


def histData(totalSteps, firstDrugStep, secondDrugStep, numTrials):
    '''
    Helper method which returns a list of final virus population for each trial (patient)
    At drugStep, it administers the drug guttagonol.
    If final virus population at each trial is less than or equal to 50, increments cured
    counter by 1. Cured is set to 0 by default at the start of each function call.
    Returns a tuple: list of final virus populations, cured patients
    '''
    # Virus Characteristics.
    maxPop = 1000
    numViruses = 100
    maxBirthProb = 0.1
    clearProb = 0.05

    finalVirusList = []
    cured = 0

    for i in range(numTrials):
        viruses = virusCollection(numViruses, maxBirthProb, clearProb, ['guttagonol', 'grimpex'])
        patient = Patient(viruses, 1000)

        for step in range(totalSteps):
            if step == firstDrugStep:
                patient.addPrescription('guttagonol')
            if step == secondDrugStep:
                patient.addPrescription('grimpex')
            virusPop = patient.update()
        if virusPop <= 50:
            cured += 1       
        finalVirusList.append(virusPop)

    return finalVirusList, cured

if __name__ == '__main__':
    simulationTwoDrugsDelayedTreatment()

