import numpy
import random
import pylab
from ps7_v1 import *

#Problem 1: Implementing a Simulation With Drugs

class ResistantVirus(SimpleVirus):


    """Representation of a virus which can have drug resistance.
    """

    def __init__(self, maxBirthProb, clearProb, resistances, mutProb):
        """
        Initialize a ResistantVirus instance, saves all parameters as
        attributes of the instance.
        maxBirthProb: Maximum reproduction probability (a float between 0-1)
        clearProb: Maximum clearance probability (a float between 0-1).
        resistances: A dictionary of drug names (strings) mapping to the
        state
        of this virus particle's resistance (either True or False) to each
        drug.
        e.g. {'guttagonol':False, 'grimpex',False}, means that this virus
        particle is resistant to neither guttagonol nor grimpex.
        mutProb: Mutation probability for this virus particle (a float). This
        is
        the probability of the offspring acquiring or losing resistance to a
        drug.
        """
        # Reference : https://stackoverflow.com/questions/2797139/python-using-methods-from-other-classes
        SimpleVirus.__init__(self, maxBirthProb, clearProb)
        #self.maxBirthProb = maxBirthProb
        #self.clearProb = clearProb
        self.resistances = resistances
        self.mutProb = mutProb

    def isResistantTo(self, drug):
        """
        Get the state of this virus particle's resistance to a drug. This
        method is called by getResistPop() in Patient to determine how many virus
        particles have resistance to a drug.
        drug: The drug (a string)
        returns: True if this virus instance is resistant to the drug, False
        otherwise.
        """
        return self.resistances.get(drug, False)  # searches for the drug in self.resistances, returns false if not found

    def reproduce(self, popDensity, activeDrugs):
        """
        Stochastically determines whether this virus particle reproduces at a
        time step. Called by the update() method in the Patient class.
        A virus particle will only reproduce if it is resistant to ALL the
        drugs
        in the activeDrugs list. For example, if there are 2 drugs in the
        activeDrugs list, and the virus particle is resistant to 1 or no
        drugs,then it will NOT reproduce.
        Hence, if the virus is resistant to all drugs
        in activeDrugs, then the virus reproduces with probability:
        self.maxBirthProb * (1 - popDensity).If this virus particle reproduces, then reproduce() creates and
        returns
        the instance of the offspring ResistantVirus (which has the same
        maxBirthProb and clearProb values as its parent).
        For each drug resistance trait of the virus (i.e. each key of
        self.resistances), the offspring has probability 1-mutProb of
        inheriting that resistance trait from the parent, and probability
        mutProb of switching that resistance trait in the offspring.
        For example, if a virus particle is resistant to guttagonol but not
        grimpex, and `self.mutProb` is 0.1, then there is a 10% chance that
        that the offspring will lose resistance to guttagonol and a 90%
        chance that the offspring will be resistant to guttagonol.
        There is also a 10% chance that the offspring will gain resistance to
        grimpex and a 90% chance that the offspring will not be resistant to
        grimpex.
        popDensity: the population density (a float), defined as the current
        virus population divided by the maximum population
        activeDrugs: a list of the drug names acting on this virus particle
        (a list of strings).
        returns: a new instance of the ResistantVirus class representing the
        offspring of this virus particle. The child should have the same
        maxBirthProb and clearProb values as this virus. Raises a
        NoChildException if this virus particle does not reproduce.
        """
        # loop through activeDrugs to see if the virus is not resistant to any
        # If it is not resistant to any, break out of loop and raise NoChildException
        for drug in activeDrugs:
            if not self.isResistantTo(drug):
                raise NoChildException('Not resistant to all drugs. Child not created!')
                break
        
        # If resistant to all, calculate proobability if it will reproduce
        maxReproduceProb = self.maxBirthProb * (1 - popDensity)
        
        # If it doesn't reproduce, raise NoChildException
        if random.random() < maxReproduceProb:
            childResistances = {}
            # Calculate probabilities to determine if child will be
            # resistant to the drug according to the given specification
            for drug, isResistant in self.resistances.items():
                if random.random() < (1-self.mutProb):
                    childResistances[drug] = isResistant
                else:
                    childResistances[drug] = not isResistant
            #Create new instance of child and return it                  
            childOfVirus = ResistantVirus(self.maxBirthProb, self.clearProb, childResistances, self.mutProb)
            return childOfVirus
        
        else: raise NoChildException('Child not created!')


# Patient Class

class Patient(SimplePatient):


    """
    Representation of a patient. The patient is able to take drugs and
    his/her
    virus population can acquire resistance to the drugs he/she takes.
    """

    def __init__(self, viruses, maxPop):
        """
        Initialization function, saves the viruses and maxPop parameters as
        attributes. Also initializes the list of drugs being administered(which should initially include no drugs).
        viruses: The list representing the virus population (a list of
        SimpleVirus instances)
        maxPop: The maximum virus population for this patient (an integer)
        """

        SimplePatient.__init__(self, viruses, maxPop)
        self.administeredDrugs = []

    def getVirusList(self):
        return self.viruses

    def addPrescription(self, newDrug):
        """
        Administer a drug to this patient. After a prescription is added, the
        drug acts on the virus population for all subsequent time steps. If
        the
        newDrug is already prescribed to this patient, the method has no
        effect.
        newDrug: The name of the drug to administer to the patient (a
        string).
        postcondition: The list of drugs being administered to a patient is
        updated
        """

        # newDrug to be administered to the patient is added to self.administeredDrugs
        self.administeredDrugs.append(newDrug)

    def getPrescriptions(self):
        """
        Returns the drugs that are being administered to this patient.
        returns: The list of drug names (strings) being administered to this
        patient.
        """

        return self.administeredDrugs

    def getResistPop(self, drugResist):
        """
        Get the population of virus particles resistant to the drugs listed
        in drugResist.
        drugResist: Which drug resistances to include in the population (a
        list of strings - e.g. ['guttagonol'] or ['guttagonol', 'grimpex'])
        returns: The population of viruses (an integer) with resistances to
        all drugs in the drugResist list.
        """
        resistantVirusPop = 0  # initializing the population of virus particles resistant to the drugs
        for virus in self.viruses:
        # contains list of resistant drugs (for the viruses)
            noOfDrugResist = 0
            # for each resistive drug in the drugResist list
            for drug in drugResist:
                if virus.isResistantTo(drug):
                    # if the virus is resistant to the drug, we increment the resistant drug list
                    noOfDrugResist += 1
            if noOfDrugResist == len(drugResist):
                resistantVirusPop += 1

        return resistantVirusPop
    
    
    def update(self):
        """
        Update the state of the virus population in this patient for a single
        time step. update() should execute these actions in order:
        - Determine whether each virus particle survives and update the list
        of virus particles accordingly
        - The current population density is calculated. This population
        density value is used until the next call to update().
        - Determine whether each virus particle should reproduce and add
        offspring virus particles to the list of viruses in this patient.
        The list of drugs being administered should be accounted for in the
        determination of whether each virus particle reproduces.
        returns: The total virus population at the end of the update (an
        integer)
        """
        # Determine number of viruses to be cleaned, "stochastically".
        numRemoveVirus = 0
        for virus in self.viruses:
            if virus.doesClear():
                numRemoveVirus += 1

        # Remove numRemoveVirus from the patient's body.
        for virusNum in range(numRemoveVirus):
            self.viruses.pop()

        # Calculate population density. TO DO check (keep self!)
        popDensity = self.getTotalPop() / float(self.maxPop)

        if popDensity >= 1:
            print 'virus population reached maximum!'
            popDensity = 1

        # Reproduce at a single time step.
        # The list of drugs being administered is being accounted for in the determination of whether each virus particle reproduces.
        offspringViruses = []
        for virus in self.viruses:
            try:
                offspringViruses.append(virus.reproduce(
                    popDensity, self.administeredDrugs))
            except NoChildException:
                pass

        self.viruses = self.viruses + offspringViruses

        return self.getTotalPop()
    
# Problem 2
def virusCollection(numViruses, maxBirthProb, clearProb):
    viruses = []
    resistances = {}
    resistances['guttagonol'] = False
    for virusNum in range(numViruses):
        viruses.append(ResistantVirus(maxBirthProb, clearProb, resistances, 0.005))
    return viruses

def simulationWithDrug(numTrials = 30, numTimeSteps = 300):

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
        viruses = virusCollection(numViruses, maxBirthProb, clearProb)
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

    f = pylab.figure(figsize=(15, 7))

    # Ploting.
    pylab.subplot(121)
    pylab.plot(time, meanData)
    pylab.errorbar(time[selectedTime], meanData[selectedTime], stdData95_CI[selectedTime], fmt = 'o')
    pylab.grid()    
    pylab.xlabel('Time Steps')
    pylab.ylabel('Total Virus Population')
    pylab.title('Effect of Guttagonol on Virus Population being administered\nafter {} Timesteps over a total period of {} Timesteps'.format('150', '300'), fontsize='medium')

    stdDevGutVirusPop = gutResistVirusMatrix.std(0) * 2

    # Plotting 2nd graph
    pylab.subplot(122)
    pylab.plot(time, meanResistVirus)
    pylab.errorbar(time[selectedTime], meanResistVirus[selectedTime], stdDevGutVirusPop[selectedTime], fmt = 'x')
    pylab.grid()
    pylab.xlabel('Time Steps')
    pylab.ylabel('Total Guttagonol-Resistant Virus Population')
    pylab.title('Total Number of Gluttagonol-Resistant Virus Population after {} Timesteps\nDrug administered after {} Timesteps'.format('300', '150'), fontsize='medium')
    pylab.show()
    

simulationWithDrug()


    