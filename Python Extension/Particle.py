class Particle:

    import numpy as np

    def __init__(self, energy, waveFunction):
        
        self.energy = energy

        self.waveFunction = waveFunction
        
    def setWaveFunction(self, waveFunction):
        self.waveFunction = waveFunction


    def getWaveFunction(self):
        return self.waveFunction

    def setEnergy(self, energy):
        self.energy = energy

    def getEnergy(self):
       return self.energy