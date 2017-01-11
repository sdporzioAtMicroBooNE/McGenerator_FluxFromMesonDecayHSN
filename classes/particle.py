from fourmomentum import *

import json, math
import numpy as np
from collections import OrderedDict as OD
const = json.loads(open('constants/constants.json').read())
pdgKey = OD([(y,x) for (x,y) in const['pdg'].items()])
c = const['constants']['c']

def DecayFunction(decayConstant,t):
    'Exponential decay function'
    return np.exp(-decayConstant*t)

class Point:
    'Particle-point corresponding to a specific moment of the particle lifetime'
    def __init__(self,m,x,y,z,px,py,pz,t):
        self.x = x
        self.y = y
        self.z = z
        self.t = t
        p = math.sqrt(pow(px,2) + pow(py,2) + pow(pz,2))
        e = Energy(m,p)
        self.fourm = FourVector(e,px,py,pz)
    def Print(self):
        print self.t, self.x, self.y, self.z
    def GetX(self):
        return self.x
    def GetY(self):
        return self.y
    def GetZ(self):
        return self.z
    def GetFourM(self):
        return self.fourm
    def GetT(self):
        return self.t

class Particle:
    'Common base class for all particles'
    def __init__(self,pdg,x,y,z,px,py,pz):
        self.pdg = pdg
        self.m = const['mass'][pdgKey[self.pdg]]
        self.startPoint = Point(self.m,x,y,z,px,py,pz,0)
        self.decayStatus = False
    def GetPDG(self):
        return self.pdg
    def GetMass(self):
        return self.m
    def GetStartPoint(self):
        return self.startPoint
    def GetEndPoint(self):
        if self.decayStatus == True:
            return self.endPoint
        else:
            raise Exception('Particle has not been decayed yet!')
    def GetDecayStatus(self):
        return self.decayStatus
    def GenerateEndPoint(self):
        'Simulate a random decay and generates end-point for particle'
        meanLife = const['mean_life'][pdgKey[self.pdg]]
        decayConst = 1/float(meanLife)
        while not self.decayStatus:
            t = np.random.rand()*meanLife*10.
            p = np.random.rand()
            if p < DecayFunction(decayConst,t):
                self.decayStatus = True
        # Determine end coordinates for decayed particle and create end-point
        endX = self.startPoint.x + self.startPoint.fourm.px*c*t/float(self.m)
        endY = self.startPoint.y + self.startPoint.fourm.py*c*t/float(self.m)
        endZ = self.startPoint.z + self.startPoint.fourm.pz*c*t/float(self.m)
        self.endPoint = Point(self.m,endX,endY,endZ,self.startPoint.fourm.px,self.startPoint.fourm.py,self.startPoint.fourm.pz,t)
    def DecayParticle(self,pdg1,pdg2):
        'Generate particles from parent meson decay'
        # Check if end-point has been generated already and generate if not
        if self.decayStatus == False:
            self.GenerateEndPoint()
        # Generate random direction in space for decay particle trajectories (isotropic back-to-back)'
        theta1 = np.random.rand()*np.pi
        phi1 = np.random.rand()*np.pi*2.
        theta2 = theta1 + np.pi
        phi2 = phi1
        # Determine masses from pdg code
        mass1 = const['mass'][pdgKey[pdg1]]
        mass2 = const['mass'][pdgKey[pdg2]]
        # Determine kinematically energies and momenta for the decay particles
        e1 = (pow(self.m,2) + pow(mass1,2) - pow(mass2,2))/float(2.*self.m)
        e2 = (pow(self.m,2) + pow(mass2,2) - pow(mass1,2))/float(2.*self.m)
        p1 = Momentum(e1,mass1)
        p2 = Momentum(e2,mass2)
        # Assign cartesian components of momentum
        px1 = p1*np.sin(theta1)*np.cos(phi1)
        py1 = p1*np.sin(theta1)*np.sin(phi1)
        pz1 = p1*np.cos(theta1)
        px2 = p2*np.sin(theta2)*np.cos(phi2)
        py2 = p2*np.sin(theta2)*np.sin(phi2)
        pz2 = p2*np.cos(theta2)
        # Register decay particles (in rest frame of parent)
        part1 = Particle(pdg1,self.endPoint.x,self.endPoint.y,self.endPoint.z,px1,py1,pz1)
        part2 = Particle(pdg2,self.endPoint.x,self.endPoint.y,self.endPoint.z,px2,py2,pz2)
        # Boost decay particles to lab frame
        part1.GetStartPoint().GetFourM().LorentzBoost(self.startPoint.fourm)
        part2.GetStartPoint().GetFourM().LorentzBoost(self.startPoint.fourm)
        return part1, part2
