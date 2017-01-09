import json, math
import numpy as np
from collections import OrderedDict as OD
const = json.loads(open('constants/constants.json').read())
pdgKey = OD([(y,x) for (x,y) in const['pdg'].items()])
c = const['constants']['c']

def Energy(mass,momentum):
    return math.sqrt(pow(mass,2)+pow(momentum,2))
def Mass(energy,momentum):
    return math.sqrt(pow(energy,2)-pow(momentum,2))
def Momentum(energy,mass):
    return math.sqrt(pow(energy,2)-pow(mass,2))

class FourVector:
    'Four-Momentum vector'
    def __init__(self,e,px,py,pz):
        self.e = e
        self.px = px
        self.py = py
        self.pz = pz
        self.p = math.sqrt(pow(px,2) + pow(py,2) + pow(pz,2))
        self.beta = self.p/float(self.e)
        self.betax = self.px/float(self.e)
        self.betay = self.py/float(self.e)
        self.betaz = self.pz/float(self.e)
        m = Mass(self.e,self.p)
        if m!=0:
            self.gamma = self.e/float(m)
        else:
            self.gamma = float('nan')
    def Print(self):
        print self.e, self.px, self.py, self.pz
    def GetE(self):
        return self.e
    def GetPx(self):
        return self.px
    def GetPy(self):
        return self.py
    def GetPz(self):
        return self.pz
    def GetP(self):
        return self.p
    def GetGamma(self):
        return self.gamma
    def GetBeta(self):
        return self.beta
    def GetBetaX(self):
        return self.betax
    def GetBetaY(self):
        return self.betay
    def GetBetaZ(self):
        return self.betaz
    def LorentzBoost(self,fourVec):
        gamma = fourVec.GetGamma()
        beta = fourVec.GetBeta()
        betax = fourVec.GetBetaX()
        betay = fourVec.GetBetaY()
        betaz = fourVec.GetBetaZ()

        colVec = np.matrix([
                [self.e],
                [self.px],
                [self.py],
                [self.pz]
            ])
        lorentzMat = np.matrix([
                [gamma, -betax*gamma, -betay*gamma, -betaz*gamma],
                [-betax*gamma,1+(gamma-1)*(betax*betax)/(beta*beta),(gamma-1)*(betax*betay)/(beta*beta),(gamma-1)*(betax*betaz)/(beta*beta)],
                [-betay*gamma,(gamma-1)*(betay*betax)/(beta*beta),1+(gamma-1)*(betay*betay)/(beta*beta),(gamma-1)*(betay*betaz)/(beta*beta)],
                [-betaz*gamma,(gamma-1)*(betaz*betax)/(beta*beta),(gamma-1)*(betaz*betay)/(beta*beta),1+(gamma-1)*(betaz*betaz)/(beta*beta)]
            ])

        resultVec = lorentzMat*colVec
        self.e = resultVec[0,0]
        self.px = resultVec[1,0]
        self.py = resultVec[2,0]
        self.pz = resultVec[3,0]
