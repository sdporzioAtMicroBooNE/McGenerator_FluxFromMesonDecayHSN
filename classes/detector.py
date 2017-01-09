import json, math
import numpy as np
from collections import OrderedDict as OD
const = json.loads(open('constants/constants.json').read())
pdgKey = OD([(y,x) for (x,y) in const['pdg'].items()])
c = const['constants']['c']

class Detector:
    '''Builds 6 faces for the detector that will be checked for intersection with trajectory.
    Requires length, width, and height of the detector and the coordinates in the world frame of the center
    of the detector.'''
    def __init__(self,length,width,height,xc,yc,zc):
        self.length = length
        self.width = width
        self.height = height
        self.xc = xc
        self.yc = yc
        self.zc = zc
        self.bound = [[xc-width/2.,yc-height/2.,zc-length/2.],
                 [xc+width/2.,yc+height/2.,zc+length/2.]]
    def IntersectTrajectory(self,particle):
        'Determine if a particle will intersect with the detector'
        coord = [particle.startPoint.x,particle.startPoint.y,particle.startPoint.z]
        speed = [particle.startPoint.fourm.betax*c,particle.startPoint.fourm.betay*c,particle.startPoint.fourm.betaz*c]

        otherInd = [[1,2],[0,2],[0,1]]
        intersBool = []
        for i in range(2): # Looping over low and high boundaries
            for j in range(3): # Looping over each dimension
                if speed[j] == 0:
                    intersBool.append(False)
                else:
                    intersTime = (self.bound[i][j] - coord[j])/float(speed[j])
                    coord1 = coord[otherInd[j][0]] + speed[otherInd[j][0]]*intersTime
                    coord2 = coord[otherInd[j][1]] + speed[otherInd[j][1]]*intersTime
                    bool1_low = coord1 > self.bound[0][otherInd[j][0]]
                    bool1_high = coord1 < self.bound[1][otherInd[j][0]]
                    bool2_low = coord2 > self.bound[0][otherInd[j][1]]
                    bool2_high = coord2 < self.bound[1][otherInd[j][1]]
                    result = bool1_low and bool1_high and bool2_low and bool2_high
                    intersBool.append(result)
        if any(face == True for face in intersBool):
            return True
        else:
            return False
