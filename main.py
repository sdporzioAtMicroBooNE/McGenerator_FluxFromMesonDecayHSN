# Import awesome stuff
from classes.detector import Detector
from classes.particle import Particle
import json, datetime
import pandas as pd
pdf = pd.read_csv('meson_data/mesonSpectrum.csv')
geom = json.loads(open('geometry/detectors.json').read())

# Create detector
detName = 'uboone'
detector = Detector(geom[detName]['length'],geom[detName]['width'],geom[detName]['height'],geom[detName]['xd'],geom[detName]['yd'],geom[detName]['zd'])

# Variables for verbose
ldf = [] # List dataframe for neutrinos
divisions = 100. # Number of checkpoints for verbose
nEvents = len(pdf) # Number of events to process
nEventsCheck = int(nEvents/divisions) # Number of events per checkpoint
nNuK = 0 # Number of neutrinos from kaons
nNuKCheck = 0 # Number of neutrinos from kaons per checkpoint
nNuPi = 0 # Number of neutrinos from pions
nNuPiCheck = 0 # Number of neutrinos from pions per checkpoint
iCheck = 0 # Checkpoint index
timeTotal = 0 # Total elapsed time
startCheck = datetime.datetime.now() # Time at starts of loop
prevCheck = startCheck # Time at previous loop instance
currCheck = startCheck # Time at current loop instance

# Loop through each event
for p in range(nEvents):
    meson = Particle(pdf['Pdg'][p],pdf['Vx'][p],pdf['Vy'][p],pdf['Vz'][p],pdf['Px'][p],pdf['Py'][p],pdf['Pz'][p])
    mu, nu = meson.DecayParticle(13,14)
    if detector.IntersectTrajectory(nu):
        ldf.append([meson.GetPDG(),nu.GetPDG(),
        nu.GetStartPoint().GetX(),nu.GetStartPoint().GetY(),nu.GetStartPoint().GetZ(),
        nu.GetStartPoint().GetFourM().GetPx(),nu.GetStartPoint().GetFourM().GetPy(),nu.GetStartPoint().GetFourM().GetPz()
        ])
        if meson.GetPDG()==211:
            nNuPiCheck += 1
        if meson.GetPDG()==321:
            nNuKCheck += 1

    # Periodic check
    if p%nEventsCheck == 0 and p!=0:
        iCheck += 1
        remainCheck = divisions-iCheck
        currCheck = datetime.datetime.now()
        timeCheck = int((currCheck-prevCheck).total_seconds())
        timeTotal += timeCheck
        timeAverage = timeTotal/float(iCheck)
        timeEta = timeAverage*remainCheck
        nNuK += nNuKCheck
        nNuPi += nNuPiCheck
        print currCheck
        print '%i%% completed (%i events of %i)' %(iCheck*int(100/divisions),p,nEvents)
        print '%i new neutrinos reached the detector: %i Pi+ | %i K+' %(nNuKCheck+nNuPiCheck,nNuPiCheck,nNuKCheck)
        print 'Current total neutrino count: %i Pi+ | %i K+' %(nNuPi,nNuK)
        print 'Last %i events were processed in %s s' %(nEventsCheck,timeCheck)
        print 'ETA: %s' %(str(datetime.timedelta(seconds=timeEta)))
        print
        prevCheck = currCheck
        nNuKCheck = 0
        nNuPiCheck = 0


    # Clean up rubbish
    del meson
    del mu
    del nu

df = pd.DataFrame(ldf,columns=['ParentPDG','PDG','Vx','Vy','Vz','Px','Py','Pz'])
df.to_csv('data.csv',index=False)
