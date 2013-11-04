from CMGTools.RootTools.physicsobjects.PhysicsObject import *

class Candidate( PhysicsObject ):
    
    def fromPV ( self ):
        '''bool to check if the candidate comes from PV'''
        return self.fromPV(),
