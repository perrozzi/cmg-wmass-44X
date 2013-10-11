import numpy as my_n
import math
from CMGTools.RootTools.fwlite.AutoHandle import AutoHandle
from CMGTools.H2TauTau.proto.analyzers.ntuple import *
from CMGTools.RootTools.physicsobjects.PhysicsObjects import Lepton, Muon, Electron
from CMGTools.RootTools.analyzers.TreeAnalyzerNumpy import TreeAnalyzerNumpy
from CMGTools.RootTools.physicsobjects.PileUpSummaryInfo import PileUpSummaryInfo

def var( tree, varName, type=float ):
    tree.var(varName, type)

def fill( tree, varName, value ):
    tree.fill( varName, value )

def bookParticle( tree, pName ):
    var(tree, '{pName}_pt'.format(pName=pName))
    var(tree, '{pName}_eta'.format(pName=pName))
    var(tree, '{pName}_phi'.format(pName=pName))
    var(tree, '{pName}_mass'.format(pName=pName))
    var(tree, '{pName}_charge'.format(pName=pName))

def bookW( tree, pName ):
    var(tree, '{pName}_pt'.format(pName=pName))
    var(tree, '{pName}_phi'.format(pName=pName))

def bookMET( tree, pName ):
    var(tree, '{pName}'.format(pName=pName))
    var(tree, '{pName}_phi'.format(pName=pName))
    var(tree, '{pName}_sumEt'.format(pName=pName))

def bookMET2( tree, pName ):
    var(tree, '{pName}'.format(pName=pName))
    var(tree, '{pName}_phi'.format(pName=pName))

def bookJet( tree, pName ):
    var(tree, '{pName}_pt'.format(pName=pName))
    var(tree, '{pName}_eta'.format(pName=pName))
    var(tree, '{pName}_phi'.format(pName=pName))

def bookMuonCovMatrix( tree, pName ):
    tree.vars['{pName}CovMatrix'.format(pName=pName)]= my_n.zeros(9, dtype=float)
    tree.tree.Branch('{pName}CovMatrix'.format(pName=pName),tree.vars['{pName}CovMatrix'.format(pName=pName)] ,'{pName}CovMatrix'.format(pName=pName)+'[9]/D' )

def fillParticle( tree, pName, particle ):
    fill(tree, '{pName}_pt'.format(pName=pName), particle.pt() )
    fill(tree, '{pName}_eta'.format(pName=pName), particle.eta() )
    fill(tree, '{pName}_phi'.format(pName=pName), particle.phi() )
    fill(tree, '{pName}_mass'.format(pName=pName), particle.mass() )
    fill(tree, '{pName}_charge'.format(pName=pName), particle.charge() )

def fillW( tree, pName, particle ):
    fill(tree, '{pName}_pt'.format(pName=pName), particle.Pt() )
    fill(tree, '{pName}_phi'.format(pName=pName), particle.Phi() )

def fillMET( tree, pName, particle ):
    fill(tree, '{pName}'.format(pName=pName), particle.pt() )
    fill(tree, '{pName}_phi'.format(pName=pName), particle.phi() )
    fill(tree, '{pName}_sumEt'.format(pName=pName), particle.sumEt() )

def fillMET2( tree, pName, particle ):
    fill(tree, '{pName}'.format(pName=pName), particle.pt() )
    fill(tree, '{pName}_phi'.format(pName=pName), particle.phi() )
    
def fillJet( tree, pName, particle ):
    fill(tree, '{pName}_pt'.format(pName=pName), particle.pt() )
    fill(tree, '{pName}_eta'.format(pName=pName), particle.eta() )
    fill(tree, '{pName}_phi'.format(pName=pName), particle.phi() )

def fillMuonCovMatrix( tree, pName, covMatrix,event ):
    for i in range(0,9):
        tree.vars['{pName}CovMatrix'.format(pName=pName)][i] = covMatrix[i]



class WTreeProducer( TreeAnalyzerNumpy ):
    
    MuonClass = Muon 

    
    def declareHandles(self):
      super(WTreeProducer, self).declareHandles()
    
      self.handles['pfMet'] = AutoHandle('cmgPFMET','std::vector<cmg::BaseMET>')
      self.handles['pfMetraw'] = AutoHandle('cmgPFMETRaw','std::vector<cmg::BaseMET>')
      self.handles['pfMetSignificance'] = AutoHandle('pfMetSignificance','cmg::METSignificance')
      self.handles['nopuMet'] = AutoHandle('nopuMet','std::vector<reco::PFMET>')
      self.handles['pucMet'] = AutoHandle('pcMet','std::vector<reco::PFMET>')
      self.handles['pfMetForRegression'] = AutoHandle('pfMetForRegression','std::vector<reco::PFMET>')
      self.handles['puMet'] = AutoHandle('puMet','std::vector<reco::PFMET>')
      self.handles['tkMet'] = AutoHandle('tkMet','std::vector<reco::PFMET>')
      
      self.handles['muons'] = AutoHandle(
            'cmgMuonSel',
            'std::vector<cmg::Muon>'
            )
      self.handles['electrons'] = AutoHandle(
            'cmgElectronSel',
            'std::vector<cmg::Electron>'
            )
      self.mchandles['genParticles'] = AutoHandle( 'genParticlesPruned',
            'std::vector<reco::GenParticle>' )
    
      self.handles['vertices'] =  AutoHandle(
          'offlinePrimaryVertices',
          'std::vector<reco::Vertex>'
          )
      self.mchandles['pusi'] =  AutoHandle(
          'addPileupInfo',
          'std::vector<PileupSummaryInfo>' 
          ) 
      self.mchandles['generator'] = AutoHandle(
          'generator','GenEventInfoProduct' 
          )
    
    def declareVariables(self):
      tr = self.tree

      var(tr, 'scalePDF')
      var(tr, 'parton1_pdgId')
      var(tr, 'parton1_x')
      var(tr, 'parton2_pdgId')
      var(tr, 'parton2_x')

      var( tr, 'run', int)
      var( tr, 'lumi', int)
      var( tr, 'evt', int)
      var( tr, 'nvtx', int)
      var( tr, 'njets', int)
      var( tr, 'npu', int)
      var( tr, 'evtHasGoodVtx', int)
      var( tr, 'Vtx_ndof', int)
      var( tr, 'firstVtxIsGood', int)
      var( tr, 'evtHasTrg', int)
      var( tr, 'evtWSel', int)
      
      var( tr, 'nMuons', int)
      var( tr, 'nTrgMuons', int)
      var( tr, 'noTrgMuonsLeadingPt', int)

      bookMET( tr, 'pfmet')
      bookMET( tr, 'nopumet')
      bookMET( tr, 'tkmet')
      bookMET( tr, 'pumet')
      bookMET( tr, 'pucmet')
      bookMET2( tr, 'pfMetForRegression')
      bookMET2( tr, 'pfmetraw')


      bookW( tr, 'W')
      var( tr, 'W_mt')
      bookW( tr, 'WGen')
      var( tr, 'WGen_rap')
      var( tr, 'WGen_m')
      var( tr, 'WGen_mt')
      var( tr, 'u1')
      var( tr, 'u2')

      bookParticle(tr, 'Mu')
      var(tr, 'Mu_dxy')
      var(tr, 'Mu_dz')
      var(tr, 'MuIsTightAndIso', int)
      var(tr, 'MuIsTight', int)
      var(tr, 'MuRelIso')
      var(tr, 'pt_vis')
      var(tr, 'phi_vis')
      bookParticle(tr, 'MuGen')
      bookParticle(tr, 'MuGenStatus1')
      var(tr, 'MuDRGenP')
      
      bookParticle(tr, 'NuGen')
      var(tr, 'FSRWeight')
      
      bookMuonCovMatrix(tr,'Mu' )
      
      bookJet(tr, 'Jet_leading')

       
    def process(self, iEvent, event):
        
        self.readCollections( iEvent )
        tr = self.tree
        tr.reset()
          
        if (event.savegenpW and self.cfg_comp.isMC):
          
          fill(tr, 'FSRWeight',event.fsrWeight)
          fillW(tr,'WGen', event.genW[0].p4())
          fill(tr, 'WGen_rap', event.genW[0].p4().Rapidity())
          fill(tr, 'WGen_m', event.genW[0].p4().M())
          fill(tr, 'WGen_mt', event.genW_mt)
          fillParticle(tr, 'MuGen',event.genMu[0])
          fillParticle(tr, 'MuGenStatus1', event.genMuStatus1[0])      
          fill(tr, 'MuDRGenP',event.muGenDeltaRgenP)
          fillParticle(tr, 'NuGen', event.genNu[0])

        if event.WGoodEvent == True :
                      
          fillW( tr, 'W',event.W4V)
          fill(tr, 'W_mt', event.W4V_mt)
          fill(tr, 'u1', event.u1)
          fill(tr, 'u2', event.u2)

          fillParticle(tr, 'Mu', event.selMuons[0])
          if ( event.selMuons[0].isGlobalMuon() or event.selMuons[0].isTrackerMuon() ) and event.passedVertexAnalyzer:
            fill(tr, 'Mu_dxy', math.fabs(event.selMuons[0].dxy()))
            fill(tr, 'Mu_dz', math.fabs(event.selMuons[0].dz()))
          fill(tr, 'MuIsTightAndIso', event.selMuonIsTightAndIso)
          fill(tr, 'MuIsTight', event.selMuonIsTight)
          fill(tr, 'MuRelIso', event.selMuons[0].relIso(0.5))
          fill(tr, 'pt_vis', event.selMuons[0].pt())
          fill(tr, 'phi_vis', event.selMuons[0].phi())
          
          fillMuonCovMatrix( tr,'Mu' ,event.covMatrixMuon ,event)
          
        if (event.savegenpW and self.cfg_comp.isMC) or event.WGoodEvent:
          fill( tr, 'run', event.run) 
          fill( tr, 'lumi',event.lumi)
          fill( tr, 'evt', event.eventId)
          fill( tr, 'nvtx', len(self.handles['vertices'].product()))          
          fill( tr, 'njets', len(event.selJets))
          if (self.cfg_comp.isMC) :
            event.pileUpInfo = map( PileUpSummaryInfo,
                                    self.mchandles['pusi'].product() )
            for puInfo in event.pileUpInfo:
              if puInfo.getBunchCrossing()==0:
                fill( tr, 'npu', puInfo.nTrueInteractions())
                # print 'puInfo.nTrueInteractions()= ',puInfo.nTrueInteractions()
              # else:
                # print 'NO INFO FOR puInfo.getBunchCrossing()==0 !!!!'
          if (self.cfg_comp.isMC) :
              event.generator = self.mchandles['generator'].product()
              # print 'WTreeProducer.py: ',event.generator.pdf().scalePDF,' ',event.generator.pdf().id.first,' ',event.generator.pdf().x.first,' ',event.generator.pdf().id.second,' ',event.generator.pdf().x.second
              fill(tr, 'scalePDF',float(event.generator.pdf().scalePDF))
              fill(tr, 'parton1_pdgId',float(event.generator.pdf().id.first))
              fill(tr, 'parton1_x',float(event.generator.pdf().x.first))
              fill(tr, 'parton2_pdgId',float(event.generator.pdf().id.second))
              fill(tr, 'parton2_x',float(event.generator.pdf().x.second))
            
          fill( tr, 'nMuons', event.nMuons)
          fill( tr, 'nTrgMuons', len(event.selMuons))
          if len(event.NoTriggeredMuonsLeadingPt) > 0 :
            fill( tr, 'noTrgMuonsLeadingPt', event.NoTriggeredMuonsLeadingPt[0].pt())
          fill( tr, 'evtHasGoodVtx', event.passedVertexAnalyzer)
          fill( tr, 'Vtx_ndof', event.goodVertices[0].ndof())
          # fill( tr, 'firstVtxIsGood', event.firstVtxIsGoodVertices) # REQUIRES DEFINITION IN CMGTools/RootTools/python/analyzers/VertexAnalyzer.py
          fill( tr, 'evtHasTrg', event.passedTriggerAnalyzer)
          fill( tr, 'evtWSel', event.WGoodEvent)
          fillMET(tr, 'pfmet', event.pfmet)
          event.pfmetraw = self.handles['pfMetraw'].product()[0]
          event.pfMetSignificance = self.handles['pfMetSignificance'].product()
          event.nopumet = self.handles['nopuMet'].product()[0]
          event.pucmet = self.handles['pucMet'].product()[0]
          event.pfMetForRegression = self.handles['pfMetForRegression'].product()[0]
          event.pumet = self.handles['puMet'].product()[0]
          event.tkmet = self.handles['tkMet'].product()[0]
          fillMET2(tr, 'nopumet', event.nopumet.p4())
          fillMET2(tr, 'tkmet', event.tkmet.p4())
          fillMET2(tr, 'pucmet', event.pucmet.p4())
          fillMET2(tr, 'pumet', event.pumet.p4())
          fillMET2(tr, 'pfMetForRegression', event.pfMetForRegression.p4())
          fillMET2(tr, 'pfmetraw', event.pfmetraw.p4())

          if len(event.selJets)>0:
            fillJet(tr, 'Jet_leading', event.selJets[0])
            
          self.tree.tree.Fill()
       
