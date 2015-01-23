import FWCore.ParameterSet.Config as cms
import math
import sys
import random
NEvents = 10 #number of events to simulate
det1 = 203.8   #position of first tracker detector
det2 = 215.
tof  = 215.5   #position of time of flight detector
trklen = 10 # default tracker length
phys_process = "GG"
Qtype = 0
OutputFile = "Test2_GG.root"
InputFile = ""
min_mass = 300
max_mass = 2000
higgs_decay = 0
higgs_mass  = 125.7
top_mass = 173.34
phi_min = -math.pi
phi_max =  math.pi
t_min   = 0.
t_max   = 2.
xi_min  = 0.
xi_max  = 0.2
pileup  = False
hit_smear= True
Ang_smear= True
E_smear = True
vtx_smear= True
run_with_CR = True
beta_star = 0.60
det1xoffset = 0.
det2xoffset = 0.

ecms = 13000.

for arg in sys.argv[2:]:
    if arg=="no-pileup":
         pileup = False 
         continue
    if arg=="no-smearing":
         vtx_smear = False
         Ang_smear = False
         E_smear   = False
         hit_smear = False
    if arg=="no-vertexsmear":
         vtx_smear = False
         continue
    if arg=="no-anglesmear":
         Ang_smear = False
         continue
    if arg=="no-energysmear":
         E_smear = False
         continue
    if arg=="no-beamsmear":
         Ang_smear = False
         E_smear   = False
         continue
    if arg=="no-crossingangle":
         run_with_CR = False
         continue
    opt = arg[0:arg.find("=")]

    if opt.lower()=="nevents":
         NEvents = int(arg[arg.find("=")+1:])
    elif opt=="output":
         OutputFile = arg[arg.find("=")+1:]
    elif opt=="input":
         InputFile = arg[arg.find("=")+1:]
    elif opt=="det1":
         det1=float(arg[arg.find("=")+1:])
    elif opt=="det2":
         det2=float(arg[arg.find("=")+1:])
    elif opt=="tof":
         tof=float(arg[arg.find("=")+1:])
    elif opt=="min_mass":
         min_mass=float(arg[arg.find("=")+1:])
    elif opt=="max_mass":
         max_mass=float(arg[arg.find("=")+1:])
    elif opt=="hdecay":
         hdecay=arg[arg.find("=")+1:]
         if hdecay=="b":
            higgs_decay = 5
         elif hdecay=="c":
            higgs_decay = 4
         elif hdecay=="W":
            higgs_decay = 24
         elif hdecay=="G":
            higgs_decay = 21
         elif hdecay=="tau":
            higgs_decay = 15
         elif hdecay=="Z":
            higgs_decay = 23
         elif hdecay=="g":
            higgs_decay = 22
         else:
            print "Unknown decay mode"
            sys.exit(2)
    elif opt=="higgs_mass":
         higgs_mass = float(arg[arg.find("=")+1:])  
    elif opt=="min_t":
         t_min = float(arg[arg.find("=")+1:])
    elif opt=="max_t":
         t_max = float(arg[arg.find("=")+1:])
    elif opt=="min_xi":
         xi_min = float(arg[arg.find("=")+1:])
    elif opt=="max_xi":
         xi_max = float(arg[arg.find("=")+1:])
    elif opt=="ecms":
         ecms = float(arg[arg.find("=")+1:])
    elif opt=="Q":
         if (arg[arg.find("=")+1:]=="top"):
            Qtype=6
         if (arg[arg.find("=")+1:]=="bottom"):
            Qtype=5
    elif opt=="beta_star":
         beta_star=float(arg[arg.find("=")+1:])
         if beta_star!=0.55 and beta_star!=0.60 and beta_star!=90:
            print "Unknown BetaStar (0.55, 0.60 or 90)"
            sys.exit(2)
    elif opt=="det1Offset":
         det1xoffset = float(arg[arg.find("=")+1:])
    elif opt=="det2Offset":
         det2xoffset = float(arg[arg.find("=")+1:])
    else:
         print opt
         print "Usage: cmsRun ",sys.argv[1],"process=[GG]","NEvent=# events","det1=det positon","tof=ToF position"

print det1xoffset,det2xoffset

if beta_star==90: run_with_CR=False

if det2>0:
   trklen=det2-det1
else:
   det2=det1+trklen

if run_with_CR:
   useCR     = True
   beamfile1 = '/afs/cern.ch/work/p/polme/public/PPS/FastSim/CMSSW_7_4_0_pre2/src/LHCB1_Beta0.60_6.5TeV_CR142.5_v6.503.tfs'
   beamfile2 = '/afs/cern.ch/work/p/polme/public/PPS/FastSim/CMSSW_7_4_0_pre2/src/LHCB2_Beta0.60_6.5TeV_CR142.5_v6.503.tfs'


print "Number of events : ",NEvents
print "Physics process  : ",phys_process
print "sqrt(s)          : ",ecms
print "1st Det. position: ",det1
print "2nd Det. position: ",det2
print "Tracker length   : ",trklen
print "ToF det. position: ",tof
print "Minimum mass     : ",min_mass
print "Maximum mass     : ",max_mass
print "Minimum t        : ",t_min
print "Maximum t        : ",t_max
print "Minimum xi       : ",xi_min
print "Maximum xi       : ",xi_max
print "Minimum phi      : ",phi_min
print "Maximum phi      : ",phi_max
print "Higgs mass       : ",higgs_mass
print "Higgs decay      : ",higgs_decay
print "Output File      : ",OutputFile
print "PileUp           : ",pileup
print "Angle smearing   : ",Ang_smear
print "Energy smearing  : ",E_smear
print "Vertex smearing  : ",vtx_smear
print "Use Xangle       : ",useCR
print "Quark type       : ",Qtype
print "BeamLine file 1  : ",beamfile1
print "BeamLine file 2  : ",beamfile2


process = cms.Process("PPS")

# import of standard configurations
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load("Geometry.PPSCommonData.cmsPPSGeometryXML_cfi")
process.load("SimG4Core.Application.g4SimHits_cfi")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')

# Pileup
if pileup:
   process.load('SimGeneral.MixingModule.mix_E14TeV_L10E33_BX2808_cfi')
   process.mix.input.fileNames = cms.untracked.vstring(
          'file:/afs/cern.ch/work/m/mundim/public/PPS/MinBias_TuneA2MB_13TeV_CMSSW_6_2_0-pythia8_GEN-SIM-001CB469-A91E-E311-9BFE-0025907FD24A.root',
          'file:/afs/cern.ch/work/m/mundim/public/PPS/MinBias_TuneA2MB_13TeV_CMSSW_6_2_0-pythia8_GEN-SIM-009CB248-A81C-E311-ACD8-00259073E4F0.root',
          'file:/afs/cern.ch/work/m/mundim/public/PPS/MinBias_TuneA2MB_13TeV_CMSSW_6_2_0-pythia8_GEN-SIM-009F81D5-B21C-E311-966C-BCAEC50971D0.root',
          'file:/afs/cern.ch/work/m/mundim/public/PPS/MinBias_TuneA2MB_13TeV_CMSSW_6_2_0-pythia8_GEN-SIM-00B5BB8C-A91E-E311-816A-782BCB1F5E6B.root',
          'file:/afs/cern.ch/work/m/mundim/public/PPS/MinBias_TuneA2MB_13TeV_CMSSW_6_2_0-pythia8_GEN-SIM-00B8F676-BA1C-E311-BA87-0019B9CABFB6.root',
          'file:/afs/cern.ch/work/m/mundim/public/PPS/MinBias_TuneA2MB_13TeV_CMSSW_6_2_0-pythia8_GEN-SIM-022A782D-A51C-E311-9856-80000048FE80.root',
          'file:/afs/cern.ch/work/m/mundim/public/PPS/MinBias_TuneA2MB_13TeV_CMSSW_6_2_0-pythia8_GEN-SIM-02A819D3-B21C-E311-A8D0-0025907DC9D6.root',
          'file:/afs/cern.ch/work/m/mundim/public/PPS/MinBias_TuneA2MB_13TeV_CMSSW_6_2_0-pythia8_GEN-SIM-04604218-B31C-E311-93A8-0017A4770030.root'
   )
   process.mix.input.nbPileupEvents.Lumi = cms.double(20) # pileup 50
   #process.mix.input.nbPileupEvents.Lumi = cms.double(10) # pileup 25
   process.mix.minBunch = cms.int32(0)
   process.mix.maxBunch = cms.int32(0)
else:
   process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.SimulationRandomNumberGeneratorSeeds_cff")
process.RandomNumberGeneratorService.g4SimHits.initialSeed  = cms.untracked.uint32(random.randint(0,900000000))
process.RandomNumberGeneratorService.VtxSmeared.initialSeed = cms.untracked.uint32(random.randint(0,900000000))
process.RandomNumberGeneratorService.generator.initialSeed  = cms.untracked.uint32(random.randint(0,900000000))
process.RandomNumberGeneratorService.mix.initialSeed  = cms.untracked.uint32(random.randint(0,900000000))

#
process.load('PhysicsTools.HepMCCandAlgos.genParticles_cfi')

# Vertex position. Taken from EventVertexGenerators/python/VtxSmearedParameters_cfi, the closest from the MB vertex position used
if vtx_smear:
   VertexX = 0.2440
   VertexY = 0.3929
   VertexZ = 0.4145
else:
   VertexX = 0.
   VertexY = 0.
   VertexZ = 0.

process.VtxSmeared.MeanX    = cms.double(VertexX)
process.VtxSmeared.MeanY    = cms.double(VertexY)
process.VtxSmeared.MeanZ    = cms.double(VertexZ)


# import of standard configurations
from GeneratorInterface.ExhumeInterface.ExhumeParameters_cfi import ExhumeParameters as ExhumeParametersRef
ExhumeParametersRef.HiggsMass = cms.double(higgs_mass)
ExhumeParametersRef.TopMass = cms.double(top_mass)

if phys_process=="GG":
   process.generator = cms.EDFilter("ExhumeGeneratorFilter",
       PythiaParameters = cms.PSet(
          parameterSets = cms.vstring()
       ),
       ExhumeParameters = ExhumeParametersRef,
       comEnergy = cms.double(ecms),
       pythiaHepMCVerbosity = cms.untracked.bool(False),
       maxEventsToPrint = cms.untracked.int32(2),
       pythiaPylistVerbosity = cms.untracked.int32(1),
       ExhumeProcess = cms.PSet(
           MassRangeLow = cms.double(min_mass),
           MassRangeHigh = cms.double(max_mass),
           ProcessType = cms.string(phys_process),
           ThetaMin = cms.double(0.30),
           HiggsDecay = cms.int32(higgs_decay),
           QuarkType  = cms.int32(Qtype)
       )
   )
   process.source = cms.Source("EmptySource")

configurationMetadata = cms.untracked.PSet(
        version = cms.untracked.string('\$Revision: 1.1 $'),
        name = cms.untracked.string('\$Source: /local/reps/CMSSW/CMSSW/Configuration/GenProduction/python/ThirteenTeV/QCD_Pt_15to3000_Tune4C_Flat_8TeV_pythia8_cff.py,v $'),
        annotation = cms.untracked.string('Summer2013 sample with PYTHIA8: QCD dijet production, pTHat = 15 .. 3000 GeV, weighted, Tune4C')
)
process.options = cms.untracked.PSet( SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(NEvents) )

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.398 $'),
    annotation = cms.untracked.string('reco nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)

ppssim_options = cms.PSet(
                         Verbosity = cms.untracked.int32(0),
                         Beam1File = cms.string(beamfile1),
                         Beam2File = cms.string(beamfile2),
                         Beam1Direction = cms.int32(1),
                         Beam2Direction = cms.int32(1),
                         SmearEnergy    = cms.bool(E_smear),
                         SmearAngle     = cms.bool(Ang_smear),
                         BeamEnergyRMS  = cms.double(1.11e-4),
                         BeamAngleRMS   = cms.double(30.03), # in mrad
                         ShowBeamLine   = cms.untracked.bool(False),
                         SimBeamProfile = cms.untracked.bool(False),
                         VtxMeanX       = cms.untracked.double(VertexX),
                         VtxMeanY       = cms.untracked.double(VertexY),
                         VtxMeanZ       = cms.untracked.double(VertexZ),
                         CollisionPoint = cms.string("IP5"),
                         TrackerPosition = cms.double(det1),
                         TrackerLength   = cms.double(trklen),
                         TrkDet1XOffset  = cms.double(det1xoffset),
                         TrkDet2XOffset  = cms.double(det2xoffset),
                         ToFPosition     = cms.double(tof),
                         TCL4Position    = cms.untracked.double(143.0),
                         TCL5Position    = cms.untracked.double(183.8),
                         TCL6Position    = cms.untracked.double(256.7),
                         SmearHit        = cms.bool(hit_smear),
                         HitSigmaX       = cms.double(10),
                         HitSigmaY       = cms.double(10),
                         HitSigmaZ       = cms.double(0),
                         TimeSigma       = cms.double(0.01), #in ns
                         PhiMin          = cms.double(phi_min),
                         PhiMax          = cms.double(phi_max),
                         CentralMass     = cms.double(125.7),
                         CentralMassErr  = cms.double(0.4),
                         KickersOFF      = cms.bool(True),
                         CrossAngleCorr  = cms.bool(useCR),
                         CrossingAngle   = cms.double(142.5), #in mrad
                         EtaMin          = cms.double(7.0),  # min eta to be tracked by HECTOR
                         MomentumMin     = cms.double(3.000), # min mom. to be tracked by HECTOR
                         TrackImpactParameterCut = cms.double(-1), # max. imp. par. for reco tracks
                         MinThetaXatDet1 = cms.double(-500.), #min. theta x at first tracker in urad
                         MaxThetaXatDet1 = cms.double(500.),   #max. theta x at first tracker in urad
                         MinThetaYatDet1 = cms.double(-500.), #min. theta y at first tracker in urad
                         MaxThetaYatDet1 = cms.double(500.), #max. theta y at first tracker in urad
                         DetectorClosestX= cms.double(-2.),  #min. distance to the beam
                         MaxXfromBeam    = cms.double(-25),  #max. x from beam for a hit
                         MaxYfromBeam    = cms.double(10.),  #max |y| from beam for a hit
                         FilterHitMap    = cms.bool(True)    #apply geometrical (X,Y) in the hits
                     )

if phys_process=="GG":
    process.ppssim = cms.EDProducer('PPSProducer',
                         ppssim_options,
                         genSource       = cms.InputTag("genProtonsPU") if pileup else  cms.InputTag("generator")  # for HepMC event -> no pileup events
                     )

# Output definition
if phys_process=="GG":
    print "QCD output"
    process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
            outputCommands = cms.untracked.vstring('drop *',
                             'keep PPSSpectrometer_*_*_*',
                             'keep recoTracksToOnerecoTracksAssociation_*_*_*',
                             'keep edmHepMCProduct_*_*_*',
                             'keep GenRunInfoProduct_*_*_*',
                             'keep GenEventInfoProduct_*_*_*',
                             'keep TrackingRecHitsOwned_*_*_*',
                             'keep CaloTowersSorted_*_*_*',
                             'keep CastorRecHitsSorted_*_*_*',
                             'keep edmTriggerResults_*_*_*',
                             'keep recoJetIDedmValueMap_*_*_*',
                             'keep recoBeamSpot_*_*_*',
                             'keep recoBeamHaloSummary_*_*_*',
                             'keep TrajectorysToOnerecoGsfTracksAssociation_*_*_*',
                             'keep PileupSummaryInfos_*_*_*',
                             'keep TrackCandidates_*_*_*',
                             'keep TrajectorySeeds_*_*_*',
                             'keep recoBasicJets_*_*_*',
                             'keep recoCaloJets_*_*_*',
                             'keep *_genParticles_*_*',
                             'keep *_genParticlesForJets_*_*',
                             'keep *_iterativeCone5GenJets_*_*',
                             'keep *_ak5GenJets_*_*',
                             'keep *_ak7GenJets_*_*',
                             'keep recoPFCandidateedmPtredmValueMap_*_*_*',
                             'keep recoPFCandidates_*_*_*',
                             'keep recoPFClusters_*_*_*',
                             'keep recoPFJets_*_*_*',
                             'keep recoTrackJets_*_*_*',
                             'keep recoTracks_*_*_*',
                             'keep recoTrackExtras_*_*_*',
                             'keep recoVertexs_*_*_*',
                             'keep recoVertex_*_*_*',
                             'keep recoVertexCompositeCandidates_*_*_*'
                             ),
            fileName = cms.untracked.string(OutputFile),
            dataset = cms.untracked.PSet(
                      filterName = cms.untracked.string(''),
                      dataTier = cms.untracked.string('')
            )
    )
# Other statements
if pileup:
   process.genParticlesPU = cms.EDProducer("GenParticleProducer",
                            saveBarCodes = cms.untracked.bool(True),
                            #src = cms.untracked.InputTag("mix","generator"),
                            mix = cms.string("mix"),
                            abortOnUnknownPDGCode = cms.untracked.bool(False),
                            useCrossingFrame = cms.untracked.bool(True)
        )
   process.genProtonsPU = cms.EDFilter("GenParticleSelector",
                          filter = cms.bool(False),
                          src = cms.InputTag("genParticlesPU"),
                          cut = cms.string('')
        )
   process.genProtonsPU.cut = 'status = 1 & pdgId == 2212 & abs(pz) >= %f' % ( 0.5*ecms/2.0)
   outputCommandsPU = [ 'keep *_genParticlesPU_*_*', 'keep *_genProtonsPU_*_*' ]
   process.AODSIMoutput.outputCommands.extend( outputCommandsPU )


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')

# Path and EndPath definitions
process.generation_step = cms.Path(process.generator)
process.genparticle_step= cms.Path(process.genParticles)
 
if phys_process=="GG":
   process.genparticle_step= cms.Path(process.genParticles*process.genParticlesForJets+process.kt4GenJets+process.kt6GenJets+process.iterativeCone5GenJets+process.ak5GenJets+process.ak8GenJets+process.genCandidatesForMET+process.genParticlesForMETAllVisible+process.genMetCalo+process.genMetCaloAndNonPrompt+process.genMetTrue+process.genMetIC5GenJets)

process.vtxsmearing_step = cms.Path(process.VtxSmeared)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.ppssim_step = cms.Path(process.ppssim)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)

if pileup:
   process.digitisation_step.replace(process.pdigi, process.pdigi * process.genParticlesPU * process.genProtonsPU)

if vtx_smear:
   process.generation_step.replace(process.generator,process.generator * process.VtxSmeared)

# Schedule definition
process.schedule = cms.Schedule(
                                   process.generation_step,
                                   process.genparticle_step,
                                   process.genfiltersummary_step,
                                   process.simulation_step,
                                   process.digitisation_step,
                                   process.L1simulation_step,
                                   process.digi2raw_step,
                                   process.raw2digi_step,
                                   process.reconstruction_step,
                                   process.ppssim_step,
                                   process.endjob_step,
                                   process.AODSIMoutput_step
   )
