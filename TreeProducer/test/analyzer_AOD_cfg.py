import FWCore.ParameterSet.Config as cms

### CMSSW command line parameter parser
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')

## data or MC options
options.register(
	'isData',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to indicate data or MC')

options.register(
	'maxEvts',-1,VarParsing.multiplicity.singleton,VarParsing.varType.int,
	'flag to indicate max events to process')


process = cms.Process("HEXAQTREE")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/EventContent/EventContent_cff')

from Configuration.AlCa.GlobalTag import GlobalTag

if(options.isData==True):
    process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v8', '')
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvts))
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(500)


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
     'file:///user/gflouris/Analysis/SIMPS/SignalProduction/testPG/CMSSW_8_0_21/src/NoPile/SUS-RunIISummer16DR80Premix-00068.root'
    )
)

# Adaptive vertex finder
process.load('RecoVertex/AdaptiveVertexFinder/inclusiveVertexFinder_cfi')

# Analyzer
process.load("HexaAnalysis.TreeProducer.Analyzer_cfi")
process.analyzer.isData = cms.untracked.bool(False)
process.p = cms.Path(process.inclusiveVertexFinder + process.analyzer)


# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('analysis.root')
)
