import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:myfile.root'
    ),
		duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.demo = cms.EDAnalyzer('gen'
)
process.TFileService = cms.Service("TFileService",
                                fileName = cms.string('histro.root'))

process.p = cms.Path(process.demo)
