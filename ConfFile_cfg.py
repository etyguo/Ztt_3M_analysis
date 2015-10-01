import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/scratch/osg/etyguo/Ztt_3M_sim/1.root',
				'file:/scratch/osg/etyguo/Ztt_3M_sim/2.root',
				'file:/scratch/osg/etyguo/Ztt_3M_sim/3.root',
				'file:/scratch/osg/etyguo/Ztt_3M_sim/4.root',
				'file:/scratch/osg/etyguo/Ztt_3M_sim/5.root',
				'file:/scratch/osg/etyguo/Ztt_3M_sim/6.root',
				'file:/scratch/osg/etyguo/Ztt_3M_sim/7.root',
				'file:/scratch/osg/etyguo/Ztt_3M_sim/8.root',
				'file:/scratch/osg/etyguo/Ztt_3M_sim/9.root',
				'file:/scratch/osg/etyguo/Ztt_3M_sim/10.root',
    		'file:/scratch/osg/etyguo/Ztt_3M_sim/191.root',
				'file:/scratch/osg/etyguo/Ztt_3M_sim/192.root',
				'file:/scratch/osg/etyguo/Ztt_3M_sim/193.root',
				'file:/scratch/osg/etyguo/Ztt_3M_sim/194.root',
				'file:/scratch/osg/etyguo/Ztt_3M_sim/195.root',
				'file:/scratch/osg/etyguo/Ztt_3M_sim/196.root',
				'file:/scratch/osg/etyguo/Ztt_3M_sim/197.root',
				'file:/scratch/osg/etyguo/Ztt_3M_sim/198.root',
				'file:/scratch/osg/etyguo/Ztt_3M_sim/199.root',
				'file:/scratch/osg/etyguo/Ztt_3M_sim/200.root'
		),
		duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.demo = cms.EDAnalyzer('gen'
)
process.TFileService = cms.Service("TFileService",
                                fileName = cms.string('histro.root'))

process.p = cms.Path(process.demo)
