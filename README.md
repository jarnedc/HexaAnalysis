# HexaQuark Analysis

## Check-out code
```
mkdir ./Analysis/HexaQuark
cd ./Analysis/HexaQuark
cmsrel CMSSW_9_2_8
cd CMSSW_9_2_8/src
git clone https://github.com/gflouris/HexaAnalysis
cmsenv
scram b -j2
```

## Produce Particle Gun MC
```
cd ../../..
mkdir ./ParticleGun
cd ./ParticleGun

cmsrel cmsrel CMSSW_7_1_20_patch3
cp ../HexaQuark/CMSSW_9_2_8/src/HexaAnalysis/TreeProducer/scripts/mc_ParticleGun/SUS-RunIISummer15GS-00146_GENSIM_cfg.py ./CMSSW_7_1_20_patch3/src/
cd CMSSW_7_1_20_patch3/src
cmsenv
## Check the content of the .py file below before running the below. It is in this file that you specify the energy/particle(http://pdg.lbl.gov/2013/reviews/rpp2012-rev-monte-carlo-numbering.pdf) to be generated etc. This script runs the gensim, this means that you simulate the formation of the particles (Pythia) and then how they will interact with the detector (Geant).
cmsRun SUS-RunIISummer15GS-00146_GENSIM_cfg.py
cd ../..
cmsrel CMSSW_8_0_21
cp ../HexaQuark/CMSSW_9_2_8/src/HexaAnalysis/TreeProducer/scripts/mc_ParticleGun/SUS-RunIISummer16DR80Premix-00068_* CMSSW_8_0_21/src/
cd CMSSW_8_0_21/src
cmsenv
##before running the below change the input file path in the file accordingly. This script uses the above generated files and adds pile up (these are the 0088E837-9985-E611-809C-0025905A60A0.root etc files), applies the L1 trigger, does digi to raw conversion and applies the HLT.
cmsRun SUS-RunIISummer16DR80Premix-00068_Step1_cfg.py
##This script does the RAW to digi converstion and does the RECO.
cmsRun SUS-RunIISummer16DR80Premix-00068_Step2_cfg.py
cd ../../..
```

## Produce ntuples
```
cd ./Analysis/HexaQuark/CMSSW_8_0_21/src/TreeProducer/test
cmsenv
cmsRun treeproducer_AOD_MC_cfg.py
```

## Template macro
```
cd ../macros/
root -l analysis.C+
```

## Documentation
### Adaptive Vertex Reco
- https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAdaptiveVertexReconstructor
- http://cds.cern.ch/record/1166320/files/NOTE2008_033.pdf
- http://iopscience.iop.org/article/10.1088/0954-3899/34/12/N01/pdf
