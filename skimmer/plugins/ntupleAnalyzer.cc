#include <algorithm>
#include <cmath> 
#include <memory>
#include <random>
#include <vector>
#include <boost/format.hpp>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/InvariantMass.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
// #include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "TTree.h"
#include "TMath.h"

#include "NtupleContainer.hh"


class ntupleAnalyzer : public edm::EDAnalyzer { //<edm::one::WatchRuns, edm::one::SharedResources> {

public:
    explicit ntupleAnalyzer(const edm::ParameterSet&);
    ~ntupleAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions&);

protected:
    virtual void beginJob() override;
    virtual void beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) override;
    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void endJob() override;

private:
    bool getCollections(const edm::Event&);

    // pointers to TTrees created (and owned) by the EDM TFileService
    // Don't delete these at the end since we don't own them (and it's rude)
    TTree * recoT, * genT, * tauT;
    NtupleContainer nt;
    edm::Service<TFileService> fs;

    std::mt19937 m_random_generator; 

    edm::ESHandle<JetCorrectorParametersCollection> JetCorParCollHandle_;
    std::unique_ptr<JetCorrectionUncertainty> jecUnc_;

    bool isData;
    const std::string triggerProcessName_;

    template <typename T>
    struct edmDataT {
        const edm::EDGetTokenT<T> token;
        edm::Handle<T> handle;
    };

    // Reco info
    edmDataT<reco::VertexCollection> primaryVertex_;
    edmDataT<reco::PFJetCollection> recoJets_;
    edmDataT<reco::JetCorrector> jetCorrector_;
    edmDataT<reco::JetTagCollection> bTagProbb_;
    edmDataT<reco::JetTagCollection> bTagProbbb_;
    edmDataT<reco::PFMETCollection> recoPFMET_;
    edmDataT<reco::CaloMETCollection> recoCaloMET_;
    edmDataT<reco::MuonCollection> dsaRecoMu_;
    edmDataT<reco::MuonTimeExtraMap> dsaRecoMuTiming_;
    edmDataT<reco::MuonCollection> pfRecoMu_;
    edmDataT<reco::TrackCollection> muTrack1_;
    edmDataT<reco::TrackCollection> muTrack2_;
    edmDataT<reco::GsfElectronCollection> recoElectrons_;
    edmDataT<edm::ValueMap<float>> recoElectronID_;
    edmDataT<reco::PhotonCollection> recoPhotons_;
    edmDataT<edm::ValueMap<bool>> recoPhotonID_;
    edmDataT<reco::PFRecHitCollection> PFRecHitECAL_;
    edmDataT<reco::PFRecHitCollection> PFRecHitHBHE_;
    edmDataT<reco::PFRecHitCollection> PFRecHitHO_;
    edmDataT<reco::PFRecHitCollection> PFRecHitHF_;
    edmDataT<reco::PFRecHitCollection> PFRecHitPS_;
    edmDataT<reco::PFTauCollection> recoPFTaus_;

    // Trigger info
    edmDataT<edm::TriggerResults> trigResults_;
    edmDataT<trigger::TriggerEvent> trigEvent_;
    edmDataT<l1t::TauBxCollection> L1Taus_;
    edmDataT<l1t::JetBxCollection> L1Jets_;

    // Gen info
    edmDataT<double> rho_;
    edmDataT<std::vector<PileupSummaryInfo>> pileupInfos_;
    edmDataT<GenEventInfoProduct> genEventInfo_;
    edmDataT<GenLumiInfoHeader> genLumiHeader_;
    edmDataT<reco::GenParticleCollection> genParticles_;
    edmDataT<reco::GenJetCollection> genJets_;
    edmDataT<reco::GenMETCollection> genMET_;

    // Filter info
    edmDataT<bool> HBHENoiseFilterResultProducer_;
    edmDataT<bool> HBHEIsoNoiseFilterResultProducer_;
    edmDataT<bool> globalSuperTightHalo2016Filter_;
    edmDataT<bool> EcalDeadCellTriggerPrimitiveFilter_;
    edmDataT<bool> ecalBadCalibFilter_;
    edmDataT<bool> BadPFMuonFilter_;
    edmDataT<bool> muonBadTrackFilter_;
    edmDataT<int> primaryVertexFilter_;
    
    std::vector<std::string> triggerPathsWithoutVersionNum_;
    std::vector<std::string> triggerPathsWithVersionNum_;
    std::vector<bool> trigExist_;
    HLTConfigProvider hltConfig_;
};


ntupleAnalyzer::ntupleAnalyzer(const edm::ParameterSet& ps):
    isData(ps.getParameter<bool>("isData")),
    triggerProcessName_(ps.getParameter<std::string>("triggerProcessName")),

    // Reco info
    primaryVertex_{consumes<reco::VertexCollection>(ps.getParameter<edm::InputTag>("primaryVertex"))},
    recoJets_{consumes<reco::PFJetCollection>(ps.getParameter<edm::InputTag>("recoJet"))},
    jetCorrector_{consumes<reco::JetCorrector>(ps.getParameter<edm::InputTag>("jetCorrector"))},
    bTagProbb_{consumes<reco::JetTagCollection>(ps.getParameter<edm::InputTag>("bTagProbb"))},
    bTagProbbb_{consumes<reco::JetTagCollection>(ps.getParameter<edm::InputTag>("bTagProbbb"))},
    recoPFMET_{consumes<reco::PFMETCollection>(ps.getParameter<edm::InputTag>("PFMET"))},
    recoCaloMET_{consumes<reco::CaloMETCollection>(ps.getParameter<edm::InputTag>("caloMET"))},
    dsaRecoMu_{consumes<reco::MuonCollection>(ps.getParameter<edm::InputTag>("dsaRecoMu"))},
    dsaRecoMuTiming_{consumes<reco::MuonTimeExtraMap>(ps.getParameter<edm::InputTag>("dsaRecoMuTiming"))},
    pfRecoMu_{consumes<reco::MuonCollection>(ps.getParameter<edm::InputTag>("pfRecoMu"))},
    muTrack1_{consumes<reco::TrackCollection>(ps.getParameter<edm::InputTag>("muTrack1"))},
    muTrack2_{consumes<reco::TrackCollection>(ps.getParameter<edm::InputTag>("muTrack2"))},
    recoElectrons_{consumes<reco::GsfElectronCollection>(ps.getParameter<edm::InputTag>("recoElectron"))},
    recoElectronID_{consumes<edm::ValueMap<float>>(ps.getParameter<edm::InputTag>("electronID"))},
    recoPhotons_{consumes<reco::PhotonCollection>(ps.getParameter<edm::InputTag>("recoPhoton"))},
    recoPhotonID_{consumes<edm::ValueMap<bool>>(ps.getParameter<edm::InputTag>("photonID"))},
    PFRecHitECAL_{consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("PFRecHitECAL"))},
    PFRecHitHBHE_{consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("PFRecHitHBHE"))},
    PFRecHitHO_{consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("PFRecHitHO"))},
    PFRecHitHF_{consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("PFRecHitHF"))},
    PFRecHitPS_{consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("PFRecHitPS"))},
    recoPFTaus_{consumes<reco::PFTauCollection>(ps.getParameter<edm::InputTag>("PFTau"))},

    // Trigger info
    trigResults_{consumes<edm::TriggerResults>(ps.getParameter<edm::InputTag>("trigResult"))},
    trigEvent_{consumes<trigger::TriggerEvent>(ps.getParameter<edm::InputTag>("trigEvent"))},
    L1Taus_{consumes<l1t::TauBxCollection>(ps.getParameter<edm::InputTag>("L1Tau"))},
    L1Jets_{consumes<l1t::JetBxCollection>(ps.getParameter<edm::InputTag>("L1Jet"))},

    // Gen info
    rho_{consumes<double>(ps.getParameter<edm::InputTag>("rho"))},
    pileupInfos_{consumes<std::vector<PileupSummaryInfo>>(ps.getParameter<edm::InputTag>("pileups"))},
    genEventInfo_{consumes<GenEventInfoProduct>(ps.getParameter<edm::InputTag>("genEvt"))},
    genLumiHeader_{consumes<GenLumiInfoHeader,edm::InLumi>(edm::InputTag("generator",""))},
    genParticles_{consumes<reco::GenParticleCollection>(ps.getParameter<edm::InputTag>("genParticle"))},
    genJets_{consumes<reco::GenJetCollection>(ps.getParameter<edm::InputTag>("genJet"))},
    genMET_{consumes<reco::GenMETCollection>(ps.getParameter<edm::InputTag>("genMET"))},

    // Filter info
    HBHENoiseFilterResultProducer_{consumes<bool>(ps.getParameter<edm::InputTag>("HBHENoiseFilter"))},
    HBHEIsoNoiseFilterResultProducer_{consumes<bool>(ps.getParameter<edm::InputTag>("HBHEIsoNoiseFilter"))},
    globalSuperTightHalo2016Filter_{consumes<bool>(ps.getParameter<edm::InputTag>("globalSuperTightHaloFilter"))},
    EcalDeadCellTriggerPrimitiveFilter_{consumes<bool>(ps.getParameter<edm::InputTag>("EcalDeadCellTrgPrimitFilter"))},
    ecalBadCalibFilter_{consumes<bool>(ps.getParameter<edm::InputTag>("ecalBadCalibFilter"))},
    BadPFMuonFilter_{consumes<bool>(ps.getParameter<edm::InputTag>("BadPFMuonFilter"))},
    muonBadTrackFilter_{consumes<bool>(ps.getParameter<edm::InputTag>("muonBadTrackFilter"))},
    primaryVertexFilter_{consumes<int>(ps.getParameter<edm::InputTag>("primaryVertexFilter"))}
{
    // usesResource("TFileService");
    m_random_generator = std::mt19937(37428479);
}

ntupleAnalyzer::~ntupleAnalyzer() = default;

void ntupleAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    using namespace edm;

    // Specifies some of the default collections to be included in the generated .cfi 
    // Only specify tags with reasonable defaults -- check .cfg for others

    edm::ParameterSetDescription desc;
    desc.add<bool>("isData", 0);
    desc.add<std::string>("triggerProcessName", "HLT");

    desc.add<InputTag>("bTagProbb", InputTag("pfDeepCSVJetTags:probb")); 
    desc.add<InputTag>("bTagProbbb", InputTag("pfDeepCSVJetTags:probbb"));
    desc.add<InputTag>("dsaRecoMu", InputTag("muonsFromdSA"));
    desc.add<InputTag>("dsaRecoMuTiming", InputTag("muontimingFromdSA","combined"));
    desc.add<InputTag>("pfRecoMu", InputTag("muons"));
    desc.add<InputTag>("muTrack1", InputTag("displacedStandAloneMuons"));
    desc.add<InputTag>("muTrack2", InputTag("globalMuons"));
    desc.add<InputTag>("primaryVertexFilter", InputTag("myPrimaryVertexFilter"));
    desc.add<InputTag>("primaryVertex", InputTag("offlinePrimaryVertices"));
    desc.add<InputTag>("genParticle", InputTag("genParticles"));
    desc.add<InputTag>("genJet", InputTag("ak4GenJets"));
    desc.add<InputTag>("genMET", InputTag("genMetTrue"));
    desc.add<InputTag>("PFMET", InputTag("pfMetT0rtT1Txy"));
    desc.add<InputTag>("caloMET", InputTag("caloMet"));
    desc.add<InputTag>("recoJet", InputTag("ak4PFJetsCHS"));
    desc.add<InputTag>("trigResult", InputTag("TriggerResults", "", "HLT"));
    desc.add<InputTag>("trigEvent", InputTag("hltTriggerSummaryAOD", "", "HLT"));
    desc.add<InputTag>("recoElectron", InputTag("gedGsfElectrons"));
    desc.add<InputTag>("electronID"); // no default (depends on year), should throw if none provided in .cfg
    desc.add<InputTag>("recoPhoton", InputTag("gedPhotons"));
    desc.add<InputTag>("photonID"); // no default (depends on year), should throw if none provided in .cfg
    desc.add<InputTag>("pileups", InputTag("addPileupInfo"));
    desc.add<InputTag>("genEvt", InputTag("generator"));
    desc.add<InputTag>("HBHENoiseFilter", InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"));
    desc.add<InputTag>("HBHEIsoNoiseFilter", InputTag("HBHENoiseFilterResultProducer","HBHEIsoNoiseFilterResult"));
    desc.add<InputTag>("globalSuperTightHaloFilter", InputTag("globalSuperTightHalo2016Filter"));
    desc.add<InputTag>("EcalDeadCellTrgPrimitFilter", InputTag("EcalDeadCellTriggerPrimitiveFilter"));
    desc.add<InputTag>("ecalBadCalibFilter", InputTag("ecalBadCalibFilter"));
    desc.add<InputTag>("BadPFMuonFilter", InputTag("BadPFMuonFilter"));
    desc.add<InputTag>("muonBadTrackFilter", InputTag("muonBadTrackFilter"));
    desc.add<InputTag>("jetCorrector"); // no default (depends on year), should throw if none provided in .cfg
    desc.add<InputTag>("PFRecHitECAL", InputTag("particleFlowRecHitECAL", "Cleaned", "RECO"));
    desc.add<InputTag>("PFRecHitHBHE", InputTag("particleFlowRecHitHBHE", "Cleaned", "RECO"));
    desc.add<InputTag>("PFRecHitHO", InputTag("particleFlowRecHitHO", "Cleaned", "RECO"));
    desc.add<InputTag>("PFRecHitHF", InputTag("particleFlowRecHitHF", "Cleaned", "RECO"));
    desc.add<InputTag>("PFRecHitPS", InputTag("particleFlowRecHitPS", "Cleaned", "RECO"));
    desc.add<InputTag>("L1Tau", InputTag("caloStage2Digis", "Tau", "RECO"));
    desc.add<InputTag>("L1Jet", InputTag("caloStage2Digis", "Jet", "RECO"));
    desc.add<InputTag>("PFTau", InputTag("hpsPFTauProducer", "", "RECO"));
    desc.add<InputTag>("rho", InputTag("fixedGridRhoFastjetAll"));
    
    descriptions.add("ntupleAnalyzer", desc);
}

void ntupleAnalyzer::beginJob()
{
    recoT = fs->make<TTree>("recoT", "recoT");
    nt.SetRecoTree(recoT);
    if (!isData) {
        genT = fs->make<TTree>("genT", "genT");
        nt.SetGenTree(genT);

        tauT = fs->make<TTree>("tauT", "tauT");
        nt.SetTauTree(tauT);
    }
    nt.CreateTreeBranches();
}

void ntupleAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
  if (!isData) {
    iLumi.getByToken(genLumiHeader_.token, genLumiHeader_.handle);
  }

  // Fill lhe comment lines with SUSY model parameter information
  nt.lheComments = "";
  if (genLumiHeader_.handle.isValid()) {
    nt.lheComments = genLumiHeader_.handle->configDescription();
  }
  std::cout << "lheComments: " << nt.lheComments << std::endl;

}


void ntupleAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    using namespace edm;

    bool changed = true;
    if (hltConfig_.init(iRun, iSetup, triggerProcessName_, changed)) {
        if (changed) {
            LogInfo("HLTConfig") << "ntupleAnalyzer::beginRun: " << "hltConfig init for Run" << iRun.run();
            hltConfig_.dump("ProcessName");
            hltConfig_.dump("GlobalTag");
            hltConfig_.dump("TableName");
        }
    } 
    else {
        LogError("HLTConfig") << "ntupleAnalyzer::beginRun: config extraction failure with triggerProcessName -> " << triggerProcessName_;
        return;
    }

    // Add trigger paths if they exist
    triggerPathsWithoutVersionNum_.clear();
    triggerPathsWithVersionNum_.clear();
    trigExist_.clear();

    triggerPathsWithoutVersionNum_ = {
        "HLT_DoubleEle25_CaloIdL_MW",
        "HLT_DoubleEle27_CaloIdL_MW",
        "HLT_DoubleEle33_CaloIdL_MW",
        "HLT_DoubleEle24_eta2p1_WPTight_Gsf",
        "HLT_Ele27_Ele37_CaloIdL_MW",
        "HLT_DoubleMu3_Trk_Tau3mu",
        "HLT_DoubleMu3_TkMu_DsTau3Mu",
        "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1",
        "HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1",
        "HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1",
        "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1",
        "HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1",
        "HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1",
        "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1",
        "HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1",
        "HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1",
        "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1",
        "HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1",
        "HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1",
        "HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1",
        "HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS30_Trk1_eta2p1_Reg_CrossL1",
        "HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1",
        "HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1",
        "HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1",
        "HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15",
        "HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1",
        "HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15",
        "HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1",
        "HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass",
        "HLT_PFHT350MinPFJet15",
        "HLT_Photon60_R9Id90_CaloIdL_IsoL",
        "HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL",
        "HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr",
        "HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr",
        "HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1",
        "HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1",
        "HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1",
        "HLT_DoubleMediumChargedIsoPFTauHPS30_L1MaxMass_Trk1_eta2p1_Reg",
        "HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg",
        "HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg",
        "HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg",
        "HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg",
        "HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg",
        "HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg",
        "HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg",
        "HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg",
        "HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1",
        "HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1",
        "HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1"
    };
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1");
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMET120_PFMHT120_IDTight");
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMET130_PFMHT130_IDTight");
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMET140_PFMHT140_IDTight");
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight");
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight");
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight");
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMET120_PFMHT120_IDTight_PFHT60"); // 2017+2018
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60"); // 2017+2018
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMET200_HBHECleaned"); // 2017+2018
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMET200_HBHE_BeamHaloCleaned"); // 2017+2018
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned"); // 2017+2018
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_PFHT500_PFMET100_PFMHT100_IDTight"); // 2017+2018
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_PFHT700_PFMET85_PFMHT85_IDTight"); // 2017+2018
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_PFHT800_PFMET75_PFMHT75_IDTight"); // 2017+2018
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMET170_HBHECleaned"); // 2016
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_PFMET300"); // 2016
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_MET200"); // 2016
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_PFHT300_PFMET110"); // 2016
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_IsoMu27"); // For MET trigger eff. studies in data
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight"); // Alternative triggers
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_DoubleMu3_DCA_PFMET50_PFMHT60"); // Alternative triggers
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_DoubleMu3_DZ_PFMET50_PFMHT60");  // Alternative triggers
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_DoubleMu3_DZ_PFMET70_PFMHT70");  // Alternative triggers
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_DoubleMu3_DZ_PFMET90_PFMHT90");  // Alternative triggers
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_L2Mu10_NoVertex_NoBPTX");    // For dSA eff. studies in data
    // triggerPathsWithoutVersionNum_.emplace_back("HLT_L2Mu10_NoVertex_NoBPTX3BX"); // For dSA eff. studies in data
    
    const std::vector<std::string>& pathNames = hltConfig_.triggerNames();
    for (auto trigPathNoVersion : triggerPathsWithoutVersionNum_) {
        auto matchedPaths(hltConfig_.restoreVersion(pathNames, trigPathNoVersion));
        if (matchedPaths.size() == 0) {
            LogWarning("TriggerNotFound") << "Could not find matched full trigger path with --> " << trigPathNoVersion;
            triggerPathsWithVersionNum_.push_back("None");
            trigExist_.push_back(false);
        }
        else {
            trigExist_.push_back(true);
            triggerPathsWithVersionNum_.push_back(matchedPaths[0]);
            if (hltConfig_.triggerIndex(matchedPaths[0]) >= hltConfig_.size()) {
                LogError("TriggerError") << "Cannot find trigger path --> " << matchedPaths[0];
                return;
            }
        }
    }

    // JEC Uncertainty object from ESSetup
    iSetup.get<JetCorrectionsRecord>().get("AK4PFchs", JetCorParCollHandle_); 
    jecUnc_ = make_unique<JetCorrectionUncertainty>((*JetCorParCollHandle_)["Uncertainty"]);
    if (!jecUnc_)
        LogError("JECUncertainty") << "ntupleAnalyzer::beginRun: failed to get jecUnc_ object!";
}


bool ntupleAnalyzer::getCollections(const edm::Event& iEvent) {
    using namespace edm;
    char error_msg[] = "ntupleAnalyzer::GetCollections: Error in getting product %s from Event!";

    bool ret = true;
    auto getHandle = [&]<typename T>(edmDataT<T> &data, std::string name) {
        iEvent.getByToken(data.token, data.handle);
        if (!data.handle.isValid()) {
            LogError("HandleError") << boost::str(boost::format(error_msg) % name);
            ret = false;
        }
    };

    getHandle(bTagProbb_, "bTagProbb");
    getHandle(bTagProbbb_, "bTagProbbb");
    getHandle(dsaRecoMu_, "dsaRecoMu");
    getHandle(dsaRecoMuTiming_, "dsaRecoMuTiming");
    getHandle(pfRecoMu_, "pfRecoMu");
    getHandle(muTrack1_, "muTrack1");
    getHandle(muTrack2_, "muTrack2");
    getHandle(primaryVertex_, "primaryVertex");
    getHandle(recoPFMET_, "PFMET");
    getHandle(recoCaloMET_, "CaloMET");
    getHandle(recoJets_, "recoJet");
    getHandle(trigResults_, "trigResults");
    getHandle(trigEvent_, "trigEvent");
    getHandle(HBHENoiseFilterResultProducer_, "HBHENoiseFilter");
    getHandle(HBHEIsoNoiseFilterResultProducer_, "HBHEIsoNoiseFilter");
    getHandle(primaryVertexFilter_, "primaryVertexFilter");
    getHandle(globalSuperTightHalo2016Filter_, "globalSuperTightHalo2016");
    getHandle(EcalDeadCellTriggerPrimitiveFilter_, "EcalDeadCellTriggerPrimitive");
    getHandle(ecalBadCalibFilter_, "ecalBadCalibFilter");
    getHandle(BadPFMuonFilter_, "BadPFMuonFilter");
    getHandle(muonBadTrackFilter_, "muonBadTrackFilter");
    getHandle(jetCorrector_, "jetCorrector");
    getHandle(recoElectrons_, "recoElectron");
    getHandle(recoElectronID_, "recoElectronID");
    getHandle(recoPhotons_, "recoPhoton");
    getHandle(recoPhotonID_, "recoPhotonID");
    getHandle(PFRecHitECAL_, "PFRecHitECAL");
    getHandle(PFRecHitHBHE_, "PFRecHitHBH");
    getHandle(PFRecHitHO_, "PFRecHitHO");
    getHandle(PFRecHitHF_, "PFRecHitHF");
    getHandle(PFRecHitPS_, "PFRecHitPS");
    getHandle(L1Taus_, "L1Tau");
    getHandle(recoPFTaus_, "PFTau");
    getHandle(L1Jets_, "L1Jets");
    getHandle(rho_, "rho");
    if (!isData) {
        getHandle(genEventInfo_, "genEventInfo");
        getHandle(genParticles_, "genParticle");
        getHandle(genJets_, "genJet");
        getHandle(genMET_, "genMET");
        getHandle(pileupInfos_, "pileupInfos");
    }
    return ret;
}


void ntupleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using std::cout, std::vector, std::endl;

    if (!getCollections(iEvent))
        return;


    // for (auto & pfrechit : *PFRecHitECAL_.handle)
    //     cout << "pf rec hit energy " << pfrechit.energy() << ", detId " << pfrechit.detId() << ", layer " << pfrechit.layer() << endl; //pfrechit.position().x() << ", " << pfrechit.position().y() << ", " << pfrechit.position().z() << endl;


    // Clear branches before filling
    nt.ClearTreeBranches();

    // Start filling

    // Event information
    nt.eventNum_ = iEvent.id().event();
    nt.lumiSec_ = iEvent.luminosityBlock();
    nt.runNum_ = iEvent.id().run();
    nt.npv_ = *primaryVertexFilter_.handle;

    // Primary vertex
    reco::Vertex pv = *primaryVertex_.handle->begin();
    nt.pvx_ = pv.x();
    nt.pvy_ = pv.y();
    nt.pvz_ = pv.z();

    // MET filters
    if (!(*HBHENoiseFilterResultProducer_.handle)) // false means event is bad
        nt.METFiltersFailBits_ |= (1 << 0);
    if (!(*HBHEIsoNoiseFilterResultProducer_.handle)) // false means event is bad
        nt.METFiltersFailBits_ |= (1 << 1);
    if (!(*primaryVertexFilter_.handle)) // primaryVertexFilter == 0 means event is bad (number of vertices)
        nt.METFiltersFailBits_ |= (1 << 2);
    if (!(*globalSuperTightHalo2016Filter_.handle)) // false means event is bad
        nt.METFiltersFailBits_ |= (1 << 3);
    if (!(*EcalDeadCellTriggerPrimitiveFilter_.handle)) // false means event is bad
        nt.METFiltersFailBits_ |= (1 << 4);
    if (!(*ecalBadCalibFilter_.handle)) // false means event is bad
        nt.METFiltersFailBits_ |= (1 << 5);
    if (!(*BadPFMuonFilter_.handle)) // false means event is bad
        nt.METFiltersFailBits_ |= (1 << 6);
    if (!(*muonBadTrackFilter_.handle)) // false means event is bad
        nt.METFiltersFailBits_ |= (1 << 7);

    // get MET
    // assumes 0-th element of PFMET collection is largest pt (and only?) element
    // in other words, why is recoPFMET_.handle even a vector?
    reco::PFMETRef PFMETr(recoPFMET_.handle, 0);
    nt.recoPFMETPt_ = PFMETr->pt();
    nt.recoPFMETPhi_ = PFMETr->phi();
    nt.recoPFMETMuonEtFraction_ = PFMETr->MuonEtFraction();
    reco::CaloMETRef CaloMETr(recoCaloMET_.handle, 0);
    nt.recoCaloMETPt_ = CaloMETr->pt();
    nt.recoCaloMETPhi_ = CaloMETr->phi();
    
    // calculate MHT
    math::XYZTLorentzVector MHT;
    for (auto & jet : *recoJets_.handle) {
        if (jet.pt() < 30) continue;
        MHT += jet.p4();
    }
    nt.MHTPt_ = MHT.pt();
    
    // As we go along, calculate recoil with electrons, photons, and GM that passed ID
    double recoil_px = PFMETr->px();
    double recoil_py = PFMETr->py();

    // Add all electrons to ntuple, regardless of ID
    // for now "good" electron means only ID is passed
    // i.e. (IDmap % 2 == 1)
    nt.recoNElectron_ = recoElectrons_.handle->size();
    nt.recoNGoodElectron_ = 0;
    const edm::ValueMap<float> & eIDmap = *recoElectronID_.handle;
    for (size_t i = 0; i < recoElectrons_.handle->size(); i++) {
        reco::GsfElectronRef electronRef(recoElectrons_.handle, i);
        if ((int)eIDmap[electronRef] % 2 == 1 && electronRef->pt() > 10 && abs(electronRef->eta() < 2.5)) {
            nt.recoNGoodElectron_++;
            recoil_px += electronRef->px();
            recoil_py += electronRef->py();
        }
        nt.recoElectronPt_.push_back(electronRef->pt());
        nt.recoElectronEta_.push_back(electronRef->eta());
        nt.recoElectronPhi_.push_back(electronRef->phi());
        nt.recoElectronVxy_.push_back(electronRef->trackPositionAtVtx().rho());
        nt.recoElectronVz_.push_back(electronRef->trackPositionAtVtx().z());
        nt.recoElectronCharge_.push_back(electronRef->charge());
        nt.recoElectronIDResult_.push_back((int)eIDmap[electronRef]);
    }

    // Also add all photons to ntuple, regardless of ID
    // Photon ID only produces 1 or 0
    nt.recoNPhoton_ = recoPhotons_.handle->size();
    nt.recoNGoodPhoton_ = 0;
    const edm::ValueMap<bool> & phIDmap = *recoPhotonID_.handle;
    for (size_t i = 0; i < recoPhotons_.handle->size(); i++) {
        reco::PhotonRef photonRef(recoPhotons_.handle, i);
        if (phIDmap[photonRef] && photonRef->pt() > 15 && abs(photonRef->eta()) < 2.5) {
            nt.recoNGoodPhoton_++;
            recoil_px += photonRef->px();
            recoil_py += photonRef->py();
        }
        nt.recoPhotonPt_.push_back(photonRef->pt());
        nt.recoPhotonEta_.push_back(photonRef->eta());
        nt.recoPhotonPhi_.push_back(photonRef->phi());
        nt.recoPhotonIDResult_.push_back(phIDmap[photonRef]);
    }

    // Assign each trigger result to a different bit
    nt.fired_ = 0;
    for (size_t i = 0; i < triggerPathsWithVersionNum_.size(); i++) {
        if (trigExist_.at(i)) {
            std::string trigPath = triggerPathsWithVersionNum_[i];
            nt.fired_ |= (trigResults_.handle->accept(hltConfig_.triggerIndex(trigPath)) << i);
        }
        else {
            nt.fired_ |= (0 << i);
        }
    }

    // dSA muon object was created prior to running this module, from original dSA track collection
    // It contains extra timing information
    // Asking for standalone muon track recovers the original track object
    vector<reco::TrackRef> muTracks1{};
    vector<reco::MuonRef> muObjs1{};
    for (size_t i = 0; i < dsaRecoMu_.handle->size(); i++) {
        reco::MuonRef muon_i(dsaRecoMu_.handle, i);
        reco::TrackRef track_i = muon_i->standAloneMuon();
        if (track_i.isNonnull()) {
            muTracks1.emplace_back(track_i);
            muObjs1.emplace_back(muon_i);
        }
    }

    // Sort dSA muons (note that muon collection is *not* sorted by pT at first)
    nt.recoNDSA_ = muTrack1_.handle->size();

    sort(muObjs1.begin(), muObjs1.end(), [](const auto & l, const auto & r) {
            reco::TrackRef lt = l->standAloneMuon();
            reco::TrackRef rt = r->standAloneMuon();
            return lt->pt() > rt->pt();
            });

    sort(muTracks1.begin(), muTracks1.end(), [](const auto & l, const auto & r) {
            return l->pt() > r->pt();
            });

    // Sort global muons (note that muon collection is *not* sorted by pT at first)
    nt.recoNGM_ = muTrack2_.handle->size();

    vector<reco::TrackRef> muTracks2{};
    vector<reco::MuonRef> muObjs2{};
    for (size_t i = 0; i < pfRecoMu_.handle->size(); i++) {
        reco::MuonRef muon_i(pfRecoMu_.handle, i);
        reco::TrackRef track_i = muon_i->combinedMuon();
        if (track_i.isNonnull()) {
            muTracks2.emplace_back(track_i);
            muObjs2.emplace_back(muon_i);
        }
    }

    sort(muObjs2.begin(), muObjs2.end(), [](const auto & l, const auto & r) {
            reco::TrackRef lt = l->combinedMuon();
            reco::TrackRef rt = r->combinedMuon();
            return lt->pt() > rt->pt();
            });

    sort(muTracks2.begin(), muTracks2.end(), [](const auto & l, const auto & r) {
           return l->pt() > r->pt();
           });

    // Create separate collection for good quality dSA muons
    // DEPRECATED
    // Now store 4 leading muons in each collection for later processing in macros
    vector<int> muGoodTracksIdx{};
    for (size_t i = 0; i < muTracks1.size(); i++) {
        if (muTracks1[i]->hitPattern().muonStationsWithValidHits() < 2 ||
            muTracks1[i]->hitPattern().numberOfValidMuonHits() < 12 ||
            muTracks1[i]->normalizedChi2() > 10 ||
	       	muTracks1[i]->pt() < 5 ||
            abs(muTracks1[i]->eta()) > 2.4 ||
            muTracks1[i]->ptError()/muTracks1[i]->pt() > 1) {
                continue;
        }
        muGoodTracksIdx.push_back(i);
    }
    nt.recoNGoodDSA_ = muGoodTracksIdx.size();

    // Create separate collection for good quality global muons
    // TODO: this needs to change, GMs have tracker info as well
    // Cuts on muon chamber variables can be looser, while also
    // implementing cuts on tracker variables (the whole point
    // of having GMs as well)
    // DEPRECATED (only used for jet cross-cleaning)
    vector<int> muGoodTracksIdx2{};
    for (size_t i = 0; i < muTracks2.size(); i++) {
        if (muTracks2[i]->hitPattern().muonStationsWithValidHits() < 2 ||
            muTracks2[i]->hitPattern().numberOfValidMuonHits() < 12 ||
            muTracks2[i]->normalizedChi2() > 10 ||
	       	muTracks2[i]->pt() < 5 ||
            abs(muTracks2[i]->eta()) > 2.4 ||
            muTracks2[i]->ptError()/muTracks2[i]->pt() > 1) {
                continue;
        }
        muGoodTracksIdx2.push_back(i);
    }
    nt.recoNGoodGM_ = muGoodTracksIdx2.size();
    
    // Add leading 4 pt muons
    for (size_t i = 0; i < 4; i++) {
        if (i >= muTracks1.size()) break;
        reco::TrackRef mu_i = muTracks1[i];

        nt.recoDSAPt_.push_back(mu_i->pt());
        nt.recoDSAPtError_.push_back(mu_i->ptError());
        nt.recoDSAEta_.push_back(mu_i->eta());
        nt.recoDSAEtaError_.push_back(mu_i->etaError());
        nt.recoDSAPhi_.push_back(mu_i->phi());
        nt.recoDSAPhiError_.push_back(mu_i->phiError());
        nt.recoDSADxy_.push_back(mu_i->dxy(pv.position()));
        nt.recoDSADxyError_.push_back(mu_i->dxyError());
        nt.recoDSADz_.push_back(mu_i->dz(pv.position()));
        nt.recoDSADzError_.push_back(mu_i->dzError());
        nt.recoDSACharge_.push_back(mu_i->charge());
        nt.recoDSATrkChi2_.push_back(mu_i->normalizedChi2());
        nt.recoDSATrkNumPlanes_.push_back(mu_i->hitPattern().muonStationsWithValidHits());
        nt.recoDSATrkNumHits_.push_back(mu_i->hitPattern().numberOfValidMuonHits());
        nt.recoDSATrkNumDTHits_.push_back(mu_i->hitPattern().numberOfValidMuonDTHits());
        nt.recoDSATrkNumCSCHits_.push_back(mu_i->hitPattern().numberOfValidMuonCSCHits());
        
        // add muon timing info from custom-built muon object
        reco::MuonRef muon_i = muObjs1[i];

        reco::MuonTimeExtra time_info = (*dsaRecoMuTiming_.handle)[muon_i];
        nt.recoDSAInvBeta_.push_back(time_info.inverseBeta());
        nt.recoDSAInvBetaErr_.push_back(time_info.inverseBeta());
        nt.recoDSAFreeInvBeta_.push_back(time_info.freeInverseBeta());
        nt.recoDSAFreeInvBetaErr_.push_back(time_info.freeInverseBetaErr());
        nt.recoDSAtimeAtIpInOut_.push_back(time_info.timeAtIpInOut());
        nt.recoDSAtimeAtIpInOutErr_.push_back(time_info.timeAtIpInOutErr());
        nt.recoDSAtimeAtIpOutIn_.push_back(time_info.timeAtIpOutIn());
        nt.recoDSAtimeAtIpOutInErr_.push_back(time_info.timeAtIpOutInErr());
        nt.recoDSAtimingNdof_.push_back(time_info.nDof());
    }
    // Add leading 4 pt muons
    for (size_t i = 0; i < 4; i++) {
        if (i >= muTracks2.size()) break;
        reco::TrackRef mu_i = muTracks2[i];

        nt.recoGMPt_.push_back(mu_i->pt());
        nt.recoGMPtError_.push_back(mu_i->ptError());
        nt.recoGMEta_.push_back(mu_i->eta());
        nt.recoGMEtaError_.push_back(mu_i->etaError());
        nt.recoGMPhi_.push_back(mu_i->phi());
        nt.recoGMPhiError_.push_back(mu_i->phiError());
        nt.recoGMDxy_.push_back(mu_i->dxy(pv.position()));
        nt.recoGMDxyError_.push_back(mu_i->dxyError());
        nt.recoGMDz_.push_back(mu_i->dz(pv.position()));
        nt.recoGMDzError_.push_back(mu_i->dzError());
        nt.recoGMCharge_.push_back(mu_i->charge());
        nt.recoGMTrkChi2_.push_back(mu_i->normalizedChi2());
        nt.recoGMTrkNumPlanes_.push_back(mu_i->hitPattern().muonStationsWithValidHits());
        nt.recoGMTrkNumHits_.push_back(mu_i->hitPattern().numberOfValidMuonHits());
        nt.recoGMTrkNumDTHits_.push_back(mu_i->hitPattern().numberOfValidMuonDTHits());
        nt.recoGMTrkNumCSCHits_.push_back(mu_i->hitPattern().numberOfValidMuonCSCHits());

        reco::MuonRef muon_i = muObjs2[i];
        nt.recoGMIsPF_.push_back(muon_i->isPFMuon());
        nt.recoGMPFIso_.push_back((muon_i->pfIsolationR04().sumChargedHadronPt + std::max(0., 
                        muon_i->pfIsolationR04().sumNeutralHadronEt + muon_i->pfIsolationR04().sumPhotonEt
                        - 0.5*muon_i->pfIsolationR04().sumPUPt))/muon_i->pt());
        nt.recoGMTrkIso_.push_back(muon_i->isolationR03().sumPt/muon_i->pt());
        nt.recoGMTrkNumPixelHit_.push_back(muon_i->innerTrack()->hitPattern().numberOfValidPixelHits());
        nt.recoGMTrkNumTrkLayers_.push_back(muon_i->innerTrack()->hitPattern().trackerLayersWithMeasurement());

        recoil_px += mu_i->px();
        recoil_py += mu_i->py();
    }

    // Calculate recoil pt after looking at all collections
    nt.recoPFRecoilPt_ = sqrt(recoil_px*recoil_px + recoil_py*recoil_py);
    nt.recoPFRecoilPhi_ = atan2(recoil_py, recoil_px);

    // Apply Jet loose ID to jet collection, tag passes/fails on a side vector
    // Additionally mitigate HEM issue on chambers 15 and 16
    vector<bool> jetIDResults;
    for (auto const & recoJet : *recoJets_.handle) {
        bool jetIDResult = true;
        if (    recoJet.neutralHadronEnergyFraction() > 0.99 || 
                recoJet.neutralEmEnergyFraction() > 0.99 ||
                recoJet.numberOfDaughters() <= 1)
            jetIDResult = false;
        if (abs(recoJet.eta()) < 2.4)
            if (    recoJet.chargedHadronEnergyFraction() <= 0 ||
                    recoJet.chargedEmEnergyFraction() > 0.99 ||
                    recoJet.chargedMultiplicity() <= 0)
                jetIDResult = false;
        // as per monojet, also check extra cuts for leading jet
        auto index = &recoJet - &recoJets_.handle->at(0);
        if (index == 0)
            if (    recoJet.chargedHadronEnergyFraction() < 0.1 ||
                    recoJet.neutralHadronEnergyFraction() > 0.8)
                jetIDResult = false;
        jetIDResults.push_back(jetIDResult);
        // If passed jet is located in HEM region, veto the event
        // Has to happpen before jet ID, so don't check for jetIDResult
        double pt = recoJet.pt(), eta = recoJet.eta(), phi = recoJet.phi();
        if (pt > 30 && eta > -3.0 && eta < -1.4 && phi > -1.57 && phi < -0.87)
            nt.recoPFHEMFlag_ = true;
    }
    
    // Perform cross-cleaning in jet collection with good-quality GM muons
    // TODO: other option is that muon is coming out of jet, in that case
    // want a way to reject entire event since it's probably J/Psi
    // Nevermind: event selection ensures jets not collimated with MET+muons
    for (size_t i = 0; i < recoJets_.handle->size(); ++i) {
        // Skip failed ID jets
        if (!jetIDResults[i]) continue;
        reco::PFJetRef jet_i(recoJets_.handle, i);
        for (size_t j = 0; j < muGoodTracksIdx2.size(); ++j) {
            reco::TrackRef mu_j = muTracks2[muGoodTracksIdx2[j]];
            double dR = reco::deltaR(*jet_i, *mu_j);
            // if muon and jet match in dR fail the ID jet
            // because the jet is probably a muon instead
            if (dR < 0.4)
                jetIDResults[i] = false;
        }
    }
    // TODO: dSA muons dont contribute to PFJets so no need to cross-clean
    // However, possible muons from J/Psi in jets --> how to mitigate?
    // Nevermind: event selection ensures jets not collimated with MET+muons

    // Apply JEC+JER to jets that pass ID
    // Need to re-order jets by pT after this
    // vector key: index that will be changed after re-ordering
    // vector value1: corrected pT, value2: original key to refer back to
    
    JME::JetResolution resolution = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
    JME::JetResolutionScaleFactor resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, "AK4PFchs");

    vector<std::pair<std::unique_ptr<reco::PFJet>, size_t>> correctedJets;
    vector<std::pair<std::unique_ptr<reco::PFJet>, size_t>> correctedJetsJESUp;
    vector<std::pair<std::unique_ptr<reco::PFJet>, size_t>> correctedJetsJESDown;
    vector<std::pair<std::unique_ptr<reco::PFJet>, size_t>> correctedJetsJERUp;
    vector<std::pair<std::unique_ptr<reco::PFJet>, size_t>> correctedJetsJERDown;

    // Initialize MET corrections due to JES/JER
    float corr_METpx = PFMETr->px(), corr_METpy = PFMETr->py();
    float corr_METpx_JESUp = PFMETr->px(), corr_METpy_JESUp = PFMETr->py();
    float corr_METpx_JESDown = PFMETr->px(), corr_METpy_JESDown = PFMETr->py();
    float corr_METpx_JERUp = PFMETr->px(), corr_METpy_JERUp = PFMETr->py();
    float corr_METpx_JERDown = PFMETr->px(), corr_METpy_JERDown = PFMETr->py();

    for (size_t i = 0; i < recoJets_.handle->size(); ++i) {
        reco::PFJetRef jet_i(recoJets_.handle, i);

        // before JEC, check for EE noise and see if jet is in
        // critical region
        if (jet_i->pt() < 50 && abs(jet_i->eta()) > 2.65 && abs(jet_i->eta()) < 3.139) {
            nt.recoPFMETEEDeltaPx_ += jet_i->px();
            nt.recoPFMETEEDeltaPy_ += jet_i->py();
        }

        double jec = jetCorrector_.handle->correction(*jet_i);
        std::unique_ptr<reco::PFJet> corr_jet_i(jet_i->clone());
        std::unique_ptr<reco::PFJet> corr_jet_i_jes_up(jet_i->clone());
        std::unique_ptr<reco::PFJet> corr_jet_i_jes_down(jet_i->clone());
        std::unique_ptr<reco::PFJet> corr_jet_i_jer_up(jet_i->clone());
        std::unique_ptr<reco::PFJet> corr_jet_i_jer_down(jet_i->clone());
        // For MET corrections due to jet smearing
        std::unique_ptr<reco::PFJet> corr_jet_i_jes_only(jet_i->clone());

        corr_jet_i->scaleEnergy(jec);
        corr_jet_i_jes_only->scaleEnergy(jec);
        corr_jet_i_jes_up->scaleEnergy(jec);
        corr_jet_i_jes_down->scaleEnergy(jec);
        corr_jet_i_jer_up->scaleEnergy(jec);
        corr_jet_i_jer_down->scaleEnergy(jec);

        jecUnc_->setJetEta(corr_jet_i->eta());
        jecUnc_->setJetPt(corr_jet_i->pt());
        double uncUp = jecUnc_->getUncertainty(true); // true: Up direction
        jecUnc_->setJetEta(corr_jet_i->eta());
        jecUnc_->setJetPt(corr_jet_i->pt());
        double uncDown = jecUnc_->getUncertainty(false); // false: Down direction

        double jesUp = 1 + uncUp;
        double jesDown = 1 - uncDown;

        corr_jet_i_jes_up->scaleEnergy(jesUp);
        corr_jet_i_jes_down->scaleEnergy(jesDown);

        // Only do jet smearing (JES) if not working on data!
        double smearFactor = 1., smearFactor_up = 1., smearFactor_down = 1.;
        if (!isData) {
        
            double jet_resolution = resolution.getResolution({{JME::Binning::JetPt, corr_jet_i->pt()}, {JME::Binning::JetEta, corr_jet_i->eta()}, {JME::Binning::Rho, *rho_.handle}});
            double jer_sf = resolution_sf.getScaleFactor({{JME::Binning::JetPt, corr_jet_i->pt()}, {JME::Binning::JetEta, corr_jet_i->eta()}, {JME::Binning::Rho, *rho_.handle}}, Variation::NOMINAL);
            double jer_sf_up = resolution_sf.getScaleFactor({{JME::Binning::JetPt, corr_jet_i->pt()}, {JME::Binning::JetEta, corr_jet_i->eta()}, {JME::Binning::Rho, *rho_.handle}}, Variation::UP);
            double jer_sf_down = resolution_sf.getScaleFactor({{JME::Binning::JetPt, corr_jet_i->pt()}, {JME::Binning::JetEta, corr_jet_i->eta()}, {JME::Binning::Rho, *rho_.handle}}, Variation::DOWN);

            // Try matching with gen jet
            double min_dR = std::numeric_limits<double>::infinity();
            const reco::GenJet* matched_genJet = nullptr;
            for (const auto& genJet: *genJets_.handle) {
                double dR = deltaR(genJet, *corr_jet_i);
                if (dR > min_dR)
                    continue;
                if (dR < 0.2) {
                    double dPt = std::abs(genJet.pt() - corr_jet_i->pt());
                    matched_genJet = &genJet;
                }
            }
            if (matched_genJet) {
                double dPt = corr_jet_i->pt() - matched_genJet->pt();
                smearFactor = 1 + (jer_sf - 1.) * dPt / corr_jet_i->pt();
                smearFactor_up = 1 + (jer_sf_up - 1.) * dPt / corr_jet_i->pt();
                smearFactor_down = 1 + (jer_sf_down - 1.) * dPt / corr_jet_i->pt();
            }
            if (!matched_genJet && jer_sf > 1) { 
                double sigma = jet_resolution * std::sqrt(jer_sf * jer_sf - 1);
                std::normal_distribution<> d(0, sigma);
                smearFactor = 1. + d(m_random_generator);
            }
            if (!matched_genJet && jer_sf_up > 1) {
                double sigma_up = jet_resolution * std::sqrt(jer_sf_up * jer_sf_up - 1);
                std::normal_distribution<> d_up(0, sigma_up);
                smearFactor_up = 1. + d_up(m_random_generator);
            }
            if (!matched_genJet && jer_sf_down > 1) {
                double sigma_down = jet_resolution * std::sqrt(jer_sf_down * jer_sf_down - 1);
                std::normal_distribution<> d_down(0, sigma_down);
                smearFactor_down = 1. + d_down(m_random_generator);
            }
            //else
            //   std::cout << "ERROR! Impossible to smear this jet. jer_sf: " << jer_sf << std::endl;

            if (corr_jet_i->energy() * smearFactor < 0.01) {
                double newSmearFactor = 0.01 / corr_jet_i->energy();
                //std::cout << "The smearing factor (" << smearFactor << ") is either negative or too small. Fixing it to " << newSmearFactor << " to avoid change of direction." << std::endl;
                smearFactor = newSmearFactor;
            }
            if (corr_jet_i->energy() * smearFactor_up < 0.01) {
                double newSmearFactor = 0.01 / corr_jet_i->energy();
                //std::cout << "The smearing factor (" << smearFactor << ") is either negative or too small. Fixing it to " << newSmearFactor << " to avoid change of direction." << std::endl;
                smearFactor_up = newSmearFactor;
            }
            if (corr_jet_i->energy() * smearFactor_down < 0.01) {
                double newSmearFactor = 0.01 / corr_jet_i->energy();
                //std::cout << "The smearing factor (" << smearFactor << ") is either negative or too small. Fixing it to " << newSmearFactor << " to avoid change of direction." << std::endl;
                smearFactor_down = newSmearFactor;
            }
        } // end if checking for isData

        corr_jet_i->scaleEnergy(smearFactor);
        corr_jet_i_jes_up->scaleEnergy(smearFactor);
        corr_jet_i_jes_down->scaleEnergy(smearFactor);
        corr_jet_i_jer_up->scaleEnergy(smearFactor_up);
        corr_jet_i_jer_down->scaleEnergy(smearFactor_down);
        //corr_jet_i_jes_only->scaleEnergy(smearFactor);

        // Before finishing, compute the MET correction due to JER alone (MET type-I already accounts for JEC)
        //corr_METpx -= (corr_jet_i_jer_only->px() - jet_i->px());
        //corr_METpy -= (corr_jet_i_jer_only->py() - jet_i->py());
        corr_METpx -= (corr_jet_i->px() - corr_jet_i_jes_only->px());
        corr_METpy -= (corr_jet_i->py() - corr_jet_i_jes_only->py());
        
        // OBSOLETE --> temporarily "un-apply" smear factor to get JEC-corrected-only jet that was used in the original type-I calculation of MET
        //corr_jet_i->scaleEnergy(1.0/smearFactor);
        
        corr_METpx_JESUp -= (corr_jet_i_jes_up->px() - corr_jet_i_jes_only->px());
        corr_METpy_JESUp -= (corr_jet_i_jes_up->py() - corr_jet_i_jes_only->py());
        corr_METpx_JESDown -= (corr_jet_i_jes_down->px() - corr_jet_i_jes_only->px());
        corr_METpy_JESDown -= (corr_jet_i_jes_down->py() - corr_jet_i_jes_only->py());
        corr_METpx_JERUp -= (corr_jet_i_jer_up->px() - corr_jet_i_jes_only->px());
        corr_METpy_JERUp -= (corr_jet_i_jer_up->py() - corr_jet_i_jes_only->py());
        corr_METpx_JERDown -= (corr_jet_i_jer_down->px() - corr_jet_i_jes_only->px());
        corr_METpy_JERDown -= (corr_jet_i_jer_down->py() - corr_jet_i_jes_only->py());
        
        //corr_jet_i->scaleEnergy(smearFactor);

        correctedJets.push_back(std::make_pair(std::move(corr_jet_i), i));
        correctedJetsJESUp.push_back(std::make_pair(std::move(corr_jet_i_jes_up), i));
        correctedJetsJESDown.push_back(std::make_pair(std::move(corr_jet_i_jes_down), i));
        correctedJetsJERUp.push_back(std::make_pair(std::move(corr_jet_i_jer_up), i));
        correctedJetsJERDown.push_back(std::make_pair(std::move(corr_jet_i_jer_down), i));
    }

    // Before sorting the jet collections by pT, calculate MET corrections for each as well
    nt.recoPFMETSmearingOnlyPt_ = std::sqrt(corr_METpx*corr_METpx + corr_METpy*corr_METpy);
    nt.recoPFMETSmearingOnlyPhi_ = atan2(corr_METpy, corr_METpx);
    nt.recoPFMETCorrectedPt_ = nt.recoPFMETSmearingOnlyPt_;
    nt.recoPFMETCorrectedPhi_ = nt.recoPFMETSmearingOnlyPhi_;
    nt.recoPFMETJESUpPt_ = std::sqrt(corr_METpx_JESUp*corr_METpx_JESUp + corr_METpy_JESUp*corr_METpy_JESUp);
    nt.recoPFMETJESUpPhi_ = atan2(corr_METpy_JESUp, corr_METpx_JESUp);
    nt.recoPFMETJESDownPt_ = std::sqrt(corr_METpx_JESDown*corr_METpx_JESDown + corr_METpy_JESDown*corr_METpy_JESDown);
    nt.recoPFMETJESDownPhi_ = atan2(corr_METpy_JESDown, corr_METpx_JESDown);
    nt.recoPFMETJERUpPt_ = std::sqrt(corr_METpx_JERUp*corr_METpx_JERUp + corr_METpy_JERUp*corr_METpy_JERUp);
    nt.recoPFMETJERUpPhi_ = atan2(corr_METpy_JERUp, corr_METpx_JERUp);
    nt.recoPFMETJERDownPt_ = std::sqrt(corr_METpx_JERDown*corr_METpx_JERDown + corr_METpy_JERDown*corr_METpy_JERDown);
    nt.recoPFMETJERDownPhi_ = atan2(corr_METpy_JERDown, corr_METpx_JERDown);

    auto reverseSortJets = [](const auto &a, const auto &b) {
        return (a.first->pt() > b.first->pt());
    };
    sort(correctedJets.begin(), correctedJets.end(), reverseSortJets);
    sort(correctedJetsJESUp.begin(), correctedJetsJESUp.end(), reverseSortJets);
    sort(correctedJetsJESDown.begin(), correctedJetsJESDown.end(), reverseSortJets);
    sort(correctedJetsJERUp.begin(), correctedJetsJERUp.end(), reverseSortJets);
    sort(correctedJetsJERDown.begin(), correctedJetsJERDown.end(), reverseSortJets);

    // Get 10 top leading jets info, sorted by corrected pT
    // Only pick jets that have passed loose ID and cross-cleaning
    nt.recoPFNJet_ = recoJets_.handle->size(); 
    nt.recoPFNPassIDJet_ = 0;
    nt.recoPFNHighPtJet_ = 0;
    
    for (size_t i = 0; i < correctedJets.size(); i++) {
        size_t index = correctedJets[i].second;
        reco::PFJet & jet_i = *(correctedJets[i].first);

        // Exclude jets that didn't pass ID above (ID checked against uncorrected jets)
        if (!jetIDResults[index]) continue;
        nt.recoPFNPassIDJet_++;

        // Use JES+JER-corrected quantities
        if (jet_i.pt() > 30) 
            nt.recoPFNHighPtJet_++;

        if (nt.recoPFJetCorrectedPt_.size() < 10 && jet_i.pt() > 30) {
            nt.recoPFJetCorrectedPt_.push_back(jet_i.pt());
            nt.recoPFJetCorrectedEta_.push_back(jet_i.eta());
            nt.recoPFJetCorrectedPhi_.push_back(jet_i.phi());
            // TODO: figure out how BtagProbb(b) collections actually behave
            // FIXME this might be problematic with the jet corrections, keep in mind
            if (bTagProbb_.handle->size() > i && bTagProbbb_.handle->size() > i)
                nt.recoPFJetCorrectedBTag_.push_back((*bTagProbb_.handle)[index].second + (*bTagProbbb_.handle)[index].second);
            else 
                nt.recoPFJetCorrectedBTag_.push_back(-9999);
            nt.recoPFJetCorrectedCHEF_.push_back(jet_i.chargedHadronEnergyFraction());
            nt.recoPFJetCorrectedNHEF_.push_back(jet_i.neutralHadronEnergyFraction());
            nt.recoPFJetCorrectedCEEF_.push_back(jet_i.chargedEmEnergyFraction());
            nt.recoPFJetCorrectedNEEF_.push_back(jet_i.neutralEmEnergyFraction());
            nt.recoPFJetCorrectedNumDaughters_.push_back(jet_i.numberOfDaughters());
            nt.recoPFJetCorrectedChargedMultiplicity_.push_back(jet_i.chargedMultiplicity());
        }
        // Add pt, eta, phi information for syst. uncert. collections
        jet_i = *(correctedJetsJESUp[i].first);
        if (nt.recoPFJetCorrectedJESUpPt_.size() < 10 && jet_i.pt() > 30) {
            nt.recoPFJetCorrectedJESUpPt_.push_back(jet_i.pt());
            nt.recoPFJetCorrectedJESUpEta_.push_back(jet_i.eta());
            nt.recoPFJetCorrectedJESUpPhi_.push_back(jet_i.phi());
        }
        jet_i = *(correctedJetsJESDown[i].first);
        if (nt.recoPFJetCorrectedJESDownPt_.size() < 10 && jet_i.pt() > 30) {
            nt.recoPFJetCorrectedJESDownPt_.push_back(jet_i.pt());
            nt.recoPFJetCorrectedJESDownEta_.push_back(jet_i.eta());
            nt.recoPFJetCorrectedJESDownPhi_.push_back(jet_i.phi());
        }
        jet_i = *(correctedJetsJERUp[i].first);
        if (nt.recoPFJetCorrectedJERUpPt_.size() < 10 && jet_i.pt() > 30) {
            nt.recoPFJetCorrectedJERUpPt_.push_back(jet_i.pt());
            nt.recoPFJetCorrectedJERUpEta_.push_back(jet_i.eta());
            nt.recoPFJetCorrectedJERUpPhi_.push_back(jet_i.phi());
        }
        jet_i = *(correctedJetsJERDown[i].first);
        if (nt.recoPFJetCorrectedJERDownPt_.size() < 10 && jet_i.pt() > 30) {
            nt.recoPFJetCorrectedJERDownPt_.push_back(jet_i.pt());
            nt.recoPFJetCorrectedJERDownEta_.push_back(jet_i.eta());
            nt.recoPFJetCorrectedJERDownPhi_.push_back(jet_i.phi());
        }
        // Also add uncorrected jet info for completeness
        reco::PFJetRef jet_ii(recoJets_.handle, i);
        if (nt.recoPFJetPt_.size() < 10 && jet_ii->pt() > 30) {
            nt.recoPFJetPt_.push_back(jet_ii->pt());
            nt.recoPFJetEta_.push_back(jet_ii->eta());
            nt.recoPFJetPhi_.push_back(jet_ii->phi());
        }
    }

    // Pick pair of muons with smallest vertex chi square fit for all collection combos
    edm::ESHandle<TransientTrackBuilder> theB;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB);
    KalmanVertexFitter kvf(true);

    auto computeVertices = [&](auto coll_1, auto coll_2, std::string type) {
        for (size_t i = 0; i < 4; i++) {
            for (size_t j = 0; j < 4; j++) {
                reco::TrackRef muon_i, muon_j;
                if (i < coll_1.size())
                    muon_i = coll_1[i];
                if (j < coll_2.size())
                    muon_j = coll_2[j];

                TransientVertex tv;
                if (muon_i.isNonnull() && muon_j.isNonnull() && i != j) {
                    vector<reco::TransientTrack> transient_tracks{};
                    transient_tracks.push_back(theB->build(muon_i));
                    transient_tracks.push_back(theB->build(muon_j));
                    tv = kvf.vertex(transient_tracks);
                }

                float vxy = -9999;
                float sigma_vxy = -9999;
                float vtx_chi2 = 999999;
                float vz = -9999;
                float dr = -9999;

                if (tv.isValid()) {
                    reco::Vertex vertex = reco::Vertex(tv);
                    vxy = sqrt(vertex.x()*vertex.x() + vertex.y()*vertex.y());
                    sigma_vxy = (1/vxy)*sqrt(vertex.x()*vertex.x()*vertex.xError()*vertex.xError() +
                            vertex.y()*vertex.y()*vertex.yError()*vertex.yError());
                    //sigma_vxy = (1/vxy)*(vertex.x()*vertex.xError() + vertex.y()*vertex.yError());
                    vtx_chi2 = vertex.normalizedChi2();
                    vz = vertex.z();
                    dr = reco::deltaR(*muon_i, *muon_j);
                }

                if (type == "dsadsa") {
                    nt.dsadsa_recoVtxReducedChi2_.push_back(vtx_chi2);
                    nt.dsadsa_recoVtxVxy_.push_back(vxy);
                    nt.dsadsa_recoVtxSigmaVxy_.push_back(sigma_vxy);
                    nt.dsadsa_recoVtxVz_.push_back(vz);
                    nt.dsadsa_recoVtxDr_.push_back(dr);
                }
                else if (type == "gmgm") {
                    nt.gmgm_recoVtxReducedChi2_.push_back(vtx_chi2);
                    nt.gmgm_recoVtxVxy_.push_back(vxy);
                    nt.gmgm_recoVtxSigmaVxy_.push_back(sigma_vxy);
                    nt.gmgm_recoVtxVz_.push_back(vz);
                    nt.gmgm_recoVtxDr_.push_back(dr);
                }
                else if (type == "dsagm") {
                    nt.dsagm_recoVtxReducedChi2_.push_back(vtx_chi2);
                    nt.dsagm_recoVtxVxy_.push_back(vxy);
                    nt.dsagm_recoVtxSigmaVxy_.push_back(sigma_vxy);
                    nt.dsagm_recoVtxVz_.push_back(vz);
                    nt.dsagm_recoVtxDr_.push_back(dr);
                }
            }
        }
    };
    // dSA-dSA
    computeVertices(muTracks1, muTracks1, "dsadsa");
    // GM-GM
    computeVertices(muTracks2, muTracks2, "gmgm");
    // dSA-GM
    computeVertices(muTracks1, muTracks2, "dsagm");

    // Do muon matching in dR (0.1) between dSA and GM collections
    for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
            reco::TrackRef dsa_i, gm_j;

            if (i < muTracks1.size())
                dsa_i = muTracks1[i];
            if (j < muTracks2.size())
                gm_j = muTracks2[j];

            float dr = -9999;
            if (dsa_i.isNonnull() && gm_j.isNonnull())
                dr = reco::deltaR(dsa_i->outerEta(), dsa_i->outerPhi(),
                        gm_j->outerEta(), gm_j->outerPhi());

            nt.recoGMdSAdR_.push_back(dr);
            if (dr < -9998)
                nt.recoGMdSAmatch_.push_back(-1);
            else if (dr < 0.1)
                nt.recoGMdSAmatch_.push_back(1);
            else
                nt.recoGMdSAmatch_.push_back(0);
        }
    }

    nt.recoPFMETCorrectedPt_ = std::sqrt(corr_METpx*corr_METpx + corr_METpy*corr_METpy);
    nt.recoPFMETCorrectedPhi_ = atan2(corr_METpy, corr_METpx);
    nt.recoPFMETJESUpPt_ = std::sqrt(corr_METpx_JESUp*corr_METpx_JESUp + corr_METpy_JESUp*corr_METpy_JESUp);
    nt.recoPFMETJESUpPhi_ = atan2(corr_METpy_JESUp, corr_METpx_JESUp);
    nt.recoPFMETJESDownPt_ = std::sqrt(corr_METpx_JESDown*corr_METpx_JESDown + corr_METpy_JESDown*corr_METpy_JESDown);
    nt.recoPFMETJESDownPhi_ = atan2(corr_METpy_JESDown, corr_METpx_JESDown);
    nt.recoPFMETJERUpPt_ = std::sqrt(corr_METpx_JERUp*corr_METpx_JERUp + corr_METpy_JERUp*corr_METpy_JERUp);
    nt.recoPFMETJERUpPhi_ = atan2(corr_METpy_JERUp, corr_METpx_JERUp);
    nt.recoPFMETJERDownPt_ = std::sqrt(corr_METpx_JERDown*corr_METpx_JERDown + corr_METpy_JERDown*corr_METpy_JERDown);
    nt.recoPFMETJERDownPhi_ = atan2(corr_METpy_JERDown, corr_METpx_JERDown);

    //recoNMatchedGBMvDSA_ = -1;
    //int nDoubleMatched = 0;
    //if (fFoundValidVertex) {
    //    recoNMatchedGBMvDSA_ = 0;
    //    for (size_t i = 0; i < muGoodTracksIdx.size(); i++) {
    //        if ((int)i != dSA1Idx && (int)i != dSA2Idx) continue;
    //        bool alreadyMatched = false;
    //        auto & track1 = muTracks[muGoodTracksIdx[i]];
    //        const reco::HitPattern &p1 = track1->hitPattern();
    //        //std::cout << "p1 number of valid muon hits: " << p1.numberOfValidMuonHits() << std::endl;
    //        //std::cout << "p1 number of stations: " << p1.muonStationsWithValidHits() << std::endl;
    //        //std::cout << "p1 pattern: "; p1.print(reco::HitPattern::TRACK_HITS);
    //        std::vector<int> stations1;
    //        // loop over the hits of track1 
    //        for (int k = 0; k < p1.numberOfAllHits(reco::HitPattern::TRACK_HITS); k++) {
    //            uint32_t hit1 = p1.getHitPattern(reco::HitPattern::TRACK_HITS, k);
    //            if (!p1.validHitFilter(hit1)) continue;
    //            stations1.push_back(hit1);
    //        }
    //        stations1.erase(std::unique(stations1.begin(), stations1.end()), stations1.end());
    //        std::sort(stations1.begin(), stations1.end());
    //        for (size_t j = 0; j < muGoodTracksIdx2.size(); j++) {
    //            size_t nCommonChambers = 0;
    //            auto & track2 = muTracks2[muGoodTracksIdx2[j]];
    //            const reco::HitPattern &p2 = track2->hitPattern();
    //            //std::cout << "p2 number of valid muon hits: " << p2.numberOfValidMuonHits() << std::endl;
    //            //std::cout << "p2 number of stations: " << p2.muonStationsWithValidHits() << std::endl;
    //            //std::cout << "p2 pattern: "; p2.print(reco::HitPattern::TRACK_HITS);
    //            std::vector<int> stations2;
    //            // loop over the hits of track2
    //            for (int l = 0; l < p2.numberOfAllHits(reco::HitPattern::TRACK_HITS); l++) {
    //                uint32_t hit2 = p2.getHitPattern(reco::HitPattern::TRACK_HITS, l);
    //                if (!p2.validHitFilter(hit2)) continue;
    //                stations2.push_back(hit2);
    //            }
    //            stations2.erase(std::unique(stations2.begin(), stations2.end()), stations2.end());
    //            std::sort(stations2.begin(), stations2.end());
    //            std::vector<int> intersect;
    //            std::set_intersection(stations1.begin(),stations1.end(),stations2.begin(),stations2.end(), std::back_inserter(intersect));
    //            nCommonChambers = intersect.size();
    //            nCommonChambers_.push_back(min(stations1.size(),stations2.size()) - nCommonChambers);
    //            //std::cout << "nCommonChambers: " << nCommonChambers << std::endl;
    //            int distanceMetric = -1;
    //            // number of chambers not equal
    //            if (std::abs((int)stations1.size() - (int)stations2.size()) >= 1) {
    //                // if all chambers of smaller set match, not ok
    //                // (to make it ok, change distanceMetric back to 1 here)
    //                if (nCommonChambers == min(stations1.size(), stations2.size())) 
    //                    distanceMetric = 1;
    //                // one chamber of smaller set deosnt match, ok/not ok:
    //                else if (nCommonChambers == min(stations1.size(),stations2.size())-1)
    //                    distanceMetric = 1;
    //                // at least two chambers of smaller set doesnt match, not ok
    //                else 
    //                    distanceMetric = 2;
    //            }
    //            // number of chambers equal
    //            else {
    //                // at most one chamber different, ok
    //                if ((stations1.size() - nCommonChambers) <= 1)
    //                    distanceMetric = 1;
    //                // more than one chamber different, not ok
    //                else 
    //                    distanceMetric = 2;
    //            }
    //            if (distanceMetric == 1) {
    //            //if (nCommonStations >=  min(p1.muonStationsWithValidHits(), p2.muonStationsWithValidHits())) {
    //                recoNMatchedGBMvDSA_++;
    //                if (alreadyMatched)
    //                    nDoubleMatched++;
    //                else
    //                    alreadyMatched = true;
    //            }
    //        }
    //    }
    //}
    //std::cout << "[AF] Number of matched muons: " << recoNMatchedGBMvDSA_ << std::endl;
    //std::cout << "[AF] Number of double-matchings: " << nDoubleMatched << std::endl;
    

    /****** GEN AND TAU INFO *******/

    if (!isData) {

        nt.nGen_ = (int)genParticles_.handle->size();
        
        // Gen weight
        nt.genwgt_ = genEventInfo_.handle->weight();

        // Pile-up
        for (const auto & pileupInfo : *pileupInfos_.handle) {
            if (pileupInfo.getBunchCrossing() == 0) {
                nt.genpuobs_ = pileupInfo.getPU_NumInteractions();
                nt.genputrue_ = pileupInfo.getTrueNumInteractions();
                break;
            }
        }

        for (auto & genParticle : *genParticles_.handle) {
            if (!genParticle.isHardProcess()) continue;
            nt.genID_.push_back(genParticle.pdgId());
            //genHardProcess_.push_back(genParticle->isHardProcess());
            nt.genCharge_.push_back(genParticle.charge());
            nt.genPt_.push_back(genParticle.pt());
            nt.genEta_.push_back(genParticle.eta());
            nt.genPhi_.push_back(genParticle.phi());
            nt.genPz_.push_back(genParticle.pz());
            nt.genEn_.push_back(genParticle.energy());
            nt.genVxy_.push_back(genParticle.vertex().rho());
            nt.genVz_.push_back(genParticle.vz());
            nt.genMass_.push_back(genParticle.mass());
        }

        // all gen jets
        for (auto genJet : *genJets_.handle) {
            nt.genJetPt_.push_back(genJet.pt());
            nt.genJetEta_.push_back(genJet.eta());
            nt.genJetPhi_.push_back(genJet.phi());
        }

        // Lead gen MET
        if (genMET_.handle->size() > 0) {
            reco::GenMETRef metRef(genMET_.handle, 0);
            nt.genLeadMETPt_ = metRef->pt();
            nt.genLeadMETPhi_ = metRef->phi();
        }

        nt.FillGenTree();

        // Add taus

        // Select loosely reco taus (though for triggring purposes only)
        vector<reco::PFTauRef> selectedRecoTaus;
        for (size_t i = 0; i < recoPFTaus_.handle->size(); i++) {
            reco::PFTauRef tau(recoPFTaus_.handle, i);
            if (tau->pt() < 18)
                continue;
            selectedRecoTaus.push_back(tau);
        }

        // Select loosely l1 taus
        vector<l1t::TauRef> selectedL1Taus;
        for (size_t i = 0; i < L1Taus_.handle->size(); i++) {
            l1t::TauRef tau(L1Taus_.handle, i);
            if (tau->pt() < 15 || L1Taus_.handle->getFirstBX() != 0)
                continue;
            selectedL1Taus.push_back(tau);
        }

        // Select gen taus
        vector<reco::GenParticleRef> selectedGenTaus;
        for (size_t i = 0; i < genParticles_.handle->size(); i++) {
            reco::GenParticleRef genP(genParticles_.handle, i);
            if (abs(genP->pdgId()) != 15 || genP->status() != 2)
                continue;
            selectedGenTaus.push_back(genP);
        }

        // Compute visible p4 from gen taus
        vector<math::XYZTLorentzVector> genVisp4;
        for (auto & genTau : selectedGenTaus) {
            auto visp4 = genTau->p4();
            for (size_t i = 0; i < genTau->numberOfDaughters(); i++) {
                auto child = genTau->daughter(i);
                auto pdg = abs(child->pdgId());
                if (pdg == 12 || pdg == 14 || pdg == 16) // remove if neutrino
                    visp4 -= child->p4();
            }
            genVisp4.push_back(visp4);
        }

        // Match arbitrary collection to gen taus
        auto findBestMatch = [&](auto & coll, map<size_t, size_t> & matches) {
            vector<size_t> used_indices{};
            for (auto & genTau : selectedGenTaus) {
                size_t index_gen = &genTau - &selectedGenTaus.at(0);
                auto min_dR = 9999.f;
                size_t idx_best = 0;
                for (auto & obj : coll) {
                    size_t index_obj = &obj - &coll.at(0);
                    if (std::find(used_indices.begin(), used_indices.end(), index_obj) != used_indices.end())
                        continue;
                    auto dR = deltaR(obj->p4(), genVisp4.at(index_gen));
                    if (dR < min_dR) {
                        min_dR = dR;
                        idx_best = index_obj;
                    }
                }
                if (min_dR < 0.3) {
                    matches[index_gen] = idx_best;
                    used_indices.push_back(idx_best);
                }
            }
        };

        map<size_t, size_t> gen_reco_matches;
        findBestMatch(selectedRecoTaus, gen_reco_matches);

        map<size_t, size_t> gen_l1_matches;
        findBestMatch(selectedL1Taus, gen_l1_matches);

        std::function<bool (reco::GenParticleRef, reco::GenParticleRef)> isAncestor = [&](reco::GenParticleRef ancestor, reco::GenParticleRef particle) {
            if (ancestor == particle)
                return true;
            for (size_t i = 0; i < particle->numberOfMothers(); i++) {
                auto parent = particle->motherRef(i);
                if (isAncestor(ancestor, parent))
                    return true;
            }
            return false;
        };

        // parentID hardcoded for now -- TODO: pass as parameter
        // int tauParentID = 9900012; // assume HNL
        // int tauParentID = 9990012; // assume HNL && Dirac
        int tauParentID = 9000006; // assume S->TauTau

        for (auto & genTau : selectedGenTaus) {
            size_t idx_gen = &genTau - &selectedGenTaus.at(0);

            vector<reco::GenParticleRef> tauParents;
            for (size_t i = 0; i < genParticles_.handle->size(); i++) {
                reco::GenParticleRef genP(genParticles_.handle, i);
                if (isAncestor(genP, genTau) && abs(genP->pdgId()) == tauParentID)
                    tauParents.push_back(genP);
            }

            if (tauParents.size() == 0)
                continue;

            auto bestParent = tauParents.at(0);

            // Fill out some information
            // Each gen tau is one index of the vectors below
            // Each entry in the TTree is one event
            // This is to align genT, recoT, and tauT,
            // so can add one tree as a friend of the other
            // and use each other's branches

            nt.tau_gen_pt.push_back(genTau->pt());
            nt.tau_gen_eta.push_back(genTau->eta());
            nt.tau_gen_phi.push_back(genTau->phi());
            nt.tau_gen_charge.push_back(genTau->charge());
            // nt.tau_gen_decaymode.push_back(); // TODO add this one
            nt.tau_gen_vis_mass.push_back(genVisp4.at(idx_gen).mass());
            nt.tau_gen_vis_pt.push_back(genVisp4.at(idx_gen).pt());
            nt.tau_gen_vis_eta.push_back(genVisp4.at(idx_gen).eta());
            nt.tau_gen_vis_phi.push_back(genVisp4.at(idx_gen).phi());

            float diff_vx = genTau->vx() - bestParent->vx();
            float diff_vy = genTau->vy() - bestParent->vy();
            float diff_vz = genTau->vz() - bestParent->vz();
            
            nt.tau_gen_lxy.push_back(sqrt(pow(diff_vx, 2) + pow(diff_vy, 2)));
            nt.tau_gen_l3d.push_back(sqrt(pow(diff_vx, 2) + pow(diff_vy, 2) + pow(diff_vz, 2)));
            nt.tau_gen_vx.push_back(bestParent->vx());
            nt.tau_gen_vy.push_back(bestParent->vy());
            nt.tau_gen_vz.push_back(bestParent->vz());

            ROOT::Math::XYZVectorF P(genTau->px(), genTau->py(), 0);
            ROOT::Math::XYZVectorF L(diff_vx, diff_vy, 0);
            float cosXY = ROOT::Math::VectorUtil::CosTheta(P, L);
            nt.tau_gen_cosxy.push_back(cosXY);

            nt.tau_gen_parent_ct.push_back(TMath::C() * nt.tau_gen_l3d.back() * bestParent->mass() / bestParent->p());
            nt.tau_gen_parent_ct2d.push_back(TMath::C() * nt.tau_gen_lxy.back() * bestParent->mass() / bestParent->pt());
            nt.tau_gen_parent_mass.push_back(bestParent->mass());

            if (gen_reco_matches.find(idx_gen) != gen_reco_matches.end()) {
                auto reco_tau = selectedRecoTaus.at(gen_reco_matches[idx_gen]);
                nt.tau_reco_mass.push_back(reco_tau->mass());
                nt.tau_reco_pt.push_back(reco_tau->pt());
                nt.tau_reco_eta.push_back(reco_tau->eta());
                nt.tau_reco_phi.push_back(reco_tau->phi());
                nt.tau_reco_charge.push_back(reco_tau->charge());
                nt.tau_reco_vx.push_back(reco_tau->vx());
                nt.tau_reco_vy.push_back(reco_tau->vy());
                nt.tau_reco_vz.push_back(reco_tau->vz());
                // nt.tau_reco_decaymode.push_back(); // TODO add this one
            }
            else {
                nt.tau_reco_mass.push_back(-99.f);
                nt.tau_reco_pt.push_back(-99.f);
                nt.tau_reco_eta.push_back(-99.f);
                nt.tau_reco_phi.push_back(-99.f);
                nt.tau_reco_charge.push_back(-99);
                nt.tau_reco_vx.push_back(-1000.f);
                nt.tau_reco_vy.push_back(-1000.f);
                nt.tau_reco_vz.push_back(-1000.f);

            }

            if (gen_l1_matches.find(idx_gen) != gen_l1_matches.end()) {
                auto l1_tau = selectedL1Taus.at(gen_l1_matches[idx_gen]);
                nt.tau_l1_pt.push_back(l1_tau->pt());
                nt.tau_l1_eta.push_back(l1_tau->eta());
                nt.tau_l1_phi.push_back(l1_tau->phi());
                nt.tau_l1_charge.push_back(l1_tau->charge());
                nt.tau_l1_hwIso.push_back(l1_tau->hwIso());
            }
            else {
                nt.tau_l1_pt.push_back(-99.f);
                nt.tau_l1_eta.push_back(-99.f);
                nt.tau_l1_phi.push_back(-99.f);
                nt.tau_l1_charge.push_back(-99);
                nt.tau_l1_hwIso.push_back(-99.f);
            }

        }

        nt.FillTauTree();
    
    }

    // auto genTaus = 

    nt.FillRecoTree();

    return;
}

void ntupleAnalyzer::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}

void ntupleAnalyzer::endJob() {}

// define this as a plug-in
DEFINE_FWK_MODULE(ntupleAnalyzer);