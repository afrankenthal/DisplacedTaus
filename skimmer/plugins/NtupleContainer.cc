#include "NtupleContainer.hh"

NtupleContainer::NtupleContainer() : isData_(true) {}

NtupleContainer::~NtupleContainer() {}

void NtupleContainer::CreateTreeBranches() {

    recoT->Branch("event_num", &eventNum_);
    recoT->Branch("lumi_sec", &lumiSec_);
    recoT->Branch("run_num", &runNum_);
    recoT->Branch("npv", &npv_);
    recoT->Branch("MET_filters_fail_bits", &METFiltersFailBits_);
    recoT->Branch("trig_fired", &fired_);
    recoT->Branch("reco_n_dsa", &recoNDSA_);
    recoT->Branch("reco_n_good_dsa", &recoNGoodDSA_);
    recoT->Branch("reco_dsa_pt",  &recoDSAPt_);
    recoT->Branch("reco_dsa_pt_err",  &recoDSAPtError_);
    recoT->Branch("reco_dsa_eta", &recoDSAEta_);
    recoT->Branch("reco_dsa_eta_err", &recoDSAEtaError_);
    recoT->Branch("reco_dsa_phi", &recoDSAPhi_);
    recoT->Branch("reco_dsa_phi_err", &recoDSAPhiError_);
    recoT->Branch("reco_dsa_dxy", &recoDSADxy_);
    recoT->Branch("reco_dsa_dxy_err", &recoDSADxyError_);
    recoT->Branch("reco_dsa_dz",  &recoDSADz_);
    recoT->Branch("reco_dsa_dz_err",  &recoDSADzError_);
    recoT->Branch("reco_dsa_charge", &recoDSACharge_);
    recoT->Branch("reco_dsa_trk_chi2", &recoDSATrkChi2_);
    recoT->Branch("reco_dsa_trk_n_planes", &recoDSATrkNumPlanes_);
    recoT->Branch("reco_dsa_trk_n_hits", &recoDSATrkNumHits_);
    recoT->Branch("reco_dsa_trk_n_DT_hits", &recoDSATrkNumDTHits_);
    recoT->Branch("reco_dsa_trk_n_CSC_hits", &recoDSATrkNumCSCHits_);
    recoT->Branch("reco_dsa_inv_beta", &recoDSAInvBeta_);
    recoT->Branch("reco_dsa_inv_beta_err", &recoDSAInvBetaErr_);
    recoT->Branch("reco_dsa_free_inv_beta", &recoDSAFreeInvBeta_);
    recoT->Branch("reco_dsa_free_inv_beta_err", &recoDSAFreeInvBetaErr_);
    recoT->Branch("reco_dsa_time_at_ip_in_out", &recoDSAtimeAtIpInOut_);
    recoT->Branch("reco_dsa_time_at_ip_in_out_err", &recoDSAtimeAtIpInOutErr_);
    recoT->Branch("reco_dsa_time_at_ip_out_in", &recoDSAtimeAtIpOutIn_);
    recoT->Branch("reco_dsa_time_at_ip_out_in_err", &recoDSAtimeAtIpOutInErr_);
    recoT->Branch("reco_dsa_time_ndof", &recoDSAtimingNdof_);
    recoT->Branch("reco_dsa_idx0", &recoDSAIdx0_);
    recoT->Branch("reco_dsa_idx1", &recoDSAIdx1_);
    recoT->Branch("reco_n_gm", &recoNGM_);
    recoT->Branch("reco_n_good_gm", &recoNGoodGM_);
    recoT->Branch("reco_gm_pt",  &recoGMPt_);
    recoT->Branch("reco_gm_pt_err",  &recoGMPtError_);
    recoT->Branch("reco_gm_eta", &recoGMEta_);
    recoT->Branch("reco_gm_eta_err", &recoGMEtaError_);
    recoT->Branch("reco_gm_phi", &recoGMPhi_);
    recoT->Branch("reco_gm_phi_err", &recoGMPhiError_);
    recoT->Branch("reco_gm_dxy", &recoGMDxy_);
    recoT->Branch("reco_gm_dxy_err", &recoGMDxyError_);
    recoT->Branch("reco_gm_dz",  &recoGMDz_);
    recoT->Branch("reco_gm_dz_err",  &recoGMDzError_);
    recoT->Branch("reco_gm_charge", &recoGMCharge_);
    recoT->Branch("reco_gm_trk_chi2", &recoGMTrkChi2_);
    recoT->Branch("reco_gm_trk_n_planes", &recoGMTrkNumPlanes_);
    recoT->Branch("reco_gm_trk_n_hits", &recoGMTrkNumHits_);
    recoT->Branch("reco_gm_trk_n_DT_hits", &recoGMTrkNumDTHits_);
    recoT->Branch("reco_gm_trk_n_CSC_hits", &recoGMTrkNumCSCHits_);
    recoT->Branch("reco_gm_isPF", &recoGMIsPF_);
    recoT->Branch("reco_gm_PFIso", &recoGMPFIso_);
    recoT->Branch("reco_gm_TrkIso", &recoGMTrkIso_);
    recoT->Branch("reco_gm_trk_n_pix_hits", &recoGMTrkNumPixelHit_);
    recoT->Branch("reco_gm_trk_n_trk_layers", &recoGMTrkNumTrkLayers_);
    recoT->Branch("reco_n_gbmdsa_match", &recoNMatchedGBMvDSA_);
    recoT->Branch("reco_gbmdsa_dR", &recoGMdSAdR_);
    recoT->Branch("reco_gbmdsa_match", &recoGMdSAmatch_);
    recoT->Branch("reco_sel_mu_pt", &selectedMuonsPt_);
    recoT->Branch("reco_sel_mu_pt_err", &selectedMuonsPtError_);
    recoT->Branch("reco_sel_mu_eta", &selectedMuonsEta_);
    recoT->Branch("reco_sel_mu_eta_err", &selectedMuonsEtaError_);
    recoT->Branch("reco_sel_mu_phi", &selectedMuonsPhi_);
    recoT->Branch("reco_sel_mu_phi_err", &selectedMuonsPhiError_);
    recoT->Branch("reco_sel_mu_dxy", &selectedMuonsDxy_);
    recoT->Branch("reco_sel_mu_dxy_error", &selectedMuonsDxyError_);
    recoT->Branch("reco_sel_mu_dz", &selectedMuonsDz_);
    recoT->Branch("reco_sel_mu_dz_error", &selectedMuonsDzError_);
    recoT->Branch("reco_sel_mu_charge", &selectedMuonsCharge_);
    recoT->Branch("reco_Mmumu",  &recoMmumu_);
    recoT->Branch("reco_METmu_dphi", &recoDeltaPhiMETMu_);
    recoT->Branch("reco_corr_METmu_dphi", &recoDeltaPhiCorrectedMETMu_);
    recoT->Branch("reco_pv_vx", &pvx_);
    recoT->Branch("reco_pv_vy", &pvy_);
    recoT->Branch("reco_pv_vz", &pvz_);
    recoT->Branch("reco_vtx_dsadsa_vxy", &dsadsa_recoVtxVxy_);
    recoT->Branch("reco_vtx_dsadsa_vz",  &dsadsa_recoVtxVz_);
    recoT->Branch("reco_vtx_dsadsa_sigmavxy", &dsadsa_recoVtxSigmaVxy_);
    recoT->Branch("reco_vtx_dsadsa_reduced_chi2", &dsadsa_recoVtxReducedChi2_);
    recoT->Branch("reco_vtx_dsadsa_dR",  &dsadsa_recoVtxDr_);
    recoT->Branch("reco_vtx_gmgm_vxy", &gmgm_recoVtxVxy_);
    recoT->Branch("reco_vtx_gmgm_vz",  &gmgm_recoVtxVz_);
    recoT->Branch("reco_vtx_gmgm_sigmavxy", &gmgm_recoVtxSigmaVxy_);
    recoT->Branch("reco_vtx_gmgm_reduced_chi2", &gmgm_recoVtxReducedChi2_);
    recoT->Branch("reco_vtx_gmgm_dR",  &gmgm_recoVtxDr_);
    recoT->Branch("reco_vtx_dsagm_vxy", &dsagm_recoVtxVxy_);
    recoT->Branch("reco_vtx_dsagm_vz",  &dsagm_recoVtxVz_);
    recoT->Branch("reco_vtx_dsagm_sigmavxy", &dsagm_recoVtxSigmaVxy_);
    recoT->Branch("reco_vtx_dsagm_reduced_chi2", &dsagm_recoVtxReducedChi2_);
    recoT->Branch("reco_vtx_dsagm_dR",  &dsagm_recoVtxDr_);
    recoT->Branch("reco_n_electron", &recoNElectron_);
    recoT->Branch("reco_n_good_electron", &recoNGoodElectron_);
    recoT->Branch("reco_electron_pt",  &recoElectronPt_);
    recoT->Branch("reco_electron_eta", &recoElectronEta_);
    recoT->Branch("reco_electron_phi", &recoElectronPhi_);
    recoT->Branch("reco_electron_vxy", &recoElectronVxy_);
    recoT->Branch("reco_electron_vz",  &recoElectronVz_);
    recoT->Branch("reco_electron_charge", &recoElectronCharge_);
    recoT->Branch("reco_electron_id_result", &recoElectronIDResult_);
    recoT->Branch("reco_n_photon", &recoNPhoton_);
    recoT->Branch("reco_n_good_photon", &recoNGoodPhoton_);
    recoT->Branch("reco_photon_pt",  &recoPhotonPt_);
    recoT->Branch("reco_photon_eta", &recoPhotonEta_);
    recoT->Branch("reco_photon_phi", &recoPhotonPhi_);
    recoT->Branch("reco_photon_id_result", &recoPhotonIDResult_);
    recoT->Branch("reco_PF_MET_pt", &recoPFMETPt_);
    recoT->Branch("reco_PF_MET_phi", &recoPFMETPhi_);
    recoT->Branch("reco_PF_MET_smearing_pt", &recoPFMETSmearingOnlyPt_);
    recoT->Branch("reco_PF_MET_smearing_phi", &recoPFMETSmearingOnlyPhi_);
    recoT->Branch("reco_PF_MET_corr_pt", &recoPFMETCorrectedPt_);
    recoT->Branch("reco_PF_MET_corr_phi", &recoPFMETCorrectedPhi_);
    recoT->Branch("reco_PF_MET_EE_delta_px", &recoPFMETEEDeltaPx_);
    recoT->Branch("reco_PF_MET_EE_delta_py", &recoPFMETEEDeltaPy_);
    recoT->Branch("reco_PF_MET_corr_JESUp_pt", &recoPFMETJESUpPt_);
    recoT->Branch("reco_PF_MET_corr_JESUp_phi", &recoPFMETJESUpPhi_);
    recoT->Branch("reco_PF_MET_corr_JESDown_pt", &recoPFMETJESDownPt_);
    recoT->Branch("reco_PF_MET_corr_JESDown_phi", &recoPFMETJESDownPhi_);
    recoT->Branch("reco_PF_MET_corr_JERUp_pt", &recoPFMETJERUpPt_);
    recoT->Branch("reco_PF_MET_corr_JERUp_phi", &recoPFMETJERUpPhi_);
    recoT->Branch("reco_PF_MET_corr_JERDown_pt", &recoPFMETJERDownPt_);
    recoT->Branch("reco_PF_MET_corr_JERDown_phi", &recoPFMETJERDownPhi_);
    recoT->Branch("reco_Calo_MET_pt", &recoCaloMETPt_);
    recoT->Branch("reco_Calo_MET_phi", &recoCaloMETPhi_);
    recoT->Branch("reco_PF_recoil_pt", &recoPFRecoilPt_);
    recoT->Branch("reco_PF_recoil_phi", &recoPFRecoilPhi_);
    recoT->Branch("reco_PF_n_jets", &recoPFNJet_);
    recoT->Branch("reco_PF_n_passID_jets", &recoPFNPassIDJet_);
    recoT->Branch("reco_PF_n_highPt_jets", &recoPFNHighPtJet_);
    recoT->Branch("reco_PF_jet_pt", &recoPFJetPt_);
    recoT->Branch("reco_PF_jet_eta", &recoPFJetEta_);
    recoT->Branch("reco_PF_jet_phi", &recoPFJetPhi_);
    recoT->Branch("reco_PF_jet_corr_pt", &recoPFJetCorrectedPt_);
    recoT->Branch("reco_PF_jet_corr_eta", &recoPFJetCorrectedEta_);
    recoT->Branch("reco_PF_jet_corr_phi", &recoPFJetCorrectedPhi_);
    recoT->Branch("reco_PF_jet_corr_BTag", &recoPFJetCorrectedBTag_);
    recoT->Branch("reco_PF_jet_corr_CHEF", &recoPFJetCorrectedCHEF_);
    recoT->Branch("reco_PF_jet_corr_NHEF", &recoPFJetCorrectedNHEF_);
    recoT->Branch("reco_PF_jet_corr_CEEF", &recoPFJetCorrectedCEEF_);
    recoT->Branch("reco_PF_jet_corr_NEEF", &recoPFJetCorrectedNEEF_);
    recoT->Branch("reco_PF_jet_corr_NumDaughters", &recoPFJetCorrectedNumDaughters_);
    recoT->Branch("reco_PF_jet_corr_ChargedMultipl", &recoPFJetCorrectedChargedMultiplicity_);
    recoT->Branch("reco_PF_jet_corr_JESUp_pt", &recoPFJetCorrectedJESUpPt_);
    recoT->Branch("reco_PF_jet_corr_JESUp_eta", &recoPFJetCorrectedJESUpEta_);
    recoT->Branch("reco_PF_jet_corr_JESUp_phi", &recoPFJetCorrectedJESUpPhi_);
    recoT->Branch("reco_PF_jet_corr_JESDown_pt", &recoPFJetCorrectedJESDownPt_);
    recoT->Branch("reco_PF_jet_corr_JESDown_eta", &recoPFJetCorrectedJESDownEta_);
    recoT->Branch("reco_PF_jet_corr_JESDown_phi", &recoPFJetCorrectedJESDownPhi_);
    recoT->Branch("reco_PF_jet_corr_JERUp_pt", &recoPFJetCorrectedJERUpPt_);
    recoT->Branch("reco_PF_jet_corr_JERUp_eta", &recoPFJetCorrectedJERUpEta_);
    recoT->Branch("reco_PF_jet_corr_JERUp_phi", &recoPFJetCorrectedJERUpPhi_);
    recoT->Branch("reco_PF_jet_corr_JERDown_pt", &recoPFJetCorrectedJERDownPt_);
    recoT->Branch("reco_PF_jet_corr_JERDown_eta", &recoPFJetCorrectedJERDownEta_);
    recoT->Branch("reco_PF_jet_corr_JERDown_phi", &recoPFJetCorrectedJERDownPhi_);
    recoT->Branch("reco_PF_HEM_flag", &recoPFHEMFlag_);
    recoT->Branch("reco_MHT_Pt", &MHTPt_);
    recoT->Branch("lheComments", "std::string", &lheComments);

    if (!isData_) {
        genT->Branch("event_num", &eventNum_);
        genT->Branch("gen_pu_obs", &genpuobs_);
        genT->Branch("gen_pu_true", &genputrue_);
        genT->Branch("gen_wgt", &genwgt_);
        genT->Branch("gen_ID", &genID_);
        //genT->Branch("gen_hard_process", &genHardProcess_);
        genT->Branch("gen_charge", &genCharge_);
        genT->Branch("gen_pt", &genPt_);
        genT->Branch("gen_eta", &genEta_);
        genT->Branch("gen_phi", &genPhi_);
        genT->Branch("gen_pz", &genPz_);
        genT->Branch("gen_energy", &genEn_);
        genT->Branch("gen_vxy", &genVxy_);
        genT->Branch("gen_vz", &genVz_);
        genT->Branch("gen_mass", &genMass_);
        genT->Branch("gen_jet_pt", &genJetPt_);
        genT->Branch("gen_jet_eta", &genJetEta_);
        genT->Branch("gen_jet_phi", &genJetPhi_);
        genT->Branch("gen_MET_pt", &genLeadMETPt_);
        genT->Branch("gen_MET_phi", &genLeadMETPhi_);
        genT->Branch("lheComments", "std::string", &lheComments);
    }

    tauT->Branch("tau_gen_pt", &tau_gen_pt);
    tauT->Branch("tau_gen_eta", &tau_gen_eta);
    tauT->Branch("tau_gen_phi", &tau_gen_phi);
    tauT->Branch("tau_gen_charge", &tau_gen_charge);
    tauT->Branch("tau_gen_vis_mass", &tau_gen_vis_mass);
    tauT->Branch("tau_gen_vis_pt", &tau_gen_vis_pt);
    tauT->Branch("tau_gen_vis_eta", &tau_gen_vis_eta);
    tauT->Branch("tau_gen_vis_phi", &tau_gen_vis_phi);
    tauT->Branch("tau_gen_lxy", &tau_gen_lxy);
    tauT->Branch("tau_gen_l3d", &tau_gen_l3d);
    tauT->Branch("tau_gen_cosxy", &tau_gen_cosxy);
    tauT->Branch("tau_gen_vx", &tau_gen_vx);
    tauT->Branch("tau_gen_vy", &tau_gen_vy);
    tauT->Branch("tau_gen_vz", &tau_gen_vz);
    tauT->Branch("tau_gen_parent_ct", &tau_gen_parent_ct);
    tauT->Branch("tau_gen_parent_ct2d", &tau_gen_parent_ct2d);
    tauT->Branch("tau_gen_parent_mass", &tau_gen_parent_mass);
    tauT->Branch("tau_reco_mass", &tau_reco_mass);
    tauT->Branch("tau_reco_pt", &tau_reco_pt);
    tauT->Branch("tau_reco_eta", &tau_reco_eta);
    tauT->Branch("tau_reco_phi", &tau_reco_phi);
    tauT->Branch("tau_reco_charge", &tau_reco_charge);
    tauT->Branch("tau_reco_vx", &tau_reco_vx);
    tauT->Branch("tau_reco_vy", &tau_reco_vy);
    tauT->Branch("tau_reco_vz", &tau_reco_vz);
    tauT->Branch("tau_l1_pt", &tau_l1_pt);
    tauT->Branch("tau_l1_eta", &tau_l1_eta);
    tauT->Branch("tau_l1_phi", &tau_l1_phi);
    tauT->Branch("tau_l1_charge", &tau_l1_charge);
    tauT->Branch("tau_l1_hwIso", &tau_l1_hwIso);

}

void NtupleContainer::ClearTreeBranches() {

    METFiltersFailBits_ = 0;

    recoDSAPt_.clear();
    recoDSAPtError_.clear();
    recoDSAEta_.clear();
    recoDSAEtaError_.clear();
    recoDSAPhi_.clear();
    recoDSAPhiError_.clear();
    recoDSADxy_.clear();
    recoDSADxyError_.clear();
    recoDSADz_.clear();
    recoDSADzError_.clear();
    recoDSACharge_.clear();
    recoDSATrkChi2_.clear();
    recoDSATrkNumPlanes_.clear();
    recoDSATrkNumHits_.clear();
    recoDSATrkNumDTHits_.clear();
    recoDSATrkNumCSCHits_.clear();
    recoDSAInvBeta_.clear();
    recoDSAInvBetaErr_.clear();
    recoDSAFreeInvBeta_.clear();
    recoDSAFreeInvBetaErr_.clear();
    recoDSAtimeAtIpInOut_.clear();
    recoDSAtimeAtIpInOutErr_.clear();
    recoDSAtimeAtIpOutIn_.clear();
    recoDSAtimeAtIpOutInErr_.clear();
    recoDSAtimingNdof_.clear();
    recoDSAIdx0_ = -9999;
    recoDSAIdx1_ = -9999;
    recoGMPt_.clear();
    recoGMPtError_.clear();
    recoGMEta_.clear();
    recoGMEtaError_.clear();
    recoGMPhi_.clear();
    recoGMPhiError_.clear();
    recoGMDxy_.clear();
    recoGMDxyError_.clear();
    recoGMDz_.clear();
    recoGMDzError_.clear();
    recoGMCharge_.clear();
    recoGMTrkChi2_.clear();
    recoGMTrkNumPlanes_.clear();
    recoGMTrkNumHits_.clear();
    recoGMTrkNumDTHits_.clear();
    recoGMTrkNumCSCHits_.clear();
    recoGMIsPF_.clear();
    recoGMPFIso_.clear();
    recoGMTrkIso_.clear();
    recoGMTrkNumPixelHit_.clear();
    recoGMTrkNumTrkLayers_.clear();
    dsadsa_recoVtxVxy_.clear();
    dsadsa_recoVtxVz_.clear();
    dsadsa_recoVtxSigmaVxy_.clear();
    dsadsa_recoVtxReducedChi2_.clear();
    dsadsa_recoVtxDr_.clear();
    gmgm_recoVtxVxy_.clear();
    gmgm_recoVtxVz_.clear();
    gmgm_recoVtxSigmaVxy_.clear();
    gmgm_recoVtxReducedChi2_.clear();
    gmgm_recoVtxDr_.clear();
    dsagm_recoVtxVxy_.clear();
    dsagm_recoVtxVz_.clear();
    dsagm_recoVtxSigmaVxy_.clear();
    dsagm_recoVtxReducedChi2_.clear();
    dsagm_recoVtxDr_.clear();
    recoGMdSAdR_.clear();
    recoGMdSAmatch_.clear();
    recoElectronPt_.clear();
    recoElectronEta_.clear();
    recoElectronPhi_.clear();
    recoElectronVxy_.clear();
    recoElectronVz_.clear();
    recoElectronCharge_.clear();
    recoElectronIDResult_.clear();
    recoPhotonPt_.clear();
    recoPhotonEta_.clear();
    recoPhotonPhi_.clear();
    recoPhotonIDResult_.clear();
    recoPFJetPt_.clear();
    recoPFJetEta_.clear();
    recoPFJetPhi_.clear();
    recoPFJetCorrectedPt_.clear();
    recoPFJetCorrectedEta_.clear();
    recoPFJetCorrectedPhi_.clear();
    recoPFJetCorrectedBTag_.clear();
    recoPFJetCorrectedCHEF_.clear();
    recoPFJetCorrectedNHEF_.clear();
    recoPFJetCorrectedCEEF_.clear();
    recoPFJetCorrectedNEEF_.clear();
    recoPFJetCorrectedNumDaughters_.clear();
    recoPFJetCorrectedChargedMultiplicity_.clear();
    recoPFJetCorrectedJESUpPt_.clear();
    recoPFJetCorrectedJESUpEta_.clear();
    recoPFJetCorrectedJESUpPhi_.clear();
    recoPFJetCorrectedJESDownPt_.clear();
    recoPFJetCorrectedJESDownEta_.clear();
    recoPFJetCorrectedJESDownPhi_.clear();
    recoPFJetCorrectedJERUpPt_.clear();
    recoPFJetCorrectedJERUpEta_.clear();
    recoPFJetCorrectedJERUpPhi_.clear();
    recoPFJetCorrectedJERDownPt_.clear();
    recoPFJetCorrectedJERDownEta_.clear();
    recoPFJetCorrectedJERDownPhi_.clear();
    recoPFHEMFlag_ = false;
    selectedMuonsPt_.clear();
    selectedMuonsPtError_.clear();
    selectedMuonsEta_.clear();
    selectedMuonsEtaError_.clear();
    selectedMuonsPhi_.clear();
    selectedMuonsPhiError_.clear();
    selectedMuonsDxy_.clear();
    selectedMuonsDxyError_.clear();
    selectedMuonsDz_.clear();
    selectedMuonsDzError_.clear();
    selectedMuonsCharge_.clear();

    recoNElectron_ = 0;
    recoNGoodElectron_ = 0;
    recoNPhoton_ = 0;
    recoNGoodPhoton_ = 0;

    fired_ = 0;
    recoPFMETPt_ = -9999;
    recoPFMETPhi_ = -9999;
    recoPFMETSmearingOnlyPt_ = -9999;
    recoPFMETSmearingOnlyPhi_ = -9999;
    recoPFMETCorrectedPt_ = -9999;
    recoPFMETCorrectedPhi_ = -9999;
    recoPFMETEEDeltaPx_ = 0.0;
    recoPFMETEEDeltaPy_ = 0.0;
    recoPFMETJESUpPt_ = -9999;
    recoPFMETJESUpPhi_ = -9999;
    recoPFMETJESDownPt_ = -9999;
    recoPFMETJESDownPhi_ = -9999;
    recoPFMETJERUpPt_ = -9999;
    recoPFMETJERUpPhi_ = -9999;
    recoPFMETJERDownPt_ = -9999;
    recoPFMETJERDownPhi_ = -9999;
    recoPFMETMuonEtFraction_ = -9999;
    recoCaloMETPt_ = -9999;
    recoCaloMETPhi_ = -9999;
    recoPFRecoilPt_ = -9999;
    recoPFRecoilPhi_ = -9999;
    recoMmumu_ = -9999;
    recoDeltaPhiMETMu_ = -9999;
    MHTPt_ = -9999;
    recoNMatchedGBMvDSA_ = -1;

    genID_.clear();
    //genHardProcess_.clear();
    genCharge_.clear();
    genPt_.clear();
    genEta_.clear();
    genPhi_.clear();
    genPz_.clear();
    genEn_.clear();
    genVxy_.clear();
    genVz_.clear();
    genMass_.clear();
    genJetPt_.clear();
    genJetEta_.clear();
    genJetPhi_.clear();

    // Pile-up and event genweight
    genpuobs_ = -9999;
    genputrue_ = -9999;
    genwgt_ = -9999;
    genLeadMETPt_ = -9999;
    genLeadMETPhi_ = -9999;

    // Tau info
    tau_gen_pt.clear();
    tau_gen_eta.clear();
    tau_gen_phi.clear();
    tau_gen_charge.clear();
    tau_gen_vis_mass.clear();
    tau_gen_vis_pt.clear();
    tau_gen_vis_eta.clear();
    tau_gen_vis_phi.clear();
    tau_gen_lxy.clear();
    tau_gen_l3d.clear();
    tau_gen_cosxy.clear();
    tau_gen_vx.clear();
    tau_gen_vy.clear();
    tau_gen_vz.clear();
    tau_gen_parent_ct.clear();
    tau_gen_parent_ct2d.clear();
    tau_gen_parent_mass.clear();
    tau_reco_mass.clear();
    tau_reco_pt.clear();
    tau_reco_eta.clear();
    tau_reco_phi.clear();
    tau_reco_charge.clear();
    tau_reco_vx.clear();
    tau_reco_vy.clear();
    tau_reco_vz.clear();
    tau_l1_pt.clear();
    tau_l1_eta.clear();
    tau_l1_phi.clear();
    tau_l1_charge.clear();
    tau_l1_hwIso.clear();

}