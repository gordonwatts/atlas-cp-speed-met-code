"""MET egamma preparation ABI shared by exporter and transport tools."""

REAL_INPUTS = """
pt eta phi mass cluster_eta cluster_phi eta_calo phi_calo pv_z
emb1_energy emb2_energy emb3_energy eme1_energy eme2_energy eme3_energy
emb1_eta emb2_eta eme1_eta eme2_eta
iso20 iso30 iso40 oldleak20 oldleak30 oldleak40 eadded_l2 eadded_l3 track_iso
track_d0 track_phi track_z0 track_vz track_var_d0 beam_sigma_x beam_sigma_y beam_sigma_xy
conv_x conv_y conv_pt1 conv_pt2 conv_trk1_q_over_p conv_trk1_theta conv_trk1_phi
conv_trk2_q_over_p conv_trk2_theta conv_trk2_phi leak_et rz_x rz_y
""".split()

INTEGER_INPUTS = """
has_cluster has_calo_moments sampling_pattern has_pv run_number random_run_number oq
author id_flag clean_flag iso_available_mask oldleak_available_mask timing_available_mask
has_track track_cov_diag_size has_conv_trk1 has_conv_trk2
conv_trk1_pixel conv_trk1_sct conv_trk2_pixel conv_trk2_sct has_conv_pts leak_key
""".split()

OUTPUTS = """
pt eta phi mass eta_bits oq_bits clean_bits pt_bits track_bits id isolated pteta_bits selected
iso20 iso30 iso40 leak20 leak30 leak40 correction_mask origin_applied
""".split()
