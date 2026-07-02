#ifndef DST_TPC_BRANCHES_HH
#define DST_TPC_BRANCHES_HH

#include <iostream>
#include <vector>

#include <TFile.h>
#include <TKey.h>
#include <TObject.h>
#include <TTree.h>
#include <Rtypes.h>

namespace dst_tpc
{
namespace detail
{
    template <class T>
    bool Set(TTree *tree, const char *name, T *address, bool required = true)
    {
        if (!tree || !tree->GetBranch(name)) {
            if (required)
                std::cerr << "DstTPCBranches: missing branch " << name << std::endl;
            return !required;
        }

        tree->SetBranchAddress(name, address);
        return true;
    }
}

struct BasicInfo {
    Int_t status = 0;
    UInt_t run_number = 0;
    UInt_t event_number = 0;
    std::vector<double> *trig_pat = nullptr;
    std::vector<std::vector<double>> *trig_flag = nullptr;
    Int_t beam_flag = 0;
    std::vector<double> *clkTpc = nullptr;

    bool SetBranchAddresses(TTree *tree, bool required = true)
    {
        bool ok = true;
        ok &= detail::Set(tree, "status", &status, required);
        ok &= detail::Set(tree, "run_number", &run_number, required);
        ok &= detail::Set(tree, "event_number", &event_number, required);
        ok &= detail::Set(tree, "trig_pat", &trig_pat, required);
        ok &= detail::Set(tree, "trig_flag", &trig_flag, required);
        ok &= detail::Set(tree, "beam_flag", &beam_flag, required);
        ok &= detail::Set(tree, "clkTpc", &clkTpc, required);
        return ok;
    }
};

struct RawHits {
    Int_t nhTpc = 0;
    std::vector<double> *raw_hitpos_x = nullptr;
    std::vector<double> *raw_hitpos_y = nullptr;
    std::vector<double> *raw_hitpos_z = nullptr;
    std::vector<double> *raw_de = nullptr;
    std::vector<int> *raw_padid = nullptr;
    std::vector<int> *raw_layer = nullptr;
    std::vector<int> *raw_row = nullptr;

    bool SetBranchAddresses(TTree *tree, bool required = true)
    {
        bool ok = true;
        ok &= detail::Set(tree, "nhTpc", &nhTpc, required);
        ok &= detail::Set(tree, "raw_hitpos_x", &raw_hitpos_x, required);
        ok &= detail::Set(tree, "raw_hitpos_y", &raw_hitpos_y, required);
        ok &= detail::Set(tree, "raw_hitpos_z", &raw_hitpos_z, required);
        ok &= detail::Set(tree, "raw_de", &raw_de, required);
        ok &= detail::Set(tree, "raw_padid", &raw_padid, required);
        ok &= detail::Set(tree, "raw_layer", &raw_layer, required);
        ok &= detail::Set(tree, "raw_row", &raw_row, required);
        return ok;
    }
};

struct Clusters {
    Int_t nclTpc = 0;
    std::vector<double> *cluster_x = nullptr;
    std::vector<double> *cluster_y = nullptr;
    std::vector<double> *cluster_z = nullptr;
    std::vector<double> *cluster_de = nullptr;
    std::vector<int> *cluster_size = nullptr;
    std::vector<int> *cluster_layer = nullptr;
    std::vector<int> *cluster_row_center = nullptr;
    std::vector<double> *cluster_mrow = nullptr;
    std::vector<double> *cluster_de_center = nullptr;
    std::vector<double> *cluster_x_center = nullptr;
    std::vector<double> *cluster_y_center = nullptr;
    std::vector<double> *cluster_z_center = nullptr;

    bool SetBranchAddresses(TTree *tree, bool required = true)
    {
        bool ok = true;
        ok &= detail::Set(tree, "nclTpc", &nclTpc, required);
        ok &= detail::Set(tree, "cluster_x", &cluster_x, required);
        ok &= detail::Set(tree, "cluster_y", &cluster_y, required);
        ok &= detail::Set(tree, "cluster_z", &cluster_z, required);
        ok &= detail::Set(tree, "cluster_de", &cluster_de, required);
        ok &= detail::Set(tree, "cluster_size", &cluster_size, required);
        ok &= detail::Set(tree, "cluster_layer", &cluster_layer, required);
        ok &= detail::Set(tree, "cluster_row_center", &cluster_row_center, required);
        ok &= detail::Set(tree, "cluster_mrow", &cluster_mrow, required);
        ok &= detail::Set(tree, "cluster_de_center", &cluster_de_center, required);
        ok &= detail::Set(tree, "cluster_x_center", &cluster_x_center, required);
        ok &= detail::Set(tree, "cluster_y_center", &cluster_y_center, required);
        ok &= detail::Set(tree, "cluster_z_center", &cluster_z_center, required);
        return ok;
    }
};

struct HelixTracks {
    Int_t ntTpc = 0;
    Int_t effective_ntTpc = 0;
    std::vector<int> *nhtrack = nullptr;
    std::vector<int> *is_beam = nullptr;
    std::vector<int> *is_accidental = nullptr;
    std::vector<double> *chisqr = nullptr;
    std::vector<double> *helix_cx = nullptr;
    std::vector<double> *helix_cy = nullptr;
    std::vector<double> *helix_z0 = nullptr;
    std::vector<double> *helix_r = nullptr;
    std::vector<double> *helix_dz = nullptr;
    std::vector<double> *helix_theta_min = nullptr;
    std::vector<double> *helix_theta_max = nullptr;
    std::vector<double> *mom0_x = nullptr;
    std::vector<double> *mom0_y = nullptr;
    std::vector<double> *mom0_z = nullptr;
    std::vector<double> *mom0 = nullptr;
    std::vector<double> *dE = nullptr;
    std::vector<double> *dEdx = nullptr;
    std::vector<int> *pid = nullptr;
    std::vector<double> *dz_factor = nullptr;
    std::vector<int> *charge = nullptr;
    std::vector<double> *path = nullptr;

    bool SetBranchAddresses(TTree *tree, bool required = true)
    {
        bool ok = true;
        ok &= detail::Set(tree, "ntTpc", &ntTpc, required);
        ok &= detail::Set(tree, "effective_ntTpc", &effective_ntTpc, required);
        ok &= detail::Set(tree, "nhtrack", &nhtrack, required);
        ok &= detail::Set(tree, "is_beam", &is_beam, required);
        ok &= detail::Set(tree, "is_accidental", &is_accidental, required);
        ok &= detail::Set(tree, "chisqr", &chisqr, required);
        ok &= detail::Set(tree, "helix_cx", &helix_cx, required);
        ok &= detail::Set(tree, "helix_cy", &helix_cy, required);
        ok &= detail::Set(tree, "helix_z0", &helix_z0, required);
        ok &= detail::Set(tree, "helix_r", &helix_r, required);
        ok &= detail::Set(tree, "helix_dz", &helix_dz, required);
        ok &= detail::Set(tree, "helix_theta_min", &helix_theta_min, required);
        ok &= detail::Set(tree, "helix_theta_max", &helix_theta_max, required);
        ok &= detail::Set(tree, "mom0_x", &mom0_x, required);
        ok &= detail::Set(tree, "mom0_y", &mom0_y, required);
        ok &= detail::Set(tree, "mom0_z", &mom0_z, required);
        ok &= detail::Set(tree, "mom0", &mom0, required);
        ok &= detail::Set(tree, "dE", &dE, required);
        ok &= detail::Set(tree, "dEdx", &dEdx, required);
        ok &= detail::Set(tree, "pid", &pid, required);
        ok &= detail::Set(tree, "dz_factor", &dz_factor, required);
        ok &= detail::Set(tree, "charge", &charge, required);
        ok &= detail::Set(tree, "path", &path, required);
        return ok;
    }
};

struct HelixPairs {
    std::vector<std::vector<double>> *combi_id = nullptr;
    std::vector<std::vector<double>> *closeDistTpc = nullptr;
    std::vector<std::vector<double>> *vtxTpc = nullptr;
    std::vector<std::vector<double>> *vtyTpc = nullptr;
    std::vector<std::vector<double>> *vtzTpc = nullptr;
    std::vector<std::vector<double>> *mom_vtx = nullptr;
    std::vector<std::vector<double>> *mom_vty = nullptr;
    std::vector<std::vector<double>> *mom_vtz = nullptr;

    bool SetBranchAddresses(TTree *tree, bool required = true)
    {
        bool ok = true;
        ok &= detail::Set(tree, "combi_id", &combi_id, required);
        ok &= detail::Set(tree, "closeDistTpc", &closeDistTpc, required);
        ok &= detail::Set(tree, "vtxTpc", &vtxTpc, required);
        ok &= detail::Set(tree, "vtyTpc", &vtyTpc, required);
        ok &= detail::Set(tree, "vtzTpc", &vtzTpc, required);
        ok &= detail::Set(tree, "mom_vtx", &mom_vtx, required);
        ok &= detail::Set(tree, "mom_vty", &mom_vty, required);
        ok &= detail::Set(tree, "mom_vtz", &mom_vtz, required);
        return ok;
    }
};

struct TrackHits {
    std::vector<std::vector<double>> *hitlayer = nullptr;
    std::vector<std::vector<double>> *hitpos_x = nullptr;
    std::vector<std::vector<double>> *hitpos_y = nullptr;
    std::vector<std::vector<double>> *hitpos_z = nullptr;
    std::vector<std::vector<double>> *calpos_x = nullptr;
    std::vector<std::vector<double>> *calpos_y = nullptr;
    std::vector<std::vector<double>> *calpos_z = nullptr;
    std::vector<std::vector<double>> *residual = nullptr;
    std::vector<std::vector<double>> *residual_x = nullptr;
    std::vector<std::vector<double>> *residual_y = nullptr;
    std::vector<std::vector<double>> *residual_z = nullptr;
    std::vector<std::vector<double>> *helix_t = nullptr;
    std::vector<std::vector<double>> *theta_diff = nullptr;
    std::vector<std::vector<double>> *pathhit = nullptr;
    std::vector<std::vector<double>> *pathhit_cor = nullptr;

    bool SetBranchAddresses(TTree *tree, bool required = true)
    {
        bool ok = true;
        ok &= detail::Set(tree, "hitlayer", &hitlayer, required);
        ok &= detail::Set(tree, "hitpos_x", &hitpos_x, required);
        ok &= detail::Set(tree, "hitpos_y", &hitpos_y, required);
        ok &= detail::Set(tree, "hitpos_z", &hitpos_z, required);
        ok &= detail::Set(tree, "calpos_x", &calpos_x, required);
        ok &= detail::Set(tree, "calpos_y", &calpos_y, required);
        ok &= detail::Set(tree, "calpos_z", &calpos_z, required);
        ok &= detail::Set(tree, "residual", &residual, required);
        ok &= detail::Set(tree, "residual_x", &residual_x, required);
        ok &= detail::Set(tree, "residual_y", &residual_y, required);
        ok &= detail::Set(tree, "residual_z", &residual_z, required);
        ok &= detail::Set(tree, "helix_t", &helix_t, required);
        ok &= detail::Set(tree, "theta_diff", &theta_diff, required);
        ok &= detail::Set(tree, "pathhit", &pathhit, required);
        ok &= detail::Set(tree, "pathhit_cor", &pathhit_cor, required);
        return ok;
    }
};

struct TrackClusters {
    std::vector<std::vector<double>> *track_cluster_de = nullptr;
    std::vector<std::vector<double>> *track_cluster_size = nullptr;
    std::vector<std::vector<double>> *track_cluster_mrow = nullptr;
    std::vector<std::vector<double>> *track_cluster_de_center = nullptr;
    std::vector<std::vector<double>> *track_cluster_x_center = nullptr;
    std::vector<std::vector<double>> *track_cluster_y_center = nullptr;
    std::vector<std::vector<double>> *track_cluster_z_center = nullptr;
    std::vector<std::vector<double>> *track_cluster_row_center = nullptr;

    bool SetBranchAddresses(TTree *tree, bool required = true)
    {
        bool ok = true;
        ok &= detail::Set(tree, "track_cluster_de", &track_cluster_de, required);
        ok &= detail::Set(tree, "track_cluster_size", &track_cluster_size, required);
        ok &= detail::Set(tree, "track_cluster_mrow", &track_cluster_mrow, required);
        ok &= detail::Set(tree, "track_cluster_de_center", &track_cluster_de_center, required);
        ok &= detail::Set(tree, "track_cluster_x_center", &track_cluster_x_center, required);
        ok &= detail::Set(tree, "track_cluster_y_center", &track_cluster_y_center, required);
        ok &= detail::Set(tree, "track_cluster_z_center", &track_cluster_z_center, required);
        ok &= detail::Set(tree, "track_cluster_row_center", &track_cluster_row_center, required);
        return ok;
    }
};

struct Lambda {
    std::vector<double> *lambda_mass = nullptr;
    std::vector<double> *lambda_close_dist = nullptr;
    std::vector<double> *lambda_vtx_x = nullptr;
    std::vector<double> *lambda_vtx_y = nullptr;
    std::vector<double> *lambda_vtx_z = nullptr;
    std::vector<double> *lambda_mom_x = nullptr;
    std::vector<double> *lambda_mom_y = nullptr;
    std::vector<double> *lambda_mom_z = nullptr;
    std::vector<double> *lambda_target_to_vtx_x = nullptr;
    std::vector<double> *lambda_target_to_vtx_y = nullptr;
    std::vector<double> *lambda_target_to_vtx_z = nullptr;
    std::vector<double> *lambda_target_to_vtx_dot_mom = nullptr;

    bool SetBranchAddresses(TTree *tree, bool required = true)
    {
        bool ok = true;
        ok &= detail::Set(tree, "lambda_mass", &lambda_mass, required);
        ok &= detail::Set(tree, "lambda_close_dist", &lambda_close_dist, required);
        ok &= detail::Set(tree, "lambda_vtx_x", &lambda_vtx_x, required);
        ok &= detail::Set(tree, "lambda_vtx_y", &lambda_vtx_y, required);
        ok &= detail::Set(tree, "lambda_vtx_z", &lambda_vtx_z, required);
        ok &= detail::Set(tree, "lambda_mom_x", &lambda_mom_x, required);
        ok &= detail::Set(tree, "lambda_mom_y", &lambda_mom_y, required);
        ok &= detail::Set(tree, "lambda_mom_z", &lambda_mom_z, required);
        ok &= detail::Set(tree, "lambda_target_to_vtx_x", &lambda_target_to_vtx_x, required);
        ok &= detail::Set(tree, "lambda_target_to_vtx_y", &lambda_target_to_vtx_y, required);
        ok &= detail::Set(tree, "lambda_target_to_vtx_z", &lambda_target_to_vtx_z, required);
        ok &= detail::Set(tree, "lambda_target_to_vtx_dot_mom", &lambda_target_to_vtx_dot_mom, required);
        return ok;
    }
};
}

struct DstTPCBranches {
    dst_tpc::BasicInfo basic;
    dst_tpc::RawHits raw_hits;
    dst_tpc::Clusters clusters;
    dst_tpc::HelixTracks tracks;
    dst_tpc::HelixPairs pairs;
    dst_tpc::TrackHits track_hits;
    dst_tpc::TrackClusters track_clusters;
    dst_tpc::Lambda lambda;

    bool SetBranchAddresses(TTree *tree)
    {
        if (!tree) {
            std::cerr << "DstTPCBranches: null TTree" << std::endl;
            return false;
        }

        bool ok = true;
        ok &= basic.SetBranchAddresses(tree);
        ok &= raw_hits.SetBranchAddresses(tree);
        ok &= clusters.SetBranchAddresses(tree);
        ok &= tracks.SetBranchAddresses(tree);
        ok &= pairs.SetBranchAddresses(tree);
        ok &= track_hits.SetBranchAddresses(tree);
        ok &= track_clusters.SetBranchAddresses(tree);
        ok &= lambda.SetBranchAddresses(tree);
        return ok;
    }
};

inline TTree *GetTree(TFile *file, const char *tree_name = "tpc")
{
    if (!file || file->IsZombie())
        return nullptr;

    if (tree_name) {
        TTree *tree = dynamic_cast<TTree*>(file->Get(tree_name));
        if (tree)
            return tree;

        std::cerr << "GetTree: tree " << tree_name << " not found; trying first TTree" << std::endl;
    }

    TIter next(file->GetListOfKeys());
    while (TKey *key = dynamic_cast<TKey*>(next())) {
        TObject *obj = key->ReadObj();
        TTree *tree = dynamic_cast<TTree*>(obj);
        if (tree)
            return tree;
        delete obj;
    }

    return nullptr;
}

#endif
