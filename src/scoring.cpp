#include <RcppEigen.h>
#include <progress.hpp>
#include <numeric>
#include <algorithm>


using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

template <typename S>
std::vector<size_t> sort_indexes(const std::vector<S> &v) {
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(),
                   [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  return idx;
}

// [[Rcpp::export]]
std::vector<double> ScoreHelper(
    Eigen::SparseMatrix<double> snn,
    Eigen::MatrixXd query_pca,
    Eigen::MatrixXd query_dists,
    Eigen::MatrixXd corrected_nns,
    int k_snn,
    bool subtract_first_nn,
    bool display_progress
) {
  std::vector<double> scores;
  // Loop over all query cells
  Progress p(snn.outerSize(), display_progress);
  for (int i=0; i < snn.outerSize(); ++i){
    p.increment();
    // create vectors to store the nonzero snn elements and their indices
    std::vector<double> nonzero;
    std::vector<size_t> nonzero_idx;
    for (Eigen::SparseMatrix<double>::InnerIterator it(snn, i); it; ++it) {
      nonzero.push_back(it.value());
      nonzero_idx.push_back(it.row());
    }
    // find the k_snn cells with the smallest non-zero edge weights to use in 
    // computing the transition probability bandwidth
    std::vector<size_t> nonzero_order = sort_indexes(nonzero);
    std::vector<double> bw_dists;
    int k_snn_i = k_snn;
    if (k_snn_i > nonzero_order.size()) k_snn_i = nonzero_order.size();
    for (int j = 0; j < nonzero_order.size(); ++j) {
      // compute euclidean distances to cells with small edge weights
      size_t cell = nonzero_idx[nonzero_order[j]];
      if(bw_dists.size() < k_snn_i || nonzero[nonzero_order[j]] == nonzero[nonzero_order[k_snn_i-1]]) {
        double res = (query_pca.col(cell) - query_pca.col(i)).norm();
        bw_dists.push_back(res);
      } else {
        break;
      }
    }
    // compute bandwidth as the mean distance of the farthest k_snn cells
    double bw;
    if (bw_dists.size() > k_snn_i) {
      std::sort(bw_dists.rbegin(), bw_dists.rend());
      bw = std::accumulate(bw_dists.begin(), bw_dists.begin() + k_snn_i, 0.0) / k_snn_i;
    } else {
      bw = std::accumulate(bw_dists.begin(), bw_dists.end(), 0.0) / bw_dists.size();
    }
    // compute transition probabilites 
    double first_neighbor_dist;
    // subtract off distance to first neighbor?
    if (subtract_first_nn) {
      first_neighbor_dist = query_dists(i, 1);
    } else {
      first_neighbor_dist = 0;
    }
    bw = bw - first_neighbor_dist;
    double q_tps = 0;
    for(int j = 0; j < query_dists.cols(); ++j) {
      q_tps += std::exp(-1 * (query_dists(i, j) - first_neighbor_dist) / bw);
    }
    q_tps = q_tps/(query_dists.cols());
    double c_tps = 0;
    for(int j = 0; j < corrected_nns.cols(); ++j) {
      c_tps += exp(-1 * ((query_pca.col(i) - query_pca.col(corrected_nns(i, j)-1)).norm() - first_neighbor_dist) / bw);
    }
    c_tps = c_tps/(corrected_nns.cols());
    scores.push_back(c_tps/q_tps);
  }
  return(scores);
}