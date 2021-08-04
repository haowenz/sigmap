#pragma once

#include <vector>

#include "nanoflann.hpp"

struct Point {
  uint64_t position;
  float value;
  Point() {}
  Point(const float value) : value(value) {}
  Point(const uint64_t position, const float value)
      : position(position), value(value) {}
  bool operator<(const Point &p) const {
    return std::tie(value, position) < std::tie(p.value, p.position);
  }
};
// ===== This example shows how to use nanoflann with these types of containers:
// =======
// typedef std::vector<std::vector<double> > my_vector_of_vectors_t;
// typedef std::vector<Eigen::VectorXd> my_vector_of_vectors_t;   // This
// requires #include <Eigen/Dense>
// =====================================================================================

/** A simple vector-of-vectors adaptor for nanoflann, without duplicating the
 * storage. The i'th vector represents a point in the state space.
 *
 *  \tparam DIM If set to >0, it specifies a compile-time fixed dimensionality
 * for the points in the data set, allowing more compiler optimizations. \tparam
 * num_t The type of the point coordinates (typically, double or float). \tparam
 * Distance The distance metric to use: nanoflann::metric_L1,
 * nanoflann::metric_L2, nanoflann::metric_L2_Simple, etc. \tparam IndexType The
 * type for indices in the KD-tree index (typically, size_t of int)
 */
template <typename num_t = float, int DIM = -1,
          class Distance = nanoflann::metric_L2, typename IndexType = size_t>
struct SigmapAdaptor {
  typedef SigmapAdaptor<num_t, DIM, Distance> self_t;
  typedef
      typename Distance::template traits<num_t, self_t>::distance_t metric_t;
  typedef nanoflann::KDTreeSingleIndexAdaptor<metric_t, self_t, DIM, IndexType>
      index_t;

  index_t *index;  //! The kd-tree index for the user to call its methods as
                   //! usual with any other FLANN index.

  size_t dims_;
  /// Constructor: takes a const ref to the vector of vectors object with the
  /// data points
  SigmapAdaptor(const size_t dims /* dimensionality */,
                const std::vector<Point> &mat, const int leaf_max_size = 10)
      : m_data(mat) {
    assert(mat.size() != 0);  // && mat[0].size() != 0);
    // const size_t dims = mat[0].size();
    if (DIM > 0 && static_cast<int>(dims) != DIM)
      throw std::runtime_error(
          "Data set dimensionality does not match the 'DIM' template argument");
    dims_ = dims;
    index =
        new index_t(static_cast<int>(dims), *this /* adaptor */,
                    nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max_size));
  }

  ~SigmapAdaptor() { delete index; }

  const std::vector<Point> &m_data;

  /** Query for the \a num_closest closest points to a given point (entered as
   * query_point[0:dim-1]). Note that this is a short-cut method for
   * index->findNeighbors(). The user can also call index->... methods as
   * desired. \note nChecks_IGNORED is ignored but kept for compatibility with
   * the original FLANN interface.
   */
  inline void query(const num_t *query_point, const size_t num_closest,
                    IndexType *out_indices, num_t *out_distances_sq,
                    const int nChecks_IGNORED = 10) const {
    nanoflann::KNNResultSet<num_t, IndexType> resultSet(num_closest);
    resultSet.init(out_indices, out_distances_sq);
    index->findNeighbors(resultSet, query_point, nanoflann::SearchParams());
  }

  /** @name Interface expected by KDTreeSingleIndexAdaptor
   * @{ */

  const self_t &derived() const { return *this; }
  self_t &derived() { return *this; }

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const {
    return m_data.size() - dims_ + 1;
  }

  // Returns the dim'th component of the idx'th point in the class:
  inline num_t kdtree_get_pt(const size_t idx, const size_t dim) const {
    // return m_data[idx][dim];
    return m_data[idx + dim].value;
  }

  // Optional bounding-box computation: return false to default to a standard
  // bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned in
  //   "bb" so it can be avoided to redo it again. Look at bb.size() to find out
  //   the expected dimensionality (e.g. 2 or 3 for point clouds)
  template <class BBOX>
  bool kdtree_get_bbox(BBOX & /*bb*/) const {
    return false;
  }

  /** @} */

};  // end of SigmapAdaptor
