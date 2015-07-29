#pragma once
#include "ps.h"
#include "base/progress.h"
namespace dmlc {
namespace linear {

using FeaID = ps::Key;
using Real = float;

class Progress : public VectorProgress {
 public:
  Progress() : VectorProgress(3, 3) {}
  virtual ~Progress() { }

  /// head string for printing
  virtual std::string HeadStr() {
    return "#example delta #ex    |w|_0       logloss     AUC    accuracy";
  }

  /// string for printing
  virtual std::string PrintStr(const IProgress* const prev) {
    Progress* const p = (Progress* const) prev;
    if (num_ex() == 0) return "";
    char buf[256];
    snprintf(buf, 256, "%7.2g  %7.2g  %11.6g  %8.6lf  %8.6lf  %8.6lf",
             (double)(p->num_ex() + num_ex()),
             (double)num_ex(),
             (double)(p->nnz_w() + nnz_w()),
             objv() / num_ex(),
             auc() / count(),
             acc() / count());
    return std::string(buf);
  }

  // mutator
  double& objv() { return fvec_[0]; }
  double& acc() { return fvec_[1]; }
  double& auc() { return fvec_[2]; }
  // double& weight2() { return fvec_[4]; }
  // double& wdelta2() { return fvec_[5]; }

  int64_t& count() { return ivec_[0]; }
  int64_t& num_ex() { return ivec_[1]; }
  int64_t& nnz_w() { return ivec_[2]; }

  int64_t num_ex() const { return ivec_[1]; }
  int64_t nnz_w() const { return ivec_[2]; }
};

} // namespace linear
} // namespace dmlc
