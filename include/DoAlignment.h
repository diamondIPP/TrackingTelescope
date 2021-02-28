#ifndef DoAlignment_h
#define DoAlignment_h

#include "Action.h"
namespace tel { class ProgressBar; }
class TH2F;
class TGraph;
class TH1;
class TProfile;

class Alignment: public Action {

public:
  Alignment(const std::string& in_file_name, const TString &run_number, uint16_t telescope_id, bool only_tel, uint16_t max_steps, float max_res, float max_angle, uint32_t max_events, int16_t sil_dut_roc);
  ~Alignment();

  /** variables */
  const uint16_t telescope_id_;
  const uint16_t n_planes_;
  const uint16_t n_telescope_planes = 4;
  bool align_only_telescope_;
  uint16_t at_step_;
  bool alignment_finished_;

  /** methods */
  int Align();
  void InitHistograms();
  void ResetHistograms();
  void SaveHistograms(unsigned, int=-1);
  void SaveAllHistograms(int=-1);
  void PrintAlignment();
  void PrintResiduals(const std::vector<uint16_t>&);
  std::vector<uint16_t> GetOrderedPlanes();
  std::vector<uint16_t> GetTelescopePlanes();
  std::vector<uint16_t> GetDiamondPlanes();
  template <typename Q>
  static void FormatHistogram(Q *, const std::string&, float, const std::string&, float, float=0, float=0, float=0, float=0);

private:
  TString const plots_dir_;
  std::string file_type_;
  float const angle_thresh_;
  float const res_thresh_;
  float const delta_fac_ = .1;
  void EventLoop(const std::vector<uint16_t>&);
  uint64_t max_event_number_;
  tel::ProgressBar * ProgressBar;
  float now_;
  uint16_t const maximum_steps_;
  int16_t sil_dut_roc_;
  std::vector<uint16_t> ordered_planes_;
  std::vector<uint16_t> inner_planes_;
  std::vector<uint16_t> planes_to_align_;
  std::vector<uint16_t> planes_under_test_; // No tracking for this planes
  std::vector<uint16_t> telescope_planes_;
  std::vector<uint16_t> dia_planes_;
  /** Means (offsets) and RMS of the residual distributions */
  std::vector<std::pair<float, float>> fdX;
  std::vector<std::pair<float, float>> fdY;
  std::vector<std::pair<float, float>> fdA;
  float max_sigma_, min_sigma_, n_sigma_;
  std::pair<float, float> last_max_res_;
  std::pair<float, float> delta_max_res_;

  void SetPlanes();
  void SetNextAlignmentStep();
  float ReduceSigma(uint16_t step) const;
  float CalcRes(uint16_t);

  /** Histograms -------------------
   *  hResidual:    x=dX / y=dY
   *  hResidualXdY: x=X  / y=dY
   *  hResidualYdX: x=Y  / y=dX  */
  std::vector<TH2F> hResidual;
  std::vector<TProfile> hResidualXdY;
  std::vector<TProfile> hResidualYdX;
  std::vector<TGraph*> g_res_mean_;
  std::vector<TGraph*> g_res_angle_;
  void CalcMaxResiduals();
  float GetMaxRes();
  float GetMaxAngle();
  static void ConfigROOT();
  void SaveGraphs();
  void InitGraphs();
};

#endif // DoAlignment_h
