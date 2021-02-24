#ifndef DoAlignment_h
#define DoAlignment_h


#include <string>
#include <TString.h>
class PSIFileReader;
namespace tel { class ProgressBar; }
class TH2F;
class TGraph;
class TH1;
class TProfile;

class Alignment {

public:
  Alignment(const std::string&, const TString&, uint16_t, bool=false, uint16_t=30, float=0.0005, float=0.005, uint64_t=0, int16_t =-1);
  ~Alignment();

  /** variables */
  const uint16_t telescope_id_;
  const uint32_t n_planes_;
  bool align_only_telescope_;
  uint16_t align_step_ = 0;
  bool track_only_telescope_ = true;
  bool alignment_finished_ = true;

  /** methods */
  int Align();
  void InitHistograms();
  void ResetHistograms();
  void SaveHistograms(unsigned, int=-1);
  void SaveAllHistograms(int=-1);
  void PrintAlignment();
  void ClearVectors();
  void PrintResiduals(const std::vector<uint16_t>&);
  std::vector<uint16_t> GetOrderedPlanes();
  std::vector<uint16_t> GetTelescopePlanes();
  std::vector<uint16_t> GetDiamondPlanes();
  template <typename Q>
  static void FormatHistogram(Q *, const std::string&, float, const std::string&, float, float=0, float=0, float=0, float=0);

private:
  std::string OutFileName;
  TString const PlotsDir;
  TString const OutDir;
  std::string FileType;
  float const AngleThreshold;
  float const ResThreshold;
  uint64_t MaxEvents;
  PSIFileReader * FR;
  PSIFileReader * InitFileReader(const std::string&);
  void EventLoop(const std::vector<uint16_t>&);
  uint64_t MaxEventNumber;
  tel::ProgressBar * ProgressBar;
  float Now;
  uint16_t const MaximumSteps;
  int16_t silDUTRoc = -1;
  std::vector<uint16_t> OrderedPlanes;
  std::vector<uint16_t> InnerPlanes;
  std::vector<uint16_t> PlanesToAlign;
  std::vector<uint16_t> PlanesUnderTest; // No tracking for this planes
  std::vector<uint16_t> TelescopePlanes;
  std::vector<uint16_t> DiaPlanes;
  /** Means (offsets) and RMS of the residual distributions */
  std::vector<std::pair<float, float>> fdX;
  std::vector<std::pair<float, float>> fdY;
  std::vector<std::pair<float, float>> fdA;
  float n_sigma_;
  std::vector<float_t> maxResiduals;
  std::vector<float_t> maxAngles;

  void SetNextAlignmentStep();
  /** Histograms
      hResidual:    x=dX / y=dY
      hResidualXdY: x=X  / y=dY
      hResidualYdX: x=Y  / y=dX  */
  std::vector<TH2F> hResidual;
//  std::vector<TH2F> hResidualXdY;
  std::vector<TProfile> hResidualXdY;
//  std::vector<TH2F> hResidualYdX;
  std::vector<TProfile> hResidualYdX;
  std::vector<TGraph*> gMeanRes;
  std::vector<TGraph*> gAngleRes;
  static float GetMaxAngle(const std::vector<std::pair<float, float>> &, std::vector<uint16_t>);
  static float GetMaximumMagRes(const std::vector<std::pair<float, float>>&, const std::vector<std::pair<float, float>>&, std::vector<uint16_t>);
};

#endif // DoAlignment_h
