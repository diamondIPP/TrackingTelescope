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
  Alignment(std::string, const TString&, short, bool=false, unsigned short=30, float=0.0005, float=0.005, unsigned long=0, int8_t=-1);
  ~Alignment();
  short const TelescopeID;
  unsigned const NPlanes;
  bool AlignOnlyTelescope = false;
    uint16_t alignStep = 0;
    int8_t silDUTRoc = -1;
    bool trackOnlyTelescope = true;
    bool alignmentFinished = false;
//  void PreAlign();
  int Align();
  void InitHistograms();
  void ResetHistograms();
  void SaveHistograms(unsigned, int=-1);
  void SaveAllHistograms(int=-1);
  void PrintAlignment();
    void ClearVectors();
  void PrintResiduals(const std::vector<unsigned short>&);
  std::vector<unsigned short> GetOrderedPlanes();
  std::vector<unsigned short> GetTelescopePlanes();
  std::vector<unsigned short> GetDiamondPlanes();
  template <typename Q>
  static void FormatHistogram(Q *, const std::string&, float, const std::string&, float, float=0, float=0, float=0, float=0);

private:
  std::string InFileName;
  std::string OutFileName;
  TString const PlotsDir;
  TString const OutDir;
  std::string FileType;
  float const AngleThreshold;
  float const ResThreshold;
    unsigned long MaxEvents;
  PSIFileReader * FR;
  PSIFileReader * InitFileReader();
  void EventLoop(const std::vector<unsigned short>&);
  unsigned long MaxEventNumber;
  tel::ProgressBar * ProgressBar;
  float Now;
  unsigned short const MaximumSteps;
  std::vector<unsigned short> OrderedPlanes;
  std::vector<unsigned short> InnerPlanes;
  std::vector<unsigned short> PlanesToAlign;
  std::vector<unsigned short> PlanesUnderTest; // No tracking for this planes
  std::vector<unsigned short> TelescopePlanes; // No tracking for this planes
  std::vector<unsigned short> DiaPlanes; // No tracking for this planes
  /** Means (offsets) and RMS of the residual distributions */
  std::vector<std::pair<float, float>> fdX;
  std::vector<std::pair<float, float>> fdY;
  std::vector<std::pair<float, float>> fdA;
    float nSigma = 10.;
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
  static float GetMaxAngle(const std::vector<std::pair<float, float>> &, std::vector<unsigned short>);
  static float GetMaximumMagRes(const std::vector<std::pair<float, float>>&, const std::vector<std::pair<float, float>>&, std::vector<unsigned short>);
};

#endif // DoAlignment_h
