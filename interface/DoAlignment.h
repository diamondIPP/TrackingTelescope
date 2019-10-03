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
  Alignment(std::string, const TString&, short, bool=false);
  ~Alignment();
  short const TelescopeID;
  unsigned const NPlanes;
  bool const AlignOnlyInnerPlanes;
  void PreAlign();
  int Align();
  void InitHistograms();
  void ResetHistograms();
  void SaveHistograms(unsigned, int=-1);
  void SaveAllHistograms(int=-1);
  void PrintAlignment();
  void PrintResiduals(const std::vector<unsigned short>&);
  std::vector<unsigned short> GetOrderedPlanes();
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
  PSIFileReader * FR;
  PSIFileReader * InitFileReader();
  void EventLoop(const std::vector<unsigned short>&);
  unsigned MaxEventNumber;
  tel::ProgressBar * ProgressBar;
  float Now;
  unsigned short const MaximumSteps;
  std::vector<unsigned short> OrderedPlanes;
  std::vector<unsigned short> InnerPlanes;
  std::vector<unsigned short> PlanesToAlign;
  /** Means (offsets) and RMS of the residual distributions */
  std::vector<std::pair<float, float>> fdX;
  std::vector<std::pair<float, float>> fdY;
  std::vector<std::pair<float, float>> fdA;
  /** Histograms
      hResidual:    x=dX / y=dY
      hResidualXdY: x=X  / y=dY
      hResidualYdX: x=Y  / y=dX  */
  std::vector<TH2F> hResidual;
//  std::vector<TH2F> hResidualXdY;
  std::vector<TProfile> hResidualXdY;
//  std::vector<TH2F> hResidualYdX;
  std::vector<TProfile> hResidualYdX;
  static float GetMaxAverageAngle(const std::vector<std::pair<float, float>> &);
  static float GetMaximumMagRes(const std::vector<std::pair<float, float>>&, const std::vector<std::pair<float, float>>&);
};

#endif // DoAlignment_h
