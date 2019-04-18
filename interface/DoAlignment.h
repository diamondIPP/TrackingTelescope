#ifndef DoAlignment_h
#define DoAlignment_h


#include <string>
#include <TString.h>
class PSIFileReader;
namespace tel { class ProgressBar; }
class TH2F;
class TGraph;
class TH1;

class Alignment {

public:
  Alignment(std::string, const TString&, short, bool=false);
  ~Alignment();
  short const TelescopeID;
  unsigned const NPlanes;
  bool const AlignOnlyInnerPlanes;
  void PreAlign();
  void InitHistograms();
  void ResetHistograms();
  void SaveHistograms(unsigned, int ind=-1);
  int Align();
  void PrintAlignment();
  std::vector<unsigned short> GetOrderedPlanes();
  static void FormatHistogram(TH1 *, const std::string&, float, const std::string&, float);
  std::pair<float, float> GetMeanResiduals(unsigned short i_plane);
  std::pair<float, float> GetRMS(unsigned short i_plane);

private:
  std::string InFileName;
  std::string OutFileName;
  TString const PlotsDir;
  TString const OutDir;
  std::string FileType;
  float const AngleThreshold;
  float const TotResThreshold;
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
  /** Histograms
      hResidual:    x=dX / y=dY
      hResidualXdY: x=X  / y=dY
      hResidualYdX: x=Y  / y=dX  */
  std::vector<TH2F> hResidual;
  std::vector<TH2F> hResidualXdY;
  std::vector<TH2F> hResidualYdX;
};

#endif // DoAlignment_h
