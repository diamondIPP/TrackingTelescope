#ifndef DoAlignment_h
#define DoAlignment_h


#include <string>
#include <TString.h>
class PSIFileReader;
namespace tel { class ProgressBar; }
class TH2F;
class TGraph;

class Alignment {

public:
  Alignment(std::string, TString, short);
  ~Alignment();
  short const TelescopeID;
  unsigned const NPlanes;
  void PreAlign();
  void InitHistograms();
  void ResetHistograms();
  void SaveHistograms(unsigned, int ind=-1);
  void SaveGraphs(unsigned);
  int Align();
  void PrintAligment();

private:
  std::string InFileName;
  std::string OutFileName;
  TString const PlotsDir;
  TString const OutDir;
  float const AngleThreshold;
  float const TotResThreshold;
  std::vector<float> XAlign, YAlign, ZAlign, RAlign;
  PSIFileReader * FR;
  PSIFileReader * InitFileReader();
  unsigned MaxEventNumber;
  tel::ProgressBar * ProgressBar;
  float Now;
  unsigned short const MaximumSteps;
  /** Histograms
      hResidual:    x=dX / y=dY
      hResidualXdY: x=X  / y=dY
      hResidualYdX: x=Y  / y=dX  */
  std::vector<TH2F> hResidual;
  std::vector<TH2F> hResidualXdY;
  std::vector<TH2F> hResidualYdX;
  std::vector<TGraph> gResidualXdY;
  std::vector<TGraph> gResidualYdX;
};

#endif // DoAlignment_h
