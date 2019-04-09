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
  void PrintAlignment();
  std::vector<unsigned short> GetOrderedPlanes();
  static void FormatHistogram(TH1 *, const std::string&, float, const std::string&, float);

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
  unsigned MaxEventNumber;
  tel::ProgressBar * ProgressBar;
  float Now;
  unsigned short const MaximumSteps;
  std::vector<unsigned short> OrderedPlanes;
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
