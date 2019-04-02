#ifndef DoAlignment_h
#define DoAlignment_h


#include <string>
#include <TString.h>
class PSIFileReader;
namespace tel { class ProgressBar; }
class TH2F;

class Alignment {

public:
  Alignment(std::string, TString, short);
  ~Alignment() {};
  short const TelescopeID;
  unsigned const NPlanes;
  void PreAlign();
  void InitHistograms();
  void ClearHistograms();
  void SaveResiduals(unsigned);

private:
  std::string InFileName;
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
  /** Histograms
      hResidual:    x=dX / y=dY
      hResidualXdY: x=X  / y=dY
      hResidualYdX: x=Y  / y=dX  */
  std::vector<TH2F> hResidual;
  std::vector<TH2F> hResidualXdY;
  std::vector<TH2F> hResidualYdX;
};

int DoAlignment (std::string, TString, int);

#endif // DoAlignment_h
