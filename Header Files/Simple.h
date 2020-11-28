#ifndef __Simple_h__
#define __Simple_h__

#include "TObject.h"
#include "TString.h"

#include <cmath>
#include <algorithm>
#include <string>

using namespace std;

/** \class Miiqtool
*/
class Miiqtool {

 public:

  unsigned int status;                      ///< nParticle()+nAntiCluster()*10+nBetaH()*100+nTrTrack()*1000+nTrRecHit()*10000+nTrdCluster()*1000000+nTofClusterH()*100000000

  int          event;
  int          utime;

  float        trk_q_inn;                   ///< Inner Tracker Charge
  float        trk_q_lay[9];                ///< Tracker Layer charge (inverted sign for bad status)
  float        trk_rig;                     ///< Choutko fit rigidities (FS) [GV]
  float        trk_chisqn[2];               ///< Choutko fit normalized Chi2 (FS | X, Y)

  float        tof_beta;                    ///< TOF Beta

  float        rich_beta;                   ///< RICH beta best estimator

  int          utime_rti;

  float        lf;
  float        cf;



  Miiqtool(){}
  virtual ~Miiqtool(){}
  bool IsStandalone();
  bool IsAnalysis();
  int IsInsideRich();

  int nParticle()    { return status%10; }
  int nAntiCluster() { return int(status/10)%10; }
  int nBetaH()       { return int(status/100)%10; }
  int nTrTrack()     { return int(status/1000)%10; }
  int nTrRecHit()    { return int(status/10000)%100; }
  int nTrdCluster()  { return int(status/1000000)%100; }
  int nTofClusterH() { return int(status/100000000)%100; }

  ClassDef(Miiqtool,1);
};



/** \class MiiqRTI
*/
class MiiqRTI {

 public:

  float        utime_rti;
  float        lf;
  float        cf;



  MiiqRTI(){}
  virtual ~MiiqRTI(){}

  ClassDef(MiiqRTI,1);
};

#endif
