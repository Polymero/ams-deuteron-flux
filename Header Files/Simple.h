#ifndef __Simple_h__
#define __Simple_h__

#include "TObject.h"
#include "TString.h"

#include <cmath>
#include <algorithm>
#include <string>

using namespace std;

// Simplified class for the NtpCompact data class
// Created on       23-10-20
// Last edited on   30-11-20
class Miiqtool {

 public:

  unsigned int status;                      ///< nParticle()+nAntiCluster()*10+nBetaH()*100+nTrTrack()*1000+nTrRecHit()*10000+nTrdCluster()*1000000+nTofClusterH()*100000000

  short int    sublvl1;                     ///< Pattern of LVL1 sub-triggers (8 bits)
  short int    trigpatt;                    ///< Pattern of trigger system members (16 bits)

  int          event;                       ///< Event
  int          utime;                       ///< JMDC unix time [s]

  float        trk_q_inn;                   ///< Inner Tracker Charge
  float        trk_q_lay[9];                ///< Tracker Layer charge (inverted sign for bad status)
  float        trk_rig;                     ///< Choutko fit rigidities (FS) [GV]
  float        trk_chisqn[2];               ///< Choutko fit normalized Chi2 (FS | X, Y)

  float        tof_beta;                    ///< TOF Beta
  float        tof_q_lay[4];                ///< TOF Charge of each layer

  float        rich_beta;                   ///< RICH beta best estimator
  short int    rich_select;                 ///< Javier selection (see Tools::RichQC)

  float        lf;                          ///< Livetime [0,1]
  float        cf;                          ///< Max geomagnetic cutoff in the field of view (Stoermer|40 degrees|+) [GV]

  Miiqtool(){}
  virtual ~Miiqtool(){}

  ClassDef(Miiqtool,1);
};



// Simplified class for the RTIInfo data class
// Created on       23-10-20
// Last edited on   28-10-20
class MiiqRTI {

 public:

  int          utime_rti;                 ///< JMDC unix time [s]
  float        lf;                        ///< Livetime [0,1]
  float        cf;                        ///< Max geomagnetic cutoff in the field of view (Stoermer|40 degrees|+) [GV]

  MiiqRTI(){}
  virtual ~MiiqRTI(){}

  ClassDef(MiiqRTI,1);
};

#endif
