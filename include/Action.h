//
// Created by micha on 28.02.21.
//
#ifndef TRACKINGTELESCOPE_ACTION_H
#define TRACKINGTELESCOPE_ACTION_H

#include <string>
#include <TString.h>
class PSIFileReader;

class Action {

public:
  Action(std::string , const TString &);
  ~Action() = default;

protected:
  const std::string in_file_name_;
  const TString run_number_;

  PSIFileReader * FR;
  PSIFileReader * InitFileReader() const;
};


#endif //TRACKINGTELESCOPE_ACTION_H
