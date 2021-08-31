//
// Created by micha on 28.02.21.
//

#include "Action.h"
#include <utility>
#include "GetNames.h"
#include "PSIRootFileReader.h"
#include "PSIBinaryFileReader.h"

using namespace std;

Action::Action(string file_name, const TString & run_number): in_file_name_(move(file_name)), run_number_(run_number), FR(nullptr) {

}

PSIFileReader * Action::InitFileReader() const {
  PSIFileReader * tmp;
  if (IsROOTFile(in_file_name_)){
    tmp = new PSIRootFileReader(in_file_name_, false, true);
  } else {
    tmp = new PSIBinaryFileReader(in_file_name_);
  }
  tmp->GetAlignment()->SetErrors(tel::Config::telescope_id_, true);
  return tmp;
}
