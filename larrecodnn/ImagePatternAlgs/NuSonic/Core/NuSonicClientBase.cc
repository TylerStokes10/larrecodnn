#include "larrecodnn/ImagePatternAlgs/NuSonic/Core/NuSonicClientBase.h"
#include "canvas/Utilities/Exception.h"

NuSonicClientBase::NuSonicClientBase(const fhicl::ParameterSet& params,
                                     const std::string& clientName)
    : allowedTries_(params.get<unsigned>("allowedTries", 0)),
      clientName_(clientName) {
}

void NuSonicClientBase::start() {
  tries_ = 0;
}

void NuSonicClientBase::finish(bool success) {
  //retries are only allowed if no exception was raised
  if (!success) {
    ++tries_;
    //if max retries has not been exceeded, call evaluate again
    if (tries_ < allowedTries_) {
      evaluate();
      //avoid calling doneWaiting() twice
      return;
    }
    //prepare an exception if exceeded
    else {
      throw cet::exception("NuSonicClient")
           << "call failed after max " << tries_ << " tries" << std::endl;
    }
  }
}
