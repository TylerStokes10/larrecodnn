#ifndef NuSonic_Core_NuSonicClientBase
#define NuSonic_Core_NuSonicClientBase

#include "fhiclcpp/ParameterSet.h"

#include <string>
//#include <exception>
//#include <memory>
//#include <optional>

class NuSonicClientBase {
public:
  //constructor
  NuSonicClientBase(const fhicl::ParameterSet& params, const std::string& clientName);

  //destructor
  virtual ~NuSonicClientBase() = default;

  const std::string& clientName() const { return clientName_; }

  //main operation
  virtual void dispatch() {
    start();
    evaluate();
  }

  //helper: does nothing by default
  virtual void reset() {}

protected:
  virtual void evaluate() = 0;

  void start();

  void finish(bool success);

  //members
  unsigned allowedTries_, tries_;

  //for logging/debugging
  std::string clientName_;

};

#endif
