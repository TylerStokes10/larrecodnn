#ifndef NuSonic_Core_NuSonicClient
#define NuSonic_Core_NuSonicClient

#include "larrecodnn/ImagePatternAlgs/NuSonic/Core/NuSonicClientBase.h"
#include "larrecodnn/ImagePatternAlgs/NuSonic/Core/NuSonicClientTypes.h"

//convenience definition for multiple inheritance (base and types)
template <typename InputT, typename OutputT = InputT>
class NuSonicClient : public NuSonicClientBase, public NuSonicClientTypes<InputT, OutputT> {
public:
  //constructor
  NuSonicClient(const fhicl::ParameterSet& params, const std::string& clientName)
      : NuSonicClientBase(params, clientName), NuSonicClientTypes<InputT, OutputT>() {}
};

#endif
