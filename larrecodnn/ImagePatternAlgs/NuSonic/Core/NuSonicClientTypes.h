#ifndef NuSonic_Core_SonicClientTypes
#define NuSonic_Core_SonicClientTypes

//this base class exists to limit the impact of dependent scope in derived classes
template <typename InputT, typename OutputT = InputT>
class NuSonicClientTypes {
public:
  //typedefs for outside accessibility
  typedef InputT Input;
  typedef OutputT Output;
  //destructor
  virtual ~NuSonicClientTypes() = default;

  //accessors
  Input& input() { return input_; }
  const Output& output() const { return output_; }

protected:
  Input input_;
  Output output_;
};

#endif
