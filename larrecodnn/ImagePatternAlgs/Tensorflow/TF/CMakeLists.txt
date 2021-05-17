include_directories( $ENV{TENSORFLOW_INC}/absl )

# the code in tf::Graph::~Graph() should be inspected 
cet_add_compiler_flags(CXX -Wno-unused-result)

art_make(
          LIB_LIBRARIES
			${FHICLCPP}
			cetlib cetlib_except
			${TENSORFLOW}
			${PROTOBUF}
			${FOR_TENSORFLOW}
      #rt #workaround for larsoft issue 17960
        )

install_headers()
install_source()

