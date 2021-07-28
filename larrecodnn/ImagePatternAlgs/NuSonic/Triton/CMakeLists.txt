cet_find_library(GRPC_CLIENT NAMES grpcclient PATHS ENV TRITON_LIB REQUIRED)

art_make(
          LIB_LIBRARIES
          fhiclcpp::fhiclcpp
          cetlib::cetlib cetlib_except
          ROOT::Core
          canvas::canvas
          messagefacility::MF_MessageLogger
	  ${GRPC_CLIENT}
          ${PROTOBUF}
        )

install_headers()
install_source()
