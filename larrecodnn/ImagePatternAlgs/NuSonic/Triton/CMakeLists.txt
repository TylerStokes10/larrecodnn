cet_find_library(GRPC_CLIENT NAMES grpcclient PATHS ENV TRITON_LIB REQUIRED)

art_make(
          LIB_LIBRARIES
          ${FHICLCPP}
          cetlib cetlib_except
          ROOT::Core
          canvas
          ${MF_MESSAGELOGGER}
	  ${GRPC_CLIENT}
          ${PROTOBUF}
        )

install_headers()
install_source()
