include_directories($ENV{TRITON_INC})
include_directories($ENV{GRPC_INC})
cet_find_library(TRITON_LIB NAMES grpcclient PATHS $ENV{TRITON_LIB})

art_make(
          TOOL_LIBRARIES
          larreco_RecoAlg_ImagePatternAlgs_DataProvider
          ROOT::Core
          canvas
          ${MF_MESSAGELOGGER}
	  ${TRITON_LIB}
        )

install_source()
