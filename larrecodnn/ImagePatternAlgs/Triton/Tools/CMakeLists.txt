
cet_find_library(GRPC_CLIENT NAMES grpcclient PATHS ENV TRITON_LIB REQUIRED)

art_make(
          TOOL_LIBRARIES
          larreco_RecoAlg_ImagePatternAlgs_DataProvider
          larcorealg_Geometry
          larrecodnn_ImagePatternAlgs_NuSonic_Triton
          ROOT::Core
          canvas::canvas
          messagefacility::MF_MessageLogger
	  ${GRPC_CLIENT}
          ${PROTOBUF}
        )

install_source()
