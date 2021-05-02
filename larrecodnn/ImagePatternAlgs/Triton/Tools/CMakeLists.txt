include_directories($ENV{TRITON_INC})
cet_find_library(TRITON_LIB NAMES request PATHS $ENV{TRITON_LIB})

art_make(
          TOOL_LIBRARIES
          larreco_RecoAlg_ImagePatternAlgs_DataProvider
          ROOT::Core
          canvas
          ${MF_MESSAGELOGGER}
	  ${TRITON_LIB}
        )

install_source()
