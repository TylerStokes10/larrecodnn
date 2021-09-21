art_make(
          LIB_LIBRARIES
          fhiclcpp::fhiclcpp
          cetlib::cetlib cetlib_except
          ROOT::Core
          canvas::canvas
          messagefacility::MF_MessageLogger
          TRITON::grpcclient
        )

install_headers()
install_source()
