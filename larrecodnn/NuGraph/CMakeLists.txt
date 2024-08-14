
cet_build_plugin(NuGraphInference art::EDProducer
    LIBRARIES PRIVATE
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_BASIC}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_PERSISTENCY_PROVENANCE}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    ${ROOT_CORE}
    ${ROOT_MATH}
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${ART_PERSISTENCY_COMMON}
    ${ART_UTILITIES}
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    ${ROOT_BASIC_LIB_LIST}
    lardataobj::RecoBase
    ${TORCH_LIBRARIES}
    TorchScatter::TorchScatter
    larcorealg::Geometry
    larcore::Geometry_Geometry_service
    IMPL_TARGET_VAR NuGraphInference_module
)

#set(DICT_INSTALL_DIR "/exp/dune/data/users/tstokes/DUNE_Software/build_slf7.x86_64/dict")
#message(STATUS "DICT_INSTALL_DIR: ${DICT_INSTALL_DIR}")

# Add MyDict.cxx to the sources
#target_sources(${NuGraphInference_module} PRIVATE MyDict.cxx)

#message(STATUS "CMAKE_INSTALL_DATAROOTDIR: ${CMAKE_INSTALL_DATAROOTDIR}")


find_path(DELAUNATOR_INC delaunator.hpp HINTS ENV DELAUNATOR_INC REQUIRED DOC "Header location for delaunator_cpp")
target_include_directories(${NuGraphInference_module} PRIVATE ${DELAUNATOR_INC})

cet_build_plugin(NuSliceHitsProducer art::EDProducer
    LIBRARIES PRIVATE
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_BASIC}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_PERSISTENCY_PROVENANCE}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    ${ROOT_CORE}
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${ART_PERSISTENCY_COMMON}
    ${ART_UTILITIES}
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    ${ROOT_BASIC_LIB_LIST}
    nusimdata::SimulationBase
    lardataobj::RecoBase
)

cet_build_plugin(NuGraphAnalyzer art::EDAnalyzer
    LIBRARIES PRIVATE
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_BASIC}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_PERSISTENCY_PROVENANCE}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    ${ROOT_CORE}
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${ART_PERSISTENCY_COMMON}
    ${ART_UTILITIES}
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    ${ROOT_BASIC_LIB_LIST}
    art_root_io::TFileService_service
    ROOT::Tree
    lardata::RecoBaseProxy
)


install_headers()
install_fhicl()
install_source()

# Ensure classes_def.xml is installed
#install(FILES classes_def.xml DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/dict)
#install(FILES classes_def.xml DESTINATION ${DICT_INSTALL_DIR})

