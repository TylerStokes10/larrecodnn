# export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/uboone/app/users/cerati/products/torch_scatter_v2.1.2-e26
# export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/uboone/app/users/cerati/products/pytorch-install-v2.1.0-e26

# unsetup libtorch

find_package(TorchScatter REQUIRED)
include_directories(SYSTEM /uboone/app/users/cerati/products/torch_scatter_v2.1.2-e26/include/)
link_directories(/uboone/app/users/cerati/products/torch_scatter_v2.1.2-e26/lib64)

include_directories(SYSTEM /uboone/app/users/cerati/products/pytorch-install-v2.1.0-e26/include)
link_directories(/uboone/app/users/cerati/products/pytorch-install-v2.1.0-e26/lib)

cet_build_plugin(NuGraphInference art::EDProducer
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
    lardataobj::RecoBase
    torch
    torchscatter
)

cet_build_plugin(PandoraNuSliceHitsProducer art::EDProducer
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
)

install_headers()
install_fhicl()
install_source()
