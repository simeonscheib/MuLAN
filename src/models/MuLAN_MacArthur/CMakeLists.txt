# Add the model target
add_model(MuLAN_MA MuLAN_MA.cc)
# NOTE The target should have the same name as the model folder and the *.cc

add_unit_test(NAME integrator_test
              GROUP model_mulan_ma
              SOURCES tests/integrator_test.cc)

add_unit_test(NAME domain_test
              GROUP model_mulan_ma
              SOURCES tests/domain_test.cc)

add_unit_test(NAME helper_test
              GROUP model_mulan_ma
              SOURCES tests/helper_test.cc)

add_unit_test(NAME Organism_Manager_test
        GROUP model_mulan_ma
        SOURCES tests/Organism_Manager_test.cc)
