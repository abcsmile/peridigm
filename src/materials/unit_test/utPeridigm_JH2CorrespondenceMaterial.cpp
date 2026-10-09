//! \file utPeridigm_JH2CorrespondenceMaterial.cpp
//! \brief Tests the Peridigm DataManager wrapper for the JH2 kernel.

#include "Peridigm_JH2CorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "correspondence.h"

#include <Epetra_Map.h>
#include <Epetra_SerialComm.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <cmath>

namespace {

Teuchos::ParameterList graniteParameters() {
  Teuchos::ParameterList params;
  params.set("Material Model", "JH2 Correspondence");
  params.set("Density", 2660.0);
  params.set("Bulk Modulus", 21.9e9);
  params.set("Shear Modulus", 25.7e9);
  params.set("Hourglass Coefficient", 0.01);
  params.set("JH2 Intact Strength Constant", 1.9);
  params.set("JH2 Intact Strength Exponent", 0.62);
  params.set("JH2 Fractured Strength Constant", 0.25);
  params.set("JH2 Fractured Strength Exponent", 0.62);
  params.set("JH2 Strain Rate Coefficient", 0.005);
  params.set("JH2 Damage Constant D1", 0.005);
  params.set("JH2 Damage Exponent D2", 0.7);
  params.set("JH2 Energy Conversion Factor", 0.5);
  params.set("Pressure at HEL", 3.7e9);
  params.set("Strength at HEL", 1.2e9);
  params.set("JH2 Tensile Pressure", 54.0e6);
  params.set("JH2 Reference Strain Rate", 1.0);
  params.set("JH2 EOS K1", 21.9e9);
  params.set("JH2 EOS K2", -4500.0e9);
  params.set("JH2 EOS K3", 300000.0e9);
  return params;
}

} // anonymous namespace

TEUCHOS_UNIT_TEST(JH2CorrespondenceMaterial, DataManagerElasticUpdate) {
  PeridigmNS::JH2CorrespondenceMaterial material(graniteParameters());

  Epetra_SerialComm communicator;
  Epetra_Map scalarMap(1, 0, communicator);
  Epetra_Map vectorMap(3, 0, communicator);
  Epetra_Map emptyBondMap(0, 0, communicator);

  PeridigmNS::DataManager dataManager;
  dataManager.setMaps(Teuchos::rcp(&scalarMap, false),
                      Teuchos::rcp(&scalarMap, false),
                      Teuchos::rcp(&vectorMap, false),
                      Teuchos::rcp(&vectorMap, false),
                      Teuchos::rcp(&emptyBondMap, false));
  dataManager.allocateData(material.FieldIds());

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  const int rateFieldId = fieldManager.getFieldId("Unrotated_Rate_Of_Deformation");
  const int deformationGradientFieldId = fieldManager.getFieldId("Deformation_Gradient");
  const int stressFieldId = fieldManager.getFieldId("Unrotated_Cauchy_Stress");
  const int damageFieldId = fieldManager.getFieldId("JH2_Damage");
  const int plasticStrainFieldId = fieldManager.getFieldId("Equivalent_Plastic_Strain");
  const int bulkingPressureFieldId = fieldManager.getFieldId("JH2_Bulking_Pressure");
  const int strainRateFieldId = fieldManager.getFieldId("JH2_Equivalent_Strain_Rate");
  const int tensileFailureFieldId = fieldManager.getFieldId("JH2_Tensile_Failure");

  double *rate, *deformationGradient, *stressN, *stressNP1;
  double *damageN, *plasticStrainN, *bulkingPressureN, *strainRate;
  double *tensileFailureN;
  dataManager.getData(rateFieldId,
      PeridigmNS::PeridigmField::STEP_NONE)->ExtractView(&rate);
  dataManager.getData(deformationGradientFieldId,
      PeridigmNS::PeridigmField::STEP_NONE)->ExtractView(&deformationGradient);
  dataManager.getData(stressFieldId,
      PeridigmNS::PeridigmField::STEP_N)->ExtractView(&stressN);
  dataManager.getData(stressFieldId,
      PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&stressNP1);
  dataManager.getData(damageFieldId,
      PeridigmNS::PeridigmField::STEP_N)->ExtractView(&damageN);
  dataManager.getData(plasticStrainFieldId,
      PeridigmNS::PeridigmField::STEP_N)->ExtractView(&plasticStrainN);
  dataManager.getData(bulkingPressureFieldId,
      PeridigmNS::PeridigmField::STEP_N)->ExtractView(&bulkingPressureN);
  dataManager.getData(strainRateFieldId,
      PeridigmNS::PeridigmField::STEP_NONE)->ExtractView(&strainRate);
  dataManager.getData(tensileFailureFieldId,
      PeridigmNS::PeridigmField::STEP_N)->ExtractView(&tensileFailureN);

  for(int component = 0; component < 9; ++component) {
    rate[component] = 0.0;
    deformationGradient[component] = 0.0;
    stressN[component] = 0.0;
  }
  rate[0] = -500.0;
  deformationGradient[0] = 1.0;
  deformationGradient[4] = 1.0;
  deformationGradient[8] = 1.0;
  damageN[0] = 0.0;
  plasticStrainN[0] = 0.0;
  bulkingPressureN[0] = 0.0;
  tensileFailureN[0] = 0.0;

  material.computeCauchyStress(1.0e-6, 1, dataManager);

  const double expectedSigma11 =
      (-21.9e9 - 4.0 * 25.7e9 / 3.0) * 5.0e-4;
  const double expectedSigma22 =
      (-21.9e9 + 2.0 * 25.7e9 / 3.0) * 5.0e-4;
  TEST_FLOATING_EQUALITY(stressNP1[0], expectedSigma11, 1.0e-12);
  TEST_FLOATING_EQUALITY(stressNP1[4], expectedSigma22, 1.0e-12);
  TEST_FLOATING_EQUALITY(stressNP1[8], expectedSigma22, 1.0e-12);
  TEST_FLOATING_EQUALITY(strainRate[0], 333.3333333333333, 1.0e-12);
}

TEUCHOS_UNIT_TEST(JH2CorrespondenceMaterial, PreloadedRigidRotation) {
  PeridigmNS::JH2CorrespondenceMaterial material(graniteParameters());

  Epetra_SerialComm communicator;
  Epetra_Map scalarMap(1, 0, communicator);
  Epetra_Map vectorMap(3, 0, communicator);
  Epetra_Map emptyBondMap(0, 0, communicator);

  PeridigmNS::DataManager dataManager;
  dataManager.setMaps(Teuchos::rcp(&scalarMap, false),
                      Teuchos::rcp(&scalarMap, false),
                      Teuchos::rcp(&vectorMap, false),
                      Teuchos::rcp(&vectorMap, false),
                      Teuchos::rcp(&emptyBondMap, false));
  dataManager.allocateData(material.FieldIds());

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  const int rateFieldId = fieldManager.getFieldId(
      "Unrotated_Rate_Of_Deformation");
  const int deformationGradientFieldId = fieldManager.getFieldId(
      "Deformation_Gradient");
  const int stressFieldId = fieldManager.getFieldId(
      "Unrotated_Cauchy_Stress");
  const int damageFieldId = fieldManager.getFieldId("JH2_Damage");
  const int plasticStrainFieldId = fieldManager.getFieldId(
      "Equivalent_Plastic_Strain");
  const int bulkingPressureFieldId = fieldManager.getFieldId(
      "JH2_Bulking_Pressure");
  const int tensileFailureFieldId = fieldManager.getFieldId(
      "JH2_Tensile_Failure");

  double *rate, *deformationGradient, *stressN, *stressNP1;
  double *damageN, *damageNP1, *plasticStrainN, *plasticStrainNP1;
  double *bulkingPressureN, *bulkingPressureNP1;
  double *tensileFailureN, *tensileFailureNP1;
  dataManager.getData(rateFieldId,
      PeridigmNS::PeridigmField::STEP_NONE)->ExtractView(&rate);
  dataManager.getData(deformationGradientFieldId,
      PeridigmNS::PeridigmField::STEP_NONE)->ExtractView(&deformationGradient);
  dataManager.getData(stressFieldId,
      PeridigmNS::PeridigmField::STEP_N)->ExtractView(&stressN);
  dataManager.getData(stressFieldId,
      PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&stressNP1);
  dataManager.getData(damageFieldId,
      PeridigmNS::PeridigmField::STEP_N)->ExtractView(&damageN);
  dataManager.getData(damageFieldId,
      PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&damageNP1);
  dataManager.getData(plasticStrainFieldId,
      PeridigmNS::PeridigmField::STEP_N)->ExtractView(&plasticStrainN);
  dataManager.getData(plasticStrainFieldId,
      PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&plasticStrainNP1);
  dataManager.getData(bulkingPressureFieldId,
      PeridigmNS::PeridigmField::STEP_N)->ExtractView(&bulkingPressureN);
  dataManager.getData(bulkingPressureFieldId,
      PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&bulkingPressureNP1);
  dataManager.getData(tensileFailureFieldId,
      PeridigmNS::PeridigmField::STEP_N)->ExtractView(&tensileFailureN);
  dataManager.getData(tensileFailureFieldId,
      PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&tensileFailureNP1);

  const double angle = 3.14159265358979323846 / 6.0;
  const double cosine = std::cos(angle);
  const double sine = std::sin(angle);
  const double rotation[9] = {
      cosine, -sine, 0.0,
      sine, cosine, 0.0,
      0.0, 0.0, 1.0
  };
  for(int component = 0; component < 9; ++component) {
    rate[component] = 0.0;
    deformationGradient[component] = rotation[component];
    stressN[component] = 0.0;
  }
  stressN[0] = 10.0e6;
  stressN[4] = -5.0e6;
  stressN[8] = 2.0e6;
  damageN[0] = 0.15;
  plasticStrainN[0] = 0.002;
  bulkingPressureN[0] = 1.0e6;
  tensileFailureN[0] = 0.0;

  material.computeCauchyStress(1.0e-6, 1, dataManager);

  double rotatedStress[9] = {};
  CORRESPONDENCE::rotateCauchyStress(
      rotation, stressNP1, rotatedStress, 1);

  const double expectedXX =
      cosine * cosine * stressN[0] + sine * sine * stressN[4];
  const double expectedYY =
      sine * sine * stressN[0] + cosine * cosine * stressN[4];
  const double expectedXY = cosine * sine * (stressN[0] - stressN[4]);
  TEST_FLOATING_EQUALITY(stressNP1[0], stressN[0], 1.0e-13);
  TEST_FLOATING_EQUALITY(stressNP1[4], stressN[4], 1.0e-13);
  TEST_FLOATING_EQUALITY(stressNP1[8], stressN[8], 1.0e-13);
  TEST_FLOATING_EQUALITY(rotatedStress[0], expectedXX, 1.0e-13);
  TEST_FLOATING_EQUALITY(rotatedStress[1], expectedXY, 1.0e-13);
  TEST_FLOATING_EQUALITY(rotatedStress[3], expectedXY, 1.0e-13);
  TEST_FLOATING_EQUALITY(rotatedStress[4], expectedYY, 1.0e-13);
  TEST_FLOATING_EQUALITY(rotatedStress[8], stressN[8], 1.0e-13);
  TEST_FLOATING_EQUALITY(damageNP1[0], damageN[0], 1.0e-13);
  TEST_FLOATING_EQUALITY(plasticStrainNP1[0], plasticStrainN[0], 1.0e-13);
  TEST_FLOATING_EQUALITY(
      bulkingPressureNP1[0], bulkingPressureN[0], 1.0e-13);
  TEST_FLOATING_EQUALITY(tensileFailureNP1[0], 0.0, 1.0e-13);
}

int main(int argc, char* argv[]) {
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
