//! \file utPeridigm_JH2CorrespondenceMaterial.cpp
//! \brief Tests the Peridigm DataManager wrapper for the JH2 kernel.

#include "Peridigm_JH2CorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "correspondence.h"

#include <Epetra_BlockMap.h>
#include <Epetra_Map.h>
#include <Epetra_SerialComm.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <algorithm>
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

void setFourPointReferenceFields(PeridigmNS::DataManager& dataManager) {
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  const int horizonFieldId = fieldManager.getFieldId("Horizon");
  const int volumeFieldId = fieldManager.getFieldId("Volume");
  const int modelCoordinatesFieldId =
      fieldManager.getFieldId("Model_Coordinates");
  const int coordinatesFieldId = fieldManager.getFieldId("Coordinates");
  const int velocityFieldId = fieldManager.getFieldId("Velocity");

  double *horizon, *volume, *modelCoordinates;
  double *coordinatesN, *coordinatesNP1, *velocityN, *velocityNP1;
  dataManager.getData(horizonFieldId,
      PeridigmNS::PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(volumeFieldId,
      PeridigmNS::PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(modelCoordinatesFieldId,
      PeridigmNS::PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(coordinatesFieldId,
      PeridigmNS::PeridigmField::STEP_N)->ExtractView(&coordinatesN);
  dataManager.getData(coordinatesFieldId,
      PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&coordinatesNP1);
  dataManager.getData(velocityFieldId,
      PeridigmNS::PeridigmField::STEP_N)->ExtractView(&velocityN);
  dataManager.getData(velocityFieldId,
      PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&velocityNP1);

  const double referenceCoordinates[12] = {
      0.0, 0.0, 0.0,
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0
  };
  for(int point = 0; point < 4; ++point) {
    horizon[point] = 2.0;
    volume[point] = 1.0;
    for(int component = 0; component < 3; ++component) {
      const int index = 3 * point + component;
      modelCoordinates[index] = referenceCoordinates[index];
      coordinatesN[index] = referenceCoordinates[index];
      coordinatesNP1[index] = referenceCoordinates[index];
      velocityN[index] = 0.0;
      velocityNP1[index] = 0.0;
    }
  }
}

void setFourPointAffineKinematics(PeridigmNS::DataManager& dataManager,
                                  const double* deformationGradient,
                                  const double* velocityGradient) {
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  const int modelCoordinatesFieldId =
      fieldManager.getFieldId("Model_Coordinates");
  const int coordinatesFieldId = fieldManager.getFieldId("Coordinates");
  const int velocityFieldId = fieldManager.getFieldId("Velocity");
  double *modelCoordinates, *coordinatesNP1, *velocityNP1;
  dataManager.getData(modelCoordinatesFieldId,
      PeridigmNS::PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(coordinatesFieldId,
      PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&coordinatesNP1);
  dataManager.getData(velocityFieldId,
      PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&velocityNP1);

  for(int point = 0; point < 4; ++point) {
    const double* referencePosition = modelCoordinates + 3 * point;
    double* currentPosition = coordinatesNP1 + 3 * point;
    double* currentVelocity = velocityNP1 + 3 * point;
    for(int row = 0; row < 3; ++row) {
      currentPosition[row] = 0.0;
      for(int column = 0; column < 3; ++column)
        currentPosition[row] +=
            deformationGradient[3 * row + column] * referencePosition[column];
    }
    for(int row = 0; row < 3; ++row) {
      currentVelocity[row] = 0.0;
      for(int column = 0; column < 3; ++column)
        currentVelocity[row] +=
            velocityGradient[3 * row + column] * currentPosition[column];
    }
  }
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

TEUCHOS_UNIT_TEST(JH2CorrespondenceMaterial, FourPointAffineCompression) {
  PeridigmNS::JH2CorrespondenceMaterial material(graniteParameters());

  const int numPoints = 4;
  const int ownedIds[numPoints] = {0, 1, 2, 3};
  const int bondSizes[numPoints] = {3, 3, 3, 3};
  const int neighborhoodList[16] = {
      3, 1, 2, 3,
      3, 0, 2, 3,
      3, 0, 1, 3,
      3, 0, 1, 2
  };

  Epetra_SerialComm communicator;
  Epetra_BlockMap scalarMap(numPoints, 1, 0, communicator);
  Epetra_BlockMap vectorMap(numPoints, 3, 0, communicator);
  Epetra_BlockMap bondMap(
      numPoints, numPoints, ownedIds, bondSizes, 0, communicator);

  PeridigmNS::DataManager dataManager;
  dataManager.setMaps(Teuchos::rcp(&scalarMap, false),
                      Teuchos::rcp(&scalarMap, false),
                      Teuchos::rcp(&vectorMap, false),
                      Teuchos::rcp(&vectorMap, false),
                      Teuchos::rcp(&bondMap, false));
  dataManager.allocateData(material.FieldIds());

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  const int horizonFieldId = fieldManager.getFieldId("Horizon");
  const int volumeFieldId = fieldManager.getFieldId("Volume");
  const int modelCoordinatesFieldId =
      fieldManager.getFieldId("Model_Coordinates");
  const int coordinatesFieldId = fieldManager.getFieldId("Coordinates");
  const int velocityFieldId = fieldManager.getFieldId("Velocity");
  const int deformationGradientFieldId =
      fieldManager.getFieldId("Deformation_Gradient");
  const int rateFieldId = fieldManager.getFieldId(
      "Unrotated_Rate_Of_Deformation");
  const int rotationFieldId = fieldManager.getFieldId("Rotation_Tensor");
  const int stressFieldId = fieldManager.getFieldId(
      "Unrotated_Cauchy_Stress");
  const int forceFieldId = fieldManager.getFieldId("Force_Density");
  const int hourglassForceFieldId =
      fieldManager.getFieldId("Hourglass_Force_Density");
  const int jacobianFieldId =
      fieldManager.getFieldId("JH2_Deformation_Jacobian");
  const int damageFieldId = fieldManager.getFieldId("JH2_Damage");

  double *horizon, *volume, *modelCoordinates;
  double *coordinatesN, *coordinatesNP1, *velocityN, *velocityNP1;
  dataManager.getData(horizonFieldId,
      PeridigmNS::PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(volumeFieldId,
      PeridigmNS::PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(modelCoordinatesFieldId,
      PeridigmNS::PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(coordinatesFieldId,
      PeridigmNS::PeridigmField::STEP_N)->ExtractView(&coordinatesN);
  dataManager.getData(coordinatesFieldId,
      PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&coordinatesNP1);
  dataManager.getData(velocityFieldId,
      PeridigmNS::PeridigmField::STEP_N)->ExtractView(&velocityN);
  dataManager.getData(velocityFieldId,
      PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&velocityNP1);

  const double referenceCoordinates[3 * numPoints] = {
      0.0, 0.0, 0.0,
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0
  };
  for(int point = 0; point < numPoints; ++point) {
    horizon[point] = 2.0;
    volume[point] = 1.0;
    for(int component = 0; component < 3; ++component) {
      const int index = 3 * point + component;
      modelCoordinates[index] = referenceCoordinates[index];
      coordinatesN[index] = referenceCoordinates[index];
      coordinatesNP1[index] = referenceCoordinates[index];
      velocityN[index] = 0.0;
      velocityNP1[index] = 0.0;
    }
  }

  const double dt = 1.0e-6;
  material.initialize(
      dt, numPoints, ownedIds, neighborhoodList, dataManager);

  // Homogeneous spatial compression: F = diag(0.9995, 1, 1) and
  // L = diag(-500, 0, 0) [1/s].  Assign v = L y so Fdot F^-1 = L.
  const double axialStretch = 0.9995;
  for(int point = 0; point < numPoints; ++point) {
    coordinatesNP1[3 * point] =
        axialStretch * modelCoordinates[3 * point];
    coordinatesNP1[3 * point + 1] = modelCoordinates[3 * point + 1];
    coordinatesNP1[3 * point + 2] = modelCoordinates[3 * point + 2];
    velocityNP1[3 * point] = -500.0 * coordinatesNP1[3 * point];
    velocityNP1[3 * point + 1] = 0.0;
    velocityNP1[3 * point + 2] = 0.0;
  }

  material.computeForce(
      dt, numPoints, ownedIds, neighborhoodList, dataManager);

  double *deformationGradient, *rate, *rotation, *stressNP1;
  double *forceNP1, *hourglassForceNP1, *jacobian, *damageNP1;
  dataManager.getData(deformationGradientFieldId,
      PeridigmNS::PeridigmField::STEP_NONE)->ExtractView(
          &deformationGradient);
  dataManager.getData(rateFieldId,
      PeridigmNS::PeridigmField::STEP_NONE)->ExtractView(&rate);
  dataManager.getData(rotationFieldId,
      PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&rotation);
  dataManager.getData(stressFieldId,
      PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&stressNP1);
  dataManager.getData(forceFieldId,
      PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&forceNP1);
  dataManager.getData(hourglassForceFieldId,
      PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&hourglassForceNP1);
  dataManager.getData(jacobianFieldId,
      PeridigmNS::PeridigmField::STEP_NONE)->ExtractView(&jacobian);
  dataManager.getData(damageFieldId,
      PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&damageNP1);

  const double expectedStress[9] = {
      (-21.9e9 - 4.0 * 25.7e9 / 3.0) * 5.0e-4, 0.0, 0.0,
      0.0, (-21.9e9 + 2.0 * 25.7e9 / 3.0) * 5.0e-4, 0.0,
      0.0, 0.0, (-21.9e9 + 2.0 * 25.7e9 / 3.0) * 5.0e-4
  };
  for(int point = 0; point < numPoints; ++point) {
    for(int component = 0; component < 9; ++component) {
      const double expectedF = component == 0 ? axialStretch :
          (component == 4 || component == 8 ? 1.0 : 0.0);
      const double expectedRate = component == 0 ? -500.0 : 0.0;
      const double expectedRotation =
          (component == 0 || component == 4 || component == 8) ? 1.0 : 0.0;
      const double actualF = deformationGradient[9 * point + component];
      const double actualRate = rate[9 * point + component];
      const double actualRotation = rotation[9 * point + component];
      const double actualStress = stressNP1[9 * point + component];
      if(expectedF == 0.0) {
        TEST_ASSERT(std::abs(actualF) <= 1.0e-14);
      }
      else {
        TEST_FLOATING_EQUALITY(actualF, expectedF, 1.0e-12);
      }
      if(expectedRate == 0.0) {
        TEST_ASSERT(std::abs(actualRate) <= 1.0e-10);
      }
      else {
        TEST_FLOATING_EQUALITY(actualRate, expectedRate, 1.0e-11);
      }
      if(expectedRotation == 0.0) {
        TEST_ASSERT(std::abs(actualRotation) <= 1.0e-14);
      }
      else {
        TEST_FLOATING_EQUALITY(actualRotation, expectedRotation, 1.0e-12);
      }
      if(expectedStress[component] == 0.0) {
        TEST_ASSERT(std::abs(actualStress) <= 1.0e-6);
      }
      else {
        TEST_FLOATING_EQUALITY(
            actualStress, expectedStress[component], 1.0e-10);
      }
    }
    TEST_FLOATING_EQUALITY(jacobian[point], axialStretch, 1.0e-12);
    TEST_FLOATING_EQUALITY(damageNP1[point], 0.0, 1.0e-13);
  }

  double forceBalance[3] = {};
  double maximumForce = 0.0;
  double maximumHourglassForce = 0.0;
  for(int point = 0; point < numPoints; ++point) {
    for(int component = 0; component < 3; ++component) {
      const int index = 3 * point + component;
      forceBalance[component] += forceNP1[index];
      maximumForce = std::max(maximumForce, std::abs(forceNP1[index]));
      maximumHourglassForce = std::max(
          maximumHourglassForce, std::abs(hourglassForceNP1[index]));
    }
  }
  const double forceTolerance = 1.0e-12 * std::max(1.0, maximumForce);
  TEST_ASSERT(std::abs(forceBalance[0]) <= forceTolerance);
  TEST_ASSERT(std::abs(forceBalance[1]) <= forceTolerance);
  TEST_ASSERT(std::abs(forceBalance[2]) <= forceTolerance);
  TEST_ASSERT(maximumHourglassForce <= forceTolerance);
}

TEUCHOS_UNIT_TEST(JH2CorrespondenceMaterial,
                  MultiStepStretchAndRotation) {
  PeridigmNS::JH2CorrespondenceMaterial material(graniteParameters());

  const int numPoints = 4;
  const int ownedIds[numPoints] = {0, 1, 2, 3};
  const int bondSizes[numPoints] = {3, 3, 3, 3};
  const int neighborhoodList[16] = {
      3, 1, 2, 3,
      3, 0, 2, 3,
      3, 0, 1, 3,
      3, 0, 1, 2
  };
  Epetra_SerialComm communicator;
  Epetra_BlockMap scalarMap(numPoints, 1, 0, communicator);
  Epetra_BlockMap vectorMap(numPoints, 3, 0, communicator);
  Epetra_BlockMap bondMap(
      numPoints, numPoints, ownedIds, bondSizes, 0, communicator);
  PeridigmNS::DataManager dataManager;
  dataManager.setMaps(Teuchos::rcp(&scalarMap, false),
                      Teuchos::rcp(&scalarMap, false),
                      Teuchos::rcp(&vectorMap, false),
                      Teuchos::rcp(&vectorMap, false),
                      Teuchos::rcp(&bondMap, false));
  dataManager.allocateData(material.FieldIds());
  setFourPointReferenceFields(dataManager);

  const double dt = 1.0e-7;
  const double compressionRate = 100.0;
  const double angularVelocity = 1.0e4;
  const int numSteps = 100;
  material.initialize(
      dt, numPoints, ownedIds, neighborhoodList, dataManager);

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  const int deformationGradientFieldId =
      fieldManager.getFieldId("Deformation_Gradient");
  const int rateFieldId = fieldManager.getFieldId(
      "Unrotated_Rate_Of_Deformation");
  const int rotationFieldId = fieldManager.getFieldId("Rotation_Tensor");
  const int unrotatedStressFieldId = fieldManager.getFieldId(
      "Unrotated_Cauchy_Stress");
  const int spatialStressFieldId = fieldManager.getFieldId("Cauchy_Stress");
  const int forceFieldId = fieldManager.getFieldId("Force_Density");
  const int hourglassForceFieldId =
      fieldManager.getFieldId("Hourglass_Force_Density");
  const int damageFieldId = fieldManager.getFieldId("JH2_Damage");
  const int jacobianFieldId =
      fieldManager.getFieldId("JH2_Deformation_Jacobian");

  for(int step = 1; step <= numSteps; ++step) {
    const double time = step * dt;
    const double angle = angularVelocity * time;
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    const double stretch = std::exp(-compressionRate * time);
    const double prescribedF[9] = {
        cosine * stretch, -sine, 0.0,
        sine * stretch, cosine, 0.0,
        0.0, 0.0, 1.0
    };
    const double spatialDxx = -compressionRate * cosine * cosine;
    const double spatialDxy = -compressionRate * cosine * sine;
    const double spatialDyy = -compressionRate * sine * sine;
    const double velocityGradient[9] = {
        spatialDxx, spatialDxy - angularVelocity, 0.0,
        spatialDxy + angularVelocity, spatialDyy, 0.0,
        0.0, 0.0, 0.0
    };
    setFourPointAffineKinematics(
        dataManager, prescribedF, velocityGradient);
    material.computeForce(
        dt, numPoints, ownedIds, neighborhoodList, dataManager);

    double *computedF, *rate, *rotation, *unrotatedStress, *spatialStress;
    double *force, *hourglassForce, *damage, *jacobian;
    dataManager.getData(deformationGradientFieldId,
        PeridigmNS::PeridigmField::STEP_NONE)->ExtractView(&computedF);
    dataManager.getData(rateFieldId,
        PeridigmNS::PeridigmField::STEP_NONE)->ExtractView(&rate);
    dataManager.getData(rotationFieldId,
        PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&rotation);
    dataManager.getData(unrotatedStressFieldId,
        PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&unrotatedStress);
    dataManager.getData(spatialStressFieldId,
        PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&spatialStress);
    dataManager.getData(forceFieldId,
        PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&force);
    dataManager.getData(hourglassForceFieldId,
        PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&hourglassForce);
    dataManager.getData(damageFieldId,
        PeridigmNS::PeridigmField::STEP_NP1)->ExtractView(&damage);
    dataManager.getData(jacobianFieldId,
        PeridigmNS::PeridigmField::STEP_NONE)->ExtractView(&jacobian);

    const double strainMagnitude = compressionRate * time;
    const double axialStress =
        -(21.9e9 + 4.0 * 25.7e9 / 3.0) * strainMagnitude;
    const double lateralStress =
        -(21.9e9 - 2.0 * 25.7e9 / 3.0) * strainMagnitude;
    const double expectedSpatialXX = cosine * cosine * axialStress +
        sine * sine * lateralStress;
    const double expectedSpatialYY = sine * sine * axialStress +
        cosine * cosine * lateralStress;
    const double expectedSpatialXY =
        cosine * sine * (axialStress - lateralStress);

    double forceBalance[3] = {};
    double maximumForce = 0.0;
    double maximumHourglassForce = 0.0;
    for(int point = 0; point < numPoints; ++point) {
      for(int component = 0; component < 9; ++component) {
        const double actualF = computedF[9 * point + component];
        if(prescribedF[component] == 0.0) {
          TEST_ASSERT(std::abs(actualF) <= 1.0e-14);
        }
        else {
          TEST_FLOATING_EQUALITY(
              actualF, prescribedF[component], 1.0e-12);
        }
      }
      TEST_FLOATING_EQUALITY(rotation[9 * point], cosine, 2.0e-4);
      TEST_FLOATING_EQUALITY(rotation[9 * point + 1], -sine, 2.0e-4);
      TEST_FLOATING_EQUALITY(rotation[9 * point + 3], sine, 2.0e-4);
      TEST_FLOATING_EQUALITY(rotation[9 * point + 4], cosine, 2.0e-4);
      TEST_FLOATING_EQUALITY(rate[9 * point], -compressionRate, 3.0e-4);
      TEST_ASSERT(std::abs(rate[9 * point + 1]) <= 2.0e-2);
      TEST_ASSERT(std::abs(rate[9 * point + 3]) <= 2.0e-2);
      TEST_FLOATING_EQUALITY(unrotatedStress[9 * point], axialStress, 3.0e-4);
      TEST_FLOATING_EQUALITY(
          unrotatedStress[9 * point + 4], lateralStress, 3.0e-4);
      TEST_FLOATING_EQUALITY(spatialStress[9 * point], expectedSpatialXX, 5.0e-4);
      TEST_FLOATING_EQUALITY(
          spatialStress[9 * point + 1], expectedSpatialXY, 5.0e-4);
      TEST_FLOATING_EQUALITY(
          spatialStress[9 * point + 3], expectedSpatialXY, 5.0e-4);
      TEST_FLOATING_EQUALITY(
          spatialStress[9 * point + 4], expectedSpatialYY, 5.0e-4);
      TEST_FLOATING_EQUALITY(
          jacobian[point], stretch, 1.0e-12);
      TEST_FLOATING_EQUALITY(damage[point], 0.0, 1.0e-13);
      for(int component = 0; component < 3; ++component) {
        const int index = 3 * point + component;
        forceBalance[component] += force[index];
        maximumForce = std::max(maximumForce, std::abs(force[index]));
        maximumHourglassForce = std::max(
            maximumHourglassForce, std::abs(hourglassForce[index]));
      }
    }
    const double forceTolerance = 1.0e-12 * std::max(1.0, maximumForce);
    TEST_ASSERT(std::abs(forceBalance[0]) <= forceTolerance);
    TEST_ASSERT(std::abs(forceBalance[1]) <= forceTolerance);
    TEST_ASSERT(std::abs(forceBalance[2]) <= forceTolerance);
    TEST_ASSERT(maximumHourglassForce <= forceTolerance);

    if(step < numSteps)
      dataManager.updateState();
  }
}

int main(int argc, char* argv[]) {
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
