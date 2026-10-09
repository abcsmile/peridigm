//! \file Peridigm_JH2CorrespondenceMaterial.cpp
//! \brief JH2 constitutive model coupled to Peridigm correspondence kinematics.

#include "Peridigm_JH2CorrespondenceMaterial.hpp"

#include "Peridigm_Field.hpp"
#include <Teuchos_Assert.hpp>

#include <cmath>
#include <sstream>
#include <string>

using namespace std;

PeridigmNS::JH2CorrespondenceMaterial::JH2CorrespondenceMaterial(
    const Teuchos::ParameterList& params)
  : CorrespondenceMaterial(params),
    m_jh2DamageFieldId(-1),
    m_equivalentPlasticStrainFieldId(-1),
    m_jh2BulkingPressureFieldId(-1),
    m_jh2EquivalentStressFieldId(-1),
    m_jh2YieldStrengthFieldId(-1),
    m_jh2EOSPressureFieldId(-1),
    m_jh2EquivalentStrainRateFieldId(-1),
    m_jh2DeformationJacobianFieldId(-1),
    m_jh2PlasticStrainIncrementFieldId(-1),
    m_jh2FracturePlasticStrainFieldId(-1),
    m_jh2DissipatedEnergyFieldId(-1),
    m_jh2TensileFailureFieldId(-1)
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      params.isParameter("Thermal Expansion Coefficient"),
      "**** Error: Thermal expansion is not supported by JH2 Correspondence.");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      !std::isfinite(m_density) || m_density <= 0.0,
      "**** Error: JH2 Correspondence density must be positive and finite.");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      !std::isfinite(m_hourglassCoefficient) || m_hourglassCoefficient < 0.0,
      "**** Error: JH2 Correspondence hourglass coefficient must be nonnegative and finite.");

  m_jh2Parameters.bulkModulus = m_bulkModulus;
  m_jh2Parameters.shearModulus = m_shearModulus;
  m_jh2Parameters.intactStrengthConstant =
      params.get<double>("JH2 Intact Strength Constant");
  m_jh2Parameters.intactStrengthExponent =
      params.get<double>("JH2 Intact Strength Exponent");
  m_jh2Parameters.fracturedStrengthConstant =
      params.get<double>("JH2 Fractured Strength Constant");
  m_jh2Parameters.fracturedStrengthExponent =
      params.get<double>("JH2 Fractured Strength Exponent");
  m_jh2Parameters.strainRateCoefficient =
      params.get<double>("JH2 Strain Rate Coefficient");
  m_jh2Parameters.damageConstantD1 =
      params.get<double>("JH2 Damage Constant D1");
  m_jh2Parameters.damageExponentD2 =
      params.get<double>("JH2 Damage Exponent D2");
  m_jh2Parameters.energyConversionFactor =
      params.get<double>("JH2 Energy Conversion Factor");
  m_jh2Parameters.pressureAtHEL = params.get<double>("Pressure at HEL");
  m_jh2Parameters.strengthAtHEL = params.get<double>("Strength at HEL");
  m_jh2Parameters.tensilePressure = params.get<double>("JH2 Tensile Pressure");
  m_jh2Parameters.referenceStrainRate =
      params.get<double>("JH2 Reference Strain Rate");
  m_jh2Parameters.eosK1 = params.get<double>("JH2 EOS K1");
  m_jh2Parameters.eosK2 = params.get<double>("JH2 EOS K2");
  m_jh2Parameters.eosK3 = params.get<double>("JH2 EOS K3");

  string validationMessage;
  const bool validParameters =
      JH2::validateParameters(m_jh2Parameters, validationMessage);
  TEUCHOS_TEST_FOR_EXCEPT_MSG(!validParameters,
      "**** Error: Invalid JH2 material parameters: " + validationMessage);

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_jh2DamageFieldId = fieldManager.getFieldId(
      PeridigmField::ELEMENT, PeridigmField::SCALAR,
      PeridigmField::TWO_STEP, "JH2_Damage");
  m_equivalentPlasticStrainFieldId = fieldManager.getFieldId(
      PeridigmField::ELEMENT, PeridigmField::SCALAR,
      PeridigmField::TWO_STEP, "Equivalent_Plastic_Strain");
  m_jh2BulkingPressureFieldId = fieldManager.getFieldId(
      PeridigmField::ELEMENT, PeridigmField::SCALAR,
      PeridigmField::TWO_STEP, "JH2_Bulking_Pressure");
  m_jh2EquivalentStressFieldId = fieldManager.getFieldId(
      PeridigmField::ELEMENT, PeridigmField::SCALAR,
      PeridigmField::CONSTANT, "JH2_Equivalent_Stress");
  m_jh2YieldStrengthFieldId = fieldManager.getFieldId(
      PeridigmField::ELEMENT, PeridigmField::SCALAR,
      PeridigmField::CONSTANT, "JH2_Yield_Strength");
  m_jh2EOSPressureFieldId = fieldManager.getFieldId(
      PeridigmField::ELEMENT, PeridigmField::SCALAR,
      PeridigmField::CONSTANT, "JH2_EOS_Pressure");
  m_jh2EquivalentStrainRateFieldId = fieldManager.getFieldId(
      PeridigmField::ELEMENT, PeridigmField::SCALAR,
      PeridigmField::CONSTANT, "JH2_Equivalent_Strain_Rate");
  m_jh2DeformationJacobianFieldId = fieldManager.getFieldId(
      PeridigmField::ELEMENT, PeridigmField::SCALAR,
      PeridigmField::CONSTANT, "JH2_Deformation_Jacobian");
  m_jh2PlasticStrainIncrementFieldId = fieldManager.getFieldId(
      PeridigmField::ELEMENT, PeridigmField::SCALAR,
      PeridigmField::CONSTANT, "JH2_Plastic_Strain_Increment");
  m_jh2FracturePlasticStrainFieldId = fieldManager.getFieldId(
      PeridigmField::ELEMENT, PeridigmField::SCALAR,
      PeridigmField::CONSTANT, "JH2_Fracture_Plastic_Strain");
  m_jh2DissipatedEnergyFieldId = fieldManager.getFieldId(
      PeridigmField::ELEMENT, PeridigmField::SCALAR,
      PeridigmField::CONSTANT, "JH2_Dissipated_Energy");
  m_jh2TensileFailureFieldId = fieldManager.getFieldId(
      PeridigmField::ELEMENT, PeridigmField::SCALAR,
      PeridigmField::TWO_STEP, "JH2_Tensile_Failure");

  m_fieldIds.push_back(m_jh2DamageFieldId);
  m_fieldIds.push_back(m_equivalentPlasticStrainFieldId);
  m_fieldIds.push_back(m_jh2BulkingPressureFieldId);
  m_fieldIds.push_back(m_jh2EquivalentStressFieldId);
  m_fieldIds.push_back(m_jh2YieldStrengthFieldId);
  m_fieldIds.push_back(m_jh2EOSPressureFieldId);
  m_fieldIds.push_back(m_jh2EquivalentStrainRateFieldId);
  m_fieldIds.push_back(m_jh2DeformationJacobianFieldId);
  m_fieldIds.push_back(m_jh2PlasticStrainIncrementFieldId);
  m_fieldIds.push_back(m_jh2FracturePlasticStrainFieldId);
  m_fieldIds.push_back(m_jh2DissipatedEnergyFieldId);
  m_fieldIds.push_back(m_jh2TensileFailureFieldId);
}

PeridigmNS::JH2CorrespondenceMaterial::~JH2CorrespondenceMaterial() {}

void PeridigmNS::JH2CorrespondenceMaterial::initialize(
    const double dt,
    const int numOwnedPoints,
    const int* ownedIDs,
    const int* neighborhoodList,
    PeridigmNS::DataManager& dataManager) {
  PeridigmNS::CorrespondenceMaterial::initialize(
      dt, numOwnedPoints, ownedIDs, neighborhoodList, dataManager);

  dataManager.getData(m_jh2DamageFieldId,
      PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_jh2DamageFieldId,
      PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_equivalentPlasticStrainFieldId,
      PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_equivalentPlasticStrainFieldId,
      PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_jh2BulkingPressureFieldId,
      PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_jh2BulkingPressureFieldId,
      PeridigmField::STEP_NP1)->PutScalar(0.0);

  dataManager.getData(m_jh2EquivalentStressFieldId,
      PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_jh2YieldStrengthFieldId,
      PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_jh2EOSPressureFieldId,
      PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_jh2EquivalentStrainRateFieldId,
      PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_jh2DeformationJacobianFieldId,
      PeridigmField::STEP_NONE)->PutScalar(1.0);
  dataManager.getData(m_jh2PlasticStrainIncrementFieldId,
      PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_jh2FracturePlasticStrainFieldId,
      PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_jh2DissipatedEnergyFieldId,
      PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_jh2TensileFailureFieldId,
      PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_jh2TensileFailureFieldId,
      PeridigmField::STEP_NP1)->PutScalar(0.0);
}

void PeridigmNS::JH2CorrespondenceMaterial::computeCauchyStress(
    const double dt,
    const int numOwnedPoints,
    PeridigmNS::DataManager& dataManager) const {
  double *rateOfDeformation, *deformationGradient;
  double *stressN, *stressNP1;
  double *damageN, *damageNP1, *plasticStrainN, *plasticStrainNP1;
  double *bulkingPressureN, *bulkingPressureNP1;
  double *equivalentStress, *yieldStrength, *eosPressure;
  double *equivalentStrainRate, *deformationJacobian;
  double *plasticStrainIncrement, *fracturePlasticStrain;
  double *dissipatedEnergy, *tensileFailureN, *tensileFailureNP1;

  dataManager.getData(m_unrotatedRateOfDeformationFieldId,
      PeridigmField::STEP_NONE)->ExtractView(&rateOfDeformation);
  dataManager.getData(m_deformationGradientFieldId,
      PeridigmField::STEP_NONE)->ExtractView(&deformationGradient);
  dataManager.getData(m_unrotatedCauchyStressFieldId,
      PeridigmField::STEP_N)->ExtractView(&stressN);
  dataManager.getData(m_unrotatedCauchyStressFieldId,
      PeridigmField::STEP_NP1)->ExtractView(&stressNP1);
  dataManager.getData(m_jh2DamageFieldId,
      PeridigmField::STEP_N)->ExtractView(&damageN);
  dataManager.getData(m_jh2DamageFieldId,
      PeridigmField::STEP_NP1)->ExtractView(&damageNP1);
  dataManager.getData(m_equivalentPlasticStrainFieldId,
      PeridigmField::STEP_N)->ExtractView(&plasticStrainN);
  dataManager.getData(m_equivalentPlasticStrainFieldId,
      PeridigmField::STEP_NP1)->ExtractView(&plasticStrainNP1);
  dataManager.getData(m_jh2BulkingPressureFieldId,
      PeridigmField::STEP_N)->ExtractView(&bulkingPressureN);
  dataManager.getData(m_jh2BulkingPressureFieldId,
      PeridigmField::STEP_NP1)->ExtractView(&bulkingPressureNP1);

  dataManager.getData(m_jh2EquivalentStressFieldId,
      PeridigmField::STEP_NONE)->ExtractView(&equivalentStress);
  dataManager.getData(m_jh2YieldStrengthFieldId,
      PeridigmField::STEP_NONE)->ExtractView(&yieldStrength);
  dataManager.getData(m_jh2EOSPressureFieldId,
      PeridigmField::STEP_NONE)->ExtractView(&eosPressure);
  dataManager.getData(m_jh2EquivalentStrainRateFieldId,
      PeridigmField::STEP_NONE)->ExtractView(&equivalentStrainRate);
  dataManager.getData(m_jh2DeformationJacobianFieldId,
      PeridigmField::STEP_NONE)->ExtractView(&deformationJacobian);
  dataManager.getData(m_jh2PlasticStrainIncrementFieldId,
      PeridigmField::STEP_NONE)->ExtractView(&plasticStrainIncrement);
  dataManager.getData(m_jh2FracturePlasticStrainFieldId,
      PeridigmField::STEP_NONE)->ExtractView(&fracturePlasticStrain);
  dataManager.getData(m_jh2DissipatedEnergyFieldId,
      PeridigmField::STEP_NONE)->ExtractView(&dissipatedEnergy);
  dataManager.getData(m_jh2TensileFailureFieldId,
      PeridigmField::STEP_N)->ExtractView(&tensileFailureN);
  dataManager.getData(m_jh2TensileFailureFieldId,
      PeridigmField::STEP_NP1)->ExtractView(&tensileFailureNP1);

  for(int point = 0; point < numOwnedPoints; ++point) {
    JH2::State stateN;
    JH2::State stateNP1;
    JH2::Diagnostics diagnostics;
    for(int component = 0; component < 9; ++component)
      stateN.stress[component] = stressN[9 * point + component];
    stateN.damage = damageN[point];
    stateN.plasticStrain = plasticStrainN[point];
    stateN.bulkingPressure = bulkingPressureN[point];
    stateN.tensileFailure = tensileFailureN[point];

    const JH2::UpdateStatus status = JH2::updatePoint(
        m_jh2Parameters, dt, rateOfDeformation + 9 * point,
        deformationGradient + 9 * point, stateN, stateNP1, diagnostics);
    if(status != JH2::UPDATE_SUCCESS) {
      ostringstream message;
      message << "**** Error: JH2 update failed at owned point " << point
              << " with status " << static_cast<int>(status) << ".";
      TEUCHOS_TEST_FOR_TERMINATION(true, message.str());
    }

    for(int component = 0; component < 9; ++component)
      stressNP1[9 * point + component] = stateNP1.stress[component];
    damageNP1[point] = stateNP1.damage;
    plasticStrainNP1[point] = stateNP1.plasticStrain;
    bulkingPressureNP1[point] = stateNP1.bulkingPressure;
    equivalentStress[point] = diagnostics.equivalentStress;
    yieldStrength[point] = diagnostics.yieldStrength;
    eosPressure[point] = diagnostics.eosPressure;
    equivalentStrainRate[point] = diagnostics.equivalentStrainRate;
    deformationJacobian[point] = diagnostics.deformationJacobian;
    plasticStrainIncrement[point] = diagnostics.plasticStrainIncrement;
    fracturePlasticStrain[point] = diagnostics.fracturePlasticStrain;
    dissipatedEnergy[point] = diagnostics.dissipatedEnergy;
    tensileFailureNP1[point] = stateNP1.tensileFailure;
  }
}
