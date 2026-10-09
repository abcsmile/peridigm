//! \file jh2.cxx
//! \brief Material-point implementation of Zhu and Zhao (2021), Eqs. (19)-(35).

#include "jh2.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

const double tiny = 1.0e-30;

bool finiteTensor(const double* tensor) {
  for(int component = 0; component < 9; ++component)
    if(!std::isfinite(tensor[component]))
      return false;
  return true;
}

double tensorInnerProduct(const double* tensorA, const double* tensorB) {
  double result = 0.0;
  for(int component = 0; component < 9; ++component)
    result += tensorA[component] * tensorB[component];
  return result;
}

double meanStress(const double* stress) {
  return (stress[0] + stress[4] + stress[8]) / 3.0;
}

void deviator(const double* tensor, double* tensorDeviator) {
  const double mean = meanStress(tensor);
  for(int component = 0; component < 9; ++component)
    tensorDeviator[component] = tensor[component];
  tensorDeviator[0] -= mean;
  tensorDeviator[4] -= mean;
  tensorDeviator[8] -= mean;
}

double equivalentStress(const double* deviatoricStress) {
  return std::sqrt(std::max(0.0,
      1.5 * tensorInnerProduct(deviatoricStress, deviatoricStress)));
}

double rateFactor(const JH2::Parameters& parameters,
                  double equivalentStrainRate) {
  const double normalizedRate = std::max(
      equivalentStrainRate / parameters.referenceStrainRate, tiny);
  return std::max(0.0,
      1.0 + parameters.strainRateCoefficient * std::log(normalizedRate));
}

double intactStrength(const JH2::Parameters& parameters,
                      double normalizedPressure,
                      double strainRateFactor) {
  const double normalizedTension =
      parameters.tensilePressure / parameters.pressureAtHEL;
  const double pressureBase = std::max(
      normalizedPressure + normalizedTension, tiny);
  return parameters.intactStrengthConstant *
      std::pow(pressureBase, parameters.intactStrengthExponent) *
      strainRateFactor;
}

double fracturedStrength(const JH2::Parameters& parameters,
                         double normalizedPressure,
                         double strainRateFactor) {
  const double pressureBase = std::max(normalizedPressure, tiny);
  return parameters.fracturedStrengthConstant *
      std::pow(pressureBase, parameters.fracturedStrengthExponent) *
      strainRateFactor;
}

double interpolatedYieldStrength(const JH2::Parameters& parameters,
                                 double intactNormalizedStrength,
                                 double fracturedNormalizedStrength,
                                 double damage) {
  const double normalizedStrength = intactNormalizedStrength - damage *
      (intactNormalizedStrength - fracturedNormalizedStrength);
  return std::max(0.0, normalizedStrength * parameters.strengthAtHEL);
}

double eosPolynomial(const JH2::Parameters& parameters, double compression) {
  if(compression >= 0.0)
    return parameters.eosK1 * compression +
           parameters.eosK2 * compression * compression +
           parameters.eosK3 * compression * compression * compression;
  return parameters.eosK1 * compression;
}

double updateBulkingPressure(const JH2::Parameters& parameters,
                             double compression,
                             double oldBulkingPressure,
                             double damage,
                             double dissipatedEnergy) {
  // Zhu and Zhao (2021), Eq. (26).
  if(damage <= 0.0 || dissipatedEnergy <= 0.0)
    return oldBulkingPressure;

  const double linearPressure = parameters.eosK1 * compression;
  const double previousPressure = linearPressure + oldBulkingPressure;
  const double radicand = previousPressure * previousPressure +
      2.0 * parameters.energyConversionFactor * parameters.eosK1 *
      dissipatedEnergy;
  return -linearPressure + std::sqrt(std::max(0.0, radicand));
}

bool finiteState(const JH2::State& state) {
  return finiteTensor(state.stress) && std::isfinite(state.damage) &&
      std::isfinite(state.plasticStrain) &&
      std::isfinite(state.bulkingPressure) &&
      std::isfinite(state.tensileFailure);
}

} // anonymous namespace

bool JH2::validateParameters(const Parameters& p, std::string& message) {
  const double values[] = {
    p.bulkModulus, p.shearModulus, p.intactStrengthConstant,
    p.intactStrengthExponent, p.fracturedStrengthConstant,
    p.fracturedStrengthExponent, p.strainRateCoefficient,
    p.damageConstantD1, p.damageExponentD2, p.energyConversionFactor,
    p.pressureAtHEL, p.strengthAtHEL, p.tensilePressure,
    p.referenceStrainRate, p.eosK1, p.eosK2, p.eosK3
  };
  for(double value : values) {
    if(!std::isfinite(value)) {
      message = "All JH2 parameters must be finite.";
      return false;
    }
  }
  if(p.bulkModulus <= 0.0 || p.shearModulus <= 0.0 || p.eosK1 <= 0.0) {
    message = "Bulk modulus, shear modulus, and EOS K1 must be positive.";
    return false;
  }
  if(p.intactStrengthConstant < 0.0 || p.fracturedStrengthConstant < 0.0 ||
     p.intactStrengthExponent < 0.0 || p.fracturedStrengthExponent < 0.0) {
    message = "JH2 strength constants and exponents must be nonnegative.";
    return false;
  }
  if(p.strainRateCoefficient < 0.0 || p.damageConstantD1 <= 0.0 ||
     p.damageExponentD2 <= 0.0 ||
     p.pressureAtHEL <= 0.0 || p.strengthAtHEL <= 0.0 ||
     p.tensilePressure < 0.0 || p.referenceStrainRate <= 0.0) {
    message = "JH2 rate, damage, HEL, tensile, or reference-rate parameter is invalid.";
    return false;
  }
  if(p.energyConversionFactor < 0.0 || p.energyConversionFactor > 1.0) {
    message = "JH2 energy conversion factor must lie in [0, 1].";
    return false;
  }
  const double eosBulkRelativeDifference =
      std::abs(p.eosK1 - p.bulkModulus) / p.bulkModulus;
  if(eosBulkRelativeDifference > 1.0e-12) {
    message = "JH2 EOS K1 must equal the elastic bulk modulus.";
    return false;
  }
  message.clear();
  return true;
}

double JH2::determinant(const double* f) {
  return f[0] * (f[4] * f[8] - f[5] * f[7]) -
         f[1] * (f[3] * f[8] - f[5] * f[6]) +
         f[2] * (f[3] * f[7] - f[4] * f[6]);
}

double JH2::equivalentStrainRate(const double* d) {
  double dDeviator[9];
  deviator(d, dDeviator);
  return std::sqrt(std::max(0.0,
      (2.0 / 3.0) * tensorInnerProduct(dDeviator, dDeviator)));
}

JH2::UpdateStatus JH2::updatePoint(const Parameters& p,
                                    double dt,
                                    const double* d,
                                    const double* f,
                                    const State& stateN,
                                    State& stateNP1,
                                    Diagnostics& output) {
  if(!std::isfinite(dt) || dt <= 0.0)
    return UPDATE_INVALID_TIME_STEP;
  if(!finiteTensor(d) || !finiteTensor(f) || !finiteState(stateN))
    return UPDATE_NONFINITE_INPUT;
  if(stateN.damage < 0.0 || stateN.damage > 1.0 ||
     stateN.plasticStrain < 0.0 || stateN.bulkingPressure < 0.0 ||
     stateN.tensileFailure < 0.0 || stateN.tensileFailure > 1.0)
    return UPDATE_INVALID_STATE;

  const double jacobian = determinant(f);
  if(!std::isfinite(jacobian))
    return UPDATE_NONFINITE_INPUT;
  if(jacobian <= 0.0)
    return UPDATE_NONPOSITIVE_JACOBIAN;

  stateNP1 = stateN;
  output = Diagnostics{};
  output.deformationJacobian = jacobian;
  output.equivalentStrainRate = equivalentStrainRate(d);
  if(stateN.tensileFailure >= 0.5) {
    for(double& component : stateNP1.stress)
      component = 0.0;
    stateNP1.tensileFailure = 1.0;
    output.tensileFailure = true;
    return UPDATE_SUCCESS;
  }

  double strainIncrement[9];
  double deviatoricStrainIncrement[9];
  for(int component = 0; component < 9; ++component)
    strainIncrement[component] = d[component] * dt;
  deviator(strainIncrement, deviatoricStrainIncrement);
  const double volumeStrainIncrement =
      strainIncrement[0] + strainIncrement[4] + strainIncrement[8];

  // Elastic trial stress, Zhu and Zhao (2021), Eq. (19).
  double trialStress[9];
  for(int component = 0; component < 9; ++component)
    trialStress[component] = stateN.stress[component] +
        2.0 * p.shearModulus * deviatoricStrainIncrement[component];
  trialStress[0] += p.bulkModulus * volumeStrainIncrement;
  trialStress[4] += p.bulkModulus * volumeStrainIncrement;
  trialStress[8] += p.bulkModulus * volumeStrainIncrement;

  double trialDeviatoricStress[9];
  deviator(trialStress, trialDeviatoricStress);
  const double trialMeanStress = meanStress(trialStress);
  output.equivalentStress = equivalentStress(trialDeviatoricStress);

  const double compression = 1.0 / jacobian - 1.0;
  const double pressureBeforeUpdate =
      eosPolynomial(p, compression) + stateN.bulkingPressure;
  const double normalizedPressure = pressureBeforeUpdate / p.pressureAtHEL;
  const double normalizedTension = p.tensilePressure / p.pressureAtHEL;
  output.eosPressure = pressureBeforeUpdate;

  // Hydrostatic tensile failure condition preceding Eq. (22).
  if(normalizedPressure <= -normalizedTension) {
    for(double& component : stateNP1.stress)
      component = 0.0;
    stateNP1.damage = 1.0;
    stateNP1.tensileFailure = 1.0;
    output.tensileFailure = true;
    return UPDATE_SUCCESS;
  }

  const double currentRateFactor = rateFactor(p, output.equivalentStrainRate);
  const double intact = intactStrength(p, normalizedPressure, currentRateFactor);
  const double fractured = fracturedStrength(p, normalizedPressure, currentRateFactor);
  output.yieldStrength = interpolatedYieldStrength(
      p, intact, fractured, std::max(0.0, std::min(1.0, stateN.damage)));

  if(output.equivalentStress < output.yieldStrength) {
    for(int component = 0; component < 9; ++component)
      stateNP1.stress[component] = trialStress[component];
    return UPDATE_SUCCESS;
  }

  output.yielded = true;
  const double returnScale = output.yieldStrength /
      std::max(output.equivalentStress, tiny);
  double returnedDeviatoricStress[9];
  for(int component = 0; component < 9; ++component)
    returnedDeviatoricStress[component] =
        returnScale * trialDeviatoricStress[component];

  double oldDeviatoricStress[9];
  deviator(stateN.stress, oldDeviatoricStress);
  double projectedPlasticStrainIncrement[9];
  for(int component = 0; component < 9; ++component) {
    // Zhu and Zhao (2021), Eqs. (32)-(33).
    const double plasticComponent = deviatoricStrainIncrement[component] -
        (returnedDeviatoricStress[component] - oldDeviatoricStress[component]) /
        (2.0 * p.shearModulus);
    projectedPlasticStrainIncrement[component] = plasticComponent *
        trialDeviatoricStress[component] /
        std::max(output.equivalentStress, tiny);
  }
  output.plasticStrainIncrement = std::sqrt(std::max(0.0,
      (2.0 / 3.0) * tensorInnerProduct(projectedPlasticStrainIncrement,
                                       projectedPlasticStrainIncrement)));
  stateNP1.plasticStrain = stateN.plasticStrain +
      output.plasticStrainIncrement;

  // Damage accumulation, Zhu and Zhao (2021), Eqs. (29)-(30).
  const double fractureBase = std::max(
      normalizedPressure + normalizedTension, tiny);
  output.fracturePlasticStrain = p.damageConstantD1 *
      std::pow(fractureBase, p.damageExponentD2);
  output.fracturePlasticStrain = std::max(
      output.fracturePlasticStrain, tiny);
  stateNP1.damage = std::max(0.0, std::min(1.0,
      stateN.damage + output.plasticStrainIncrement /
      output.fracturePlasticStrain));

  const double damagedYieldStrength = interpolatedYieldStrength(
      p, intact, fractured, stateNP1.damage);
  output.dissipatedEnergy = std::max(0.0,
      (output.yieldStrength * output.yieldStrength -
       damagedYieldStrength * damagedYieldStrength) /
      (6.0 * p.shearModulus));
  stateNP1.bulkingPressure = updateBulkingPressure(
      p, compression, stateN.bulkingPressure, stateNP1.damage,
      output.dissipatedEnergy);
  output.eosPressure = eosPolynomial(p, compression) +
      stateNP1.bulkingPressure;

  // Peridigm tension-positive convention: Eq. (35)'s compression-positive
  // pressure increment reduces the mean Cauchy stress.
  const double finalMeanStress = trialMeanStress - stateNP1.bulkingPressure;
  for(int component = 0; component < 9; ++component)
    stateNP1.stress[component] = returnedDeviatoricStress[component];
  stateNP1.stress[0] += finalMeanStress;
  stateNP1.stress[4] += finalMeanStress;
  stateNP1.stress[8] += finalMeanStress;

  return UPDATE_SUCCESS;
}
