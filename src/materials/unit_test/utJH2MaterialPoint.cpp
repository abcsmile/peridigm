//! \file utJH2MaterialPoint.cpp
//! \brief Standalone material-point verification for the JH2 kernel.

#include "jh2.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace {

JH2::Parameters graniteParameters() {
  JH2::Parameters p{};
  p.bulkModulus = 21.9e9;
  p.shearModulus = 25.7e9;
  p.intactStrengthConstant = 1.9;
  p.intactStrengthExponent = 0.62;
  p.fracturedStrengthConstant = 0.25;
  p.fracturedStrengthExponent = 0.62;
  p.strainRateCoefficient = 0.005;
  p.damageConstantD1 = 0.005;
  p.damageExponentD2 = 0.7;
  p.energyConversionFactor = 0.5;
  p.pressureAtHEL = 3.7e9;
  p.strengthAtHEL = 1.2e9;
  p.tensilePressure = 54.0e6;
  p.referenceStrainRate = 1.0;
  p.eosK1 = 21.9e9;
  p.eosK2 = -4500.0e9;
  p.eosK3 = 300000.0e9;
  return p;
}

bool close(double actual, double expected, double relativeTolerance = 1.0e-11) {
  const double scale = std::max(1.0, std::abs(expected));
  return std::abs(actual - expected) <= relativeTolerance * scale;
}

void identity(double* tensor) {
  std::fill(tensor, tensor + 9, 0.0);
  tensor[0] = tensor[4] = tensor[8] = 1.0;
}

bool testElasticResponse() {
  const JH2::Parameters p = graniteParameters();
  JH2::State stateN{};
  JH2::State stateNP1{};
  JH2::Diagnostics diagnostics{};
  double d[9] = {-500.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double f[9];
  identity(f);

  const JH2::UpdateStatus status = JH2::updatePoint(
      p, 1.0e-6, d, f, stateN, stateNP1, diagnostics);
  const double expectedSigma11 =
      (-p.bulkModulus - 4.0 * p.shearModulus / 3.0) * 5.0e-4;
  const double expectedSigma22 =
      (-p.bulkModulus + 2.0 * p.shearModulus / 3.0) * 5.0e-4;
  return status == JH2::UPDATE_SUCCESS && !diagnostics.yielded &&
      close(stateNP1.stress[0], expectedSigma11) &&
      close(stateNP1.stress[4], expectedSigma22) &&
      close(diagnostics.equivalentStrainRate, 333.3333333333333) &&
      close(stateNP1.damage, 0.0);
}

bool testZeroIncrementPreservesState() {
  const JH2::Parameters p = graniteParameters();
  JH2::State stateN{};
  stateN.stress[0] = -1.0e6;
  stateN.stress[4] = -1.0e6;
  stateN.stress[8] = -1.0e6;
  stateN.damage = 0.2;
  stateN.plasticStrain = 0.01;
  stateN.bulkingPressure = 2.0e6;
  JH2::State stateNP1{};
  JH2::Diagnostics diagnostics{};
  double d[9] = {};
  double f[9];
  identity(f);

  const JH2::UpdateStatus status = JH2::updatePoint(
      p, 1.0e-6, d, f, stateN, stateNP1, diagnostics);
  return status == JH2::UPDATE_SUCCESS && !diagnostics.yielded &&
      close(stateNP1.stress[0], stateN.stress[0]) &&
      close(stateNP1.stress[4], stateN.stress[4]) &&
      close(stateNP1.damage, stateN.damage) &&
      close(stateNP1.plasticStrain, stateN.plasticStrain) &&
      close(stateNP1.bulkingPressure, stateN.bulkingPressure);
}

bool testHydrostaticCompression() {
  const JH2::Parameters p = graniteParameters();
  JH2::State stateN{};
  JH2::State stateNP1{};
  JH2::Diagnostics diagnostics{};
  double d[9] = {-100.0, 0.0, 0.0,
                  0.0, -100.0, 0.0,
                  0.0, 0.0, -100.0};
  double f[9];
  identity(f);
  f[0] = f[4] = f[8] = 0.9999;

  const JH2::UpdateStatus status = JH2::updatePoint(
      p, 1.0e-6, d, f, stateN, stateNP1, diagnostics);
  const double expectedStress = -3.0e-4 * p.bulkModulus;
  return status == JH2::UPDATE_SUCCESS && !diagnostics.yielded &&
      close(stateNP1.stress[0], expectedStress) &&
      close(stateNP1.stress[4], expectedStress) &&
      close(stateNP1.stress[8], expectedStress) &&
      close(diagnostics.equivalentStress, 0.0) &&
      diagnostics.eosPressure > 0.0;
}

bool testPlasticDamageResponse() {
  const JH2::Parameters p = graniteParameters();
  JH2::State stateN{};
  JH2::State stateNP1{};
  JH2::Diagnostics diagnostics{};
  double d[9] = {-10000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double f[9];
  identity(f);

  const JH2::UpdateStatus status = JH2::updatePoint(
      p, 1.0e-6, d, f, stateN, stateNP1, diagnostics);
  double stressDeviator[9];
  const double mean =
      (stateNP1.stress[0] + stateNP1.stress[4] + stateNP1.stress[8]) / 3.0;
  for(int i = 0; i < 9; ++i)
    stressDeviator[i] = stateNP1.stress[i];
  stressDeviator[0] -= mean;
  stressDeviator[4] -= mean;
  stressDeviator[8] -= mean;
  double inner = 0.0;
  for(double value : stressDeviator)
    inner += value * value;
  const double returnedEquivalentStress = std::sqrt(1.5 * inner);

  return status == JH2::UPDATE_SUCCESS && diagnostics.yielded &&
      close(diagnostics.equivalentStress, 5.14e8) &&
      close(diagnostics.yieldStrength, 1.731587922556829e8) &&
      close(diagnostics.plasticStrainIncrement, 2.5523315569662593e-3) &&
      close(diagnostics.fracturePlasticStrain, 2.5936085510646473e-4) &&
      close(stateNP1.damage, 1.0) &&
      close(stateNP1.plasticStrain, 2.5523315569662593e-3) &&
      close(diagnostics.dissipatedEnergy, 1.9444855600160017e5) &&
      close(stateNP1.bulkingPressure, 6.525659642086035e7) &&
      close(returnedEquivalentStress, diagnostics.yieldStrength, 1.0e-10);
}

bool testDissipatedEnergyDefinitions() {
  const JH2::Parameters p = graniteParameters();
  JH2::State stateN{};
  JH2::State stateNP1{};
  JH2::Diagnostics diagnostics{};
  double d[9] = {-10000.0, 0.0, 0.0,
                   0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0};
  double f[9];
  identity(f);

  const JH2::UpdateStatus status = JH2::updatePoint(
      p, 1.0e-6, d, f, stateN, stateNP1, diagnostics);

  // Zhu and Zhao (2021), Eqs. (27)-(28): elastic energy released as
  // damage reduces the yield surface.  This point reaches D = 1 and its
  // fractured strength is numerically zero because the pressure is zero.
  const double energyDifference =
      diagnostics.yieldStrength * diagnostics.yieldStrength /
      (6.0 * p.shearModulus);

  // Alternative research diagnostic used in the companion NOSB study:
  // trapezoidal plastic work over the radial-return increment.  It is not
  // used by the production JH2 kernel because the two definitions are not
  // equivalent and the cited equations do not select this replacement.
  const double trapezoidalPlasticWork = 0.5 *
      (diagnostics.equivalentStress + diagnostics.yieldStrength) *
      diagnostics.plasticStrainIncrement;

  return status == JH2::UPDATE_SUCCESS && diagnostics.yielded &&
      close(diagnostics.dissipatedEnergy, energyDifference) &&
      close(energyDifference, 1.9444855600160017e5) &&
      close(trapezoidalPlasticWork, 8.7692853506050073e5) &&
      trapezoidalPlasticWork > 4.0 * energyDifference;
}

bool testPureShearRadialReturn() {
  const JH2::Parameters p = graniteParameters();
  JH2::State stateN{};
  JH2::State stateNP1{};
  JH2::Diagnostics diagnostics{};
  double d[9] = {};
  d[1] = d[3] = 10000.0;
  double f[9];
  identity(f);

  const JH2::UpdateStatus status = JH2::updatePoint(
      p, 1.0e-6, d, f, stateN, stateNP1, diagnostics);
  double stressDeviator[9];
  const double mean =
      (stateNP1.stress[0] + stateNP1.stress[4] + stateNP1.stress[8]) / 3.0;
  double inner = 0.0;
  for(int component = 0; component < 9; ++component) {
    stressDeviator[component] = stateNP1.stress[component];
  }
  stressDeviator[0] -= mean;
  stressDeviator[4] -= mean;
  stressDeviator[8] -= mean;
  for(double value : stressDeviator)
    inner += value * value;
  const double returnedEquivalentStress = std::sqrt(1.5 * inner);
  return status == JH2::UPDATE_SUCCESS && diagnostics.yielded &&
      close(returnedEquivalentStress, diagnostics.yieldStrength, 1.0e-10);
}

bool testRateDependence() {
  const JH2::Parameters p = graniteParameters();
  JH2::State stateN{};
  JH2::State lowRateState{};
  JH2::State highRateState{};
  JH2::Diagnostics lowRateDiagnostics{};
  JH2::Diagnostics highRateDiagnostics{};
  double lowRate[9] = {-10000.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0};
  double highRate[9] = {-100000.0, 0.0, 0.0,
                         0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0};
  double f[9];
  identity(f);

  const JH2::UpdateStatus lowStatus = JH2::updatePoint(
      p, 1.0e-6, lowRate, f, stateN, lowRateState, lowRateDiagnostics);
  const JH2::UpdateStatus highStatus = JH2::updatePoint(
      p, 1.0e-7, highRate, f, stateN, highRateState, highRateDiagnostics);
  return lowStatus == JH2::UPDATE_SUCCESS &&
      highStatus == JH2::UPDATE_SUCCESS &&
      highRateDiagnostics.equivalentStrainRate >
          lowRateDiagnostics.equivalentStrainRate &&
      highRateDiagnostics.yieldStrength > lowRateDiagnostics.yieldStrength;
}

bool testMonotonicHistory() {
  JH2::Parameters p = graniteParameters();
  p.damageConstantD1 = 0.5;
  JH2::State state0{};
  JH2::State state1{};
  JH2::State state2{};
  JH2::Diagnostics diagnostics1{};
  JH2::Diagnostics diagnostics2{};
  double d[9] = {-5000.0, 0.0, 0.0,
                  0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0};
  double f1[9];
  double f2[9];
  identity(f1);
  identity(f2);
  f1[0] = 0.995;
  f2[0] = 0.990;

  const JH2::UpdateStatus firstStatus = JH2::updatePoint(
      p, 1.0e-6, d, f1, state0, state1, diagnostics1);
  const JH2::UpdateStatus secondStatus = JH2::updatePoint(
      p, 1.0e-6, d, f2, state1, state2, diagnostics2);
  return firstStatus == JH2::UPDATE_SUCCESS &&
      secondStatus == JH2::UPDATE_SUCCESS &&
      state1.damage > 0.0 && state2.damage >= state1.damage &&
      state2.plasticStrain >= state1.plasticStrain &&
      state2.bulkingPressure >= state1.bulkingPressure;
}

bool testHydrostaticTensileFailure() {
  const JH2::Parameters p = graniteParameters();
  JH2::State stateN{};
  JH2::State stateNP1{};
  JH2::Diagnostics diagnostics{};
  double d[9] = {};
  double f[9];
  identity(f);
  f[0] = f[4] = f[8] = 1.02;

  const JH2::UpdateStatus status = JH2::updatePoint(
      p, 1.0e-6, d, f, stateN, stateNP1, diagnostics);
  double stressNorm = 0.0;
  for(double value : stateNP1.stress)
    stressNorm += value * value;
  if(status != JH2::UPDATE_SUCCESS || !diagnostics.tensileFailure ||
     !close(stateNP1.damage, 1.0) || !close(stateNP1.tensileFailure, 1.0) ||
     !close(stressNorm, 0.0))
    return false;

  // A failed point must not recover stress if compression later resumes.
  JH2::State stateAfterFailure{};
  JH2::Diagnostics diagnosticsAfterFailure{};
  double compressionRate[9] = {-10000.0, 0.0, 0.0,
                                0.0, 0.0, 0.0,
                                0.0, 0.0, 0.0};
  identity(f);
  const JH2::UpdateStatus secondStatus = JH2::updatePoint(
      p, 1.0e-6, compressionRate, f, stateNP1, stateAfterFailure,
      diagnosticsAfterFailure);
  double recoveredStressNorm = 0.0;
  for(double value : stateAfterFailure.stress)
    recoveredStressNorm += value * value;
  return secondStatus == JH2::UPDATE_SUCCESS &&
      diagnosticsAfterFailure.tensileFailure &&
      close(stateAfterFailure.tensileFailure, 1.0) &&
      close(recoveredStressNorm, 0.0);
}

bool testInvalidJacobian() {
  const JH2::Parameters p = graniteParameters();
  JH2::State stateN{};
  JH2::State stateNP1{};
  JH2::Diagnostics diagnostics{};
  double d[9] = {};
  double f[9] = {};
  return JH2::updatePoint(p, 1.0e-6, d, f, stateN, stateNP1,
                          diagnostics) == JH2::UPDATE_NONPOSITIVE_JACOBIAN;
}

bool testInvalidHistoryState() {
  const JH2::Parameters p = graniteParameters();
  JH2::State stateN{};
  stateN.damage = 1.1;
  JH2::State stateNP1{};
  JH2::Diagnostics diagnostics{};
  double d[9] = {};
  double f[9];
  identity(f);
  return JH2::updatePoint(p, 1.0e-6, d, f, stateN, stateNP1,
                          diagnostics) == JH2::UPDATE_INVALID_STATE;
}

bool testShearStrainRateContraction() {
  double d[9] = {};
  d[1] = 10.0;
  d[3] = 10.0;
  return close(JH2::equivalentStrainRate(d), std::sqrt(400.0 / 3.0));
}

bool testParameterValidation() {
  JH2::Parameters p = graniteParameters();
  std::string message;
  p.eosK1 = 0.99 * p.bulkModulus;
  const bool rejectsInconsistentK1 = !JH2::validateParameters(p, message);
  p = graniteParameters();
  p.energyConversionFactor = 1.1;
  const bool rejectsInvalidBeta = !JH2::validateParameters(p, message);
  return rejectsInconsistentK1 && rejectsInvalidBeta;
}

} // anonymous namespace

int main() {
  std::string validationMessage;
  if(!JH2::validateParameters(graniteParameters(), validationMessage)) {
    std::cerr << "Parameter validation failed: " << validationMessage << '\n';
    return 1;
  }

  const bool elastic = testElasticResponse();
  const bool zeroIncrement = testZeroIncrementPreservesState();
  const bool hydrostaticCompression = testHydrostaticCompression();
  const bool plastic = testPlasticDamageResponse();
  const bool energyDefinitions = testDissipatedEnergyDefinitions();
  const bool pureShear = testPureShearRadialReturn();
  const bool rateDependence = testRateDependence();
  const bool monotonicHistory = testMonotonicHistory();
  const bool tension = testHydrostaticTensileFailure();
  const bool jacobian = testInvalidJacobian();
  const bool history = testInvalidHistoryState();
  const bool shearRate = testShearStrainRateContraction();
  const bool parameters = testParameterValidation();

  std::cout << "JH2 elastic response: " << (elastic ? "PASS" : "FAIL") << '\n'
            << "JH2 zero increment: " << (zeroIncrement ? "PASS" : "FAIL") << '\n'
            << "JH2 hydrostatic compression: "
            << (hydrostaticCompression ? "PASS" : "FAIL") << '\n'
            << "JH2 plastic damage response: " << (plastic ? "PASS" : "FAIL") << '\n'
            << "JH2 dissipated-energy definitions: "
            << (energyDefinitions ? "PASS" : "FAIL") << '\n'
            << "JH2 pure-shear return: " << (pureShear ? "PASS" : "FAIL") << '\n'
            << "JH2 rate dependence: " << (rateDependence ? "PASS" : "FAIL") << '\n'
            << "JH2 monotonic history: " << (monotonicHistory ? "PASS" : "FAIL") << '\n'
            << "JH2 tensile failure: " << (tension ? "PASS" : "FAIL") << '\n'
            << "JH2 invalid Jacobian rejection: " << (jacobian ? "PASS" : "FAIL") << '\n'
            << "JH2 invalid history rejection: " << (history ? "PASS" : "FAIL") << '\n'
            << "JH2 shear-rate contraction: " << (shearRate ? "PASS" : "FAIL") << '\n'
            << "JH2 parameter validation: " << (parameters ? "PASS" : "FAIL") << '\n';
  return elastic && zeroIncrement && hydrostaticCompression && plastic &&
      energyDefinitions &&
      pureShear && rateDependence && monotonicHistory && tension && jacobian &&
      history && shearRate && parameters ? 0 : 1;
}
