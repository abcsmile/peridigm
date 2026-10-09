//! \file jh2.h
//! \brief Material-point kernel for the JH2 constitutive and damage model.

#ifndef PERIDIGM_JH2_H
#define PERIDIGM_JH2_H

#include <string>

namespace JH2 {

//! JH2 material constants.  Stress-like quantities use Pa and time uses s.
struct Parameters {
  double bulkModulus;
  double shearModulus;
  double intactStrengthConstant;
  double intactStrengthExponent;
  double fracturedStrengthConstant;
  double fracturedStrengthExponent;
  double strainRateCoefficient;
  double damageConstantD1;
  double damageExponentD2;
  double energyConversionFactor;
  double pressureAtHEL;
  double strengthAtHEL;
  double tensilePressure;
  double referenceStrainRate;
  double eosK1;
  double eosK2;
  double eosK3;
};

//! History variables stored at one material point.
struct State {
  double stress[9];       //!< Unrotated Cauchy stress [Pa], full tensor.
  double damage;         //!< JH2 continuous damage [-], 0 <= D <= 1.
  double plasticStrain;  //!< Accumulated equivalent plastic strain [-].
  double bulkingPressure;//!< JH2 pressure history Delta p [Pa].
  double tensileFailure; //!< Persistent hydrostatic-tension failure flag, 0 or 1.
};

//! Diagnostic values produced by one constitutive update.
struct Diagnostics {
  double equivalentStress;     //!< Trial von Mises stress [Pa], Eq. (20).
  double yieldStrength;        //!< Yield strength before return [Pa], Eq. (21).
  double eosPressure;          //!< Updated EOS pressure [Pa], Eq. (25).
  double equivalentStrainRate; //!< Effective deviatoric strain rate [1/s], Eq. (24).
  double deformationJacobian;  //!< J = det(F) [-].
  double plasticStrainIncrement;//!< Effective plastic strain increment [-].
  double fracturePlasticStrain;//!< Plastic strain to fracture [-], Eq. (30).
  double dissipatedEnergy;     //!< Damage-related energy loss [J/m^3], Eqs. (27)-(28).
  bool yielded;                //!< True if s >= sigma_y.
  bool tensileFailure;         //!< True after p* <= -T* has occurred.
};

enum UpdateStatus {
  UPDATE_SUCCESS = 0,
  UPDATE_INVALID_TIME_STEP = 1,
  UPDATE_NONPOSITIVE_JACOBIAN = 2,
  UPDATE_NONFINITE_INPUT = 3,
  UPDATE_INVALID_STATE = 4
};

//! Validate material constants before entering a time-integration loop.
bool validateParameters(const Parameters& parameters, std::string& message);

//! Compute J = det(F) for a row-major full tensor [9].
double determinant(const double* deformationGradient);

//! Compute the JH2 equivalent strain rate, Eq. (24), from d [1/s].
double equivalentStrainRate(const double* rateOfDeformation);

//! Update one material point using Zhu and Zhao (2021), Eqs. (19)-(35).
//!
//! Input and output stresses are unrotated Cauchy stresses with tension
//! positive.  JH2 pressure is positive in compression.  Stress rotation is
//! deliberately excluded because Peridigm's CorrespondenceMaterial base class
//! performs that operation after this update.
UpdateStatus updatePoint(const Parameters& parameters,
                         double timeStep,
                         const double* unrotatedRateOfDeformation,
                         const double* deformationGradient,
                         const State& stateN,
                         State& stateNP1,
                         Diagnostics& diagnostics);

} // namespace JH2

#endif // PERIDIGM_JH2_H
