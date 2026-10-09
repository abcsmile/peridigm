//! \file Peridigm_JH2CorrespondenceMaterial.hpp
//! \brief Peridigm wrapper for the JH2 correspondence material.

#ifndef PERIDIGM_JH2CORRESPONDENCEMATERIAL_HPP
#define PERIDIGM_JH2CORRESPONDENCEMATERIAL_HPP

#include "Peridigm_CorrespondenceMaterial.hpp"
#include "jh2.h"

namespace PeridigmNS {

class JH2CorrespondenceMaterial : public CorrespondenceMaterial {
public:
  explicit JH2CorrespondenceMaterial(const Teuchos::ParameterList& params);
  virtual ~JH2CorrespondenceMaterial();

  virtual std::string Name() const { return "JH2 Correspondence"; }

  virtual void initialize(const double dt,
                          const int numOwnedPoints,
                          const int* ownedIDs,
                          const int* neighborhoodList,
                          PeridigmNS::DataManager& dataManager);

  virtual void computeCauchyStress(const double dt,
                                   const int numOwnedPoints,
                                   PeridigmNS::DataManager& dataManager) const;

protected:
  JH2::Parameters m_jh2Parameters;

  int m_jh2DamageFieldId;
  int m_equivalentPlasticStrainFieldId;
  int m_jh2BulkingPressureFieldId;
  int m_jh2EquivalentStressFieldId;
  int m_jh2YieldStrengthFieldId;
  int m_jh2EOSPressureFieldId;
  int m_jh2EquivalentStrainRateFieldId;
  int m_jh2DeformationJacobianFieldId;
  int m_jh2PlasticStrainIncrementFieldId;
  int m_jh2FracturePlasticStrainFieldId;
  int m_jh2DissipatedEnergyFieldId;
  int m_jh2TensileFailureFieldId;
};

} // namespace PeridigmNS

#endif // PERIDIGM_JH2CORRESPONDENCEMATERIAL_HPP
