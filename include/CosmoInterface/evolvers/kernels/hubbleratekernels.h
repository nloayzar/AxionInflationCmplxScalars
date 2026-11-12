#ifndef COSMOINTERFACE_EVOLVERS_KERNELS_HUBBLERATEKERNELS_H
#define COSMOINTERFACE_EVOLVERS_KERNELS_HUBBLERATEKERNELS_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Nicolas Loayza, Year: 2026

#include "CosmoInterface/definitions/energies.h"
#include "CosmoInterface/evolvers/kernels/kernelstypes.h"

namespace TempLat
{

  /** @brief Kernel of the program Hubble rate when the number of e-folds N = log(a) is the time variable.
   *
   * The evolved background variable is the program Hubble rate HP = aDot / a (aDot = da/dtau, with the program time
   * dtau = omegaStar a^-alpha dt_cosmic), i.e. HP = a^alpha H / omegaStar. Using the second Friedmann equation in the
   * form dH/dt = -(rho + p) / (2 MPl^2), which does not involve the potential,
   *
   *   dHP/dN = alpha HP - a^(2 alpha) (fStar/MPl)^2 (3 K + G + 2 (E + B)) / (3 HP),
   *
   * with K, G and E + B the volume averaged kinetic, gradient and gauge (electric + magnetic) energy densities in
   * program units (rho + p = 2 K + 2/3 G + 4/3 (E + B)). For alpha = 0 this is the kernel of the old e-fold
   * implementation, dH/dN = -(fStar/MPl)^2 (3 K + G + 2 (E + B)) / (3 H).
   *
   * As in ScaleFactorKernels, for axion-U(1) models the energy of the U(1) fields only sources the expansion in the
   * non-linear regime, tMinust0 > tNonLinearAxionU1 (in e-folds).
   **/
  class HubbleRateKernels
  {
  public:
    HubbleRateKernels() = delete;

    template <class Model> static auto get(Model &model, KernelsTypes::EFolds<Model> kt)
    {
      static_assert(not Model::IsNonMinimallyCoupled,
                    "The e-fold evolution is not implemented for non-minimally coupled models.");
      using T = decltype(model.aI);

      T kinetic = 0;  // 3 K in rho + p
      T gradient = 0; // G
      T gauge = 0;    // E + B

      if constexpr (Model::Ns > 0) {
        kinetic += Energies::kineticS(model);
        gradient += Energies::gradientS(model);
      }
      if constexpr (Model::NCs > 0) {
        kinetic += Energies::kineticCS(model);
        gradient += Energies::gradientCS(model);
      }
      if constexpr (Model::NSU2Doublet > 0) {
        kinetic += Energies::kineticSU2Doublet(model);
        gradient += Energies::gradientSU2Doublet(model);
      }
      if constexpr (Model::NU1 > 0) {
        T nonLinearSwitch = 1;
        if constexpr (Model::ScalarU1AxionCouplings::active)
          nonLinearSwitch = (kt.tMinust0 > model.tNonLinearAxionU1) ? 1 : 0;
        gauge += nonLinearSwitch * (Energies::electricU1(model) + Energies::magneticU1(model));
      }
      if constexpr (Model::NSU2 > 0) gauge += Energies::electricSU2(model) + Energies::magneticSU2(model);

      const T hubble = model.aDotI / model.aI;
      return model.alpha * hubble - pow(model.aI, 2 * model.alpha) * pow<2>(model.fStar / Model::MPl) *
                                        (3 * kinetic + gradient + 2 * gauge) / (3 * hubble);
    }
  };

} // namespace TempLat

#endif
