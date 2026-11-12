#ifndef COSMOINTERFACE_EVOLVERS_KERNELS_KERNELSTYPES_H
#define COSMOINTERFACE_EVOLVERS_KERNELS_KERNELSTYPES_H
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio  Year: 2024

#include "CosmoInterface/definitions/gausslaws.h"
#include "TempLat/lattice/algebra/su2algebra/su2algebra.h"
#include "TempLat/lattice/field/collections/vectorfieldcollection.h"

namespace TempLat::KernelsTypes
{
  /** @brief A mechanism to define different kernels for the same solver.
   *
   **/
  template <typename Model> class EoM
  {
    using T = typename Model::FloatType;

  public:
    T tMinust0 = T{};

    void cache(Model &model, T tIn) { tMinust0 = tIn; }

  private:
    void cache(Model &model) {}
  };

  template <typename Model> class Diffusion
  {
    using T = typename Model::FloatType;

  public:
    T tMinust0 = T{};

    void cache(Model &model, T tIn) { tMinust0 = tIn; }

  private:
    void cache(Model &model) {}
  };

  /** @brief Equations of motion with the number of e-folds N = log(a) as time variable.
   *
   * With the program time tau (dt_cosmic = a^alpha dtau / omegaStar) the e-folds advance as dN = (aDot / a) dtau, with
   * aDot = da/dtau. Every field equation dX/dtau = K_EoM(X) therefore becomes dX/dN = (a / aDot) K_EoM(X): the e-fold
   * kernels are the standard EoM kernels rescaled by the inverse program Hubble rate a/aDot, which is cached at every
   * Runge-Kutta stage. The scale factor is not evolved, a = a(N0) exp(N - N0), and the program Hubble rate
   * aDot / a is evolved instead (see HubbleRateKernels). Within the kernels tMinust0 is the number of e-folds elapsed
   * since the initial time (so that e.g. tNonLinearAxionU1 is given in e-folds).
   **/
  template <typename Model> class EFolds
  {
    using T = typename Model::FloatType;

  public:
    T tMinust0 = T{};
    T invHubble = T{1}; // a / aDot = dtau / dN

    void cache(Model &model, T tIn)
    {
      tMinust0 = tIn;
      invHubble = model.aI / model.aDotI;
    }

    // The EoM kernel type at the same time, used to build the e-fold kernels from the standard ones.
    EoM<Model> eom() const
    {
      EoM<Model> res;
      res.tMinust0 = tMinust0;
      return res;
    }
  };

  template <class Model> auto getDt(Model &model, EoM<Model> &kt) { return model.dt; }
  template <class Model> auto setDt(Model &model, EoM<Model> &kt, typename Model::FloatType dt) { model.dt = dt; }

  template <class Model> auto getDt(Model &model, Diffusion<Model> &kt) { return model.dt; }
  template <class Model> auto getDt(Model &model, EFolds<Model> &kt) { return model.dt; }
  template <class Model> auto setDt(Model &model, EFolds<Model> &kt, typename Model::FloatType dt) { model.dt = dt; }
  template <class Model> auto setDt(Model &model, Diffusion<Model> &kt, typename Model::FloatType dt) { model.dt = dt; }
} // namespace TempLat::KernelsTypes

#endif
