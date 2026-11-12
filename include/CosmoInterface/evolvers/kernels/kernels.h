#ifndef COSMOINTERFACE_EVOLVERS_KERNELS_KERNELS_H
#define COSMOINTERFACE_EVOLVERS_KERNELS_KERNELS_H

#include "CosmoInterface/evolvers/kernels/scalarsingletkernels.h"
#include "CosmoInterface/evolvers/kernels/gwskernels.h"
#include "CosmoInterface/evolvers/kernels/scalefactorkernels.h"
#include "CosmoInterface/evolvers/kernels/complexscalarkernels.h"
#include "CosmoInterface/evolvers/kernels/su2doubletkernels.h"
#include "CosmoInterface/evolvers/kernels/u1kernels.h"
#include "CosmoInterface/evolvers/kernels/su2kernels.h"

#include "CosmoInterface/fieldsnumbering.h"

#include <type_traits>

namespace TempLat::Kernels
{

  template <class Model, class KernelType, int N>
  static inline auto get(FieldsNumbering::fldS, Model &model, Tag<N> n, KernelType kt)
  {
    return ScalarSingletKernels::get_momentum(model, n, kt);
  }

  template <class Model, class KernelType, int N>
  static inline auto get(FieldsNumbering::piS, Model &model, Tag<N> n, KernelType kt)
  {
    return ScalarSingletKernels::get(model, n, kt);
  }

  template <class Model, class KernelType, int N>
  static inline auto get(FieldsNumbering::fldCS, Model &model, Tag<N> n, KernelType kt)
  {
    return ComplexScalarKernels::get_momentum(model, n, kt);
  }

  template <class Model, class KernelType, int N>
  static inline auto get(FieldsNumbering::piCS, Model &model, Tag<N> n, KernelType kt)
  {
    return ComplexScalarKernels::get(model, n, kt);
  }

  template <class Model, class KernelType, int N>
  static inline auto get(FieldsNumbering::fldSU2Doublet, Model &model, Tag<N> n, KernelType kt)
  {
    return SU2DoubletKernels::get_momentum(model, n, kt);
  }

  template <class Model, class KernelType, int N>
  static inline auto get(FieldsNumbering::piSU2Doublet, Model &model, Tag<N> n, KernelType kt)
  {
    return SU2DoubletKernels::get(model, n, kt);
  }

  template <class Model, class KernelType, int N>
  static inline auto get(FieldsNumbering::fldU1, Model &model, Tag<N> n, KernelType kt)
  {
    return U1Kernels::get_momentum(model, n, kt);
  }

  template <class Model, class KernelType, int N>
  static inline auto get(FieldsNumbering::piU1, Model &model, Tag<N> n, KernelType kt)
  {
    return U1Kernels::get(model, n, kt);
  }

  template <class Model, class KernelType, int N>
  static inline auto get(FieldsNumbering::fldSU2, Model &model, Tag<N> n, KernelType kt)
  {
    return SU2Kernels::get_momentum(model, n, kt);
  }

  template <class Model, class KernelType, int N>
  static inline auto get(FieldsNumbering::piSU2, Model &model, Tag<N> n, KernelType kt)
  {
    return SU2Kernels::get(model, n, kt);
  }
  /** @brief Kernel used by the evolvers for the field FLD, for a given kernel type. By default the kernel type selects
   * the corresponding overload of get() above.
   **/
  template <class Model, int FLD, int N, class KernelType>
  static inline auto evolutionKernel(Tag<FLD> fld, Model &model, Tag<N> n, KernelType kt)
  {
    return get(fld, model, n, kt);
  }

  /** @brief E-fold kernels: the standard equations of motion rescaled by the inverse program Hubble rate,
   * dX/dN = (a / aDot) dX/dtau. Kernels which vanish identically (ZeroType) are kept as such, so that the evolver
   * still skips the corresponding fields.
   **/
  template <class Model, int FLD, int N>
  static inline auto evolutionKernel(Tag<FLD> fld, Model &model, Tag<N> n, KernelsTypes::EFolds<Model> kt)
  {
    auto kernel = get(fld, model, n, kt.eom());
    if constexpr (std::is_same_v<std::decay_t<decltype(kernel)>, ZeroType>)
      return kernel;
    else
      return kt.invHubble * kernel;
  }
} // namespace TempLat::Kernels

#endif
