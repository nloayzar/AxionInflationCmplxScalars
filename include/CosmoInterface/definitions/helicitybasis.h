#ifndef COSMOINTERFACE_DEFINITIONS_HELICITYBASIS_H
#define COSMOINTERFACE_DEFINITIONS_HELICITYBASIS_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Nicolas Loayza, Year: 2026

#include <complex>

#include "TempLat/util/constants.h"
#include "TempLat/lattice/algebra/coordinates/dimensioncountrecorder.h"
#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"

namespace TempLat
{

  /** @brief Fourier-space functor returning one Cartesian component of the circular polarisation vectors
   *
   *     eps_{+/-}(k) = ( e_theta +/- i e_phi ) / sqrt(2)
   *                  = 1/sqrt(2) ( cos(th)cos(ph) -/+ i sin(ph),  cos(th)sin(ph) +/- i cos(ph),  -sin(th) ),
   *
   * where (th, ph) are the polar angles of the lattice momentum ktilde_i = sin(pi n_i / N). This is the direction
   * orthogonal to which a vector field must be (after the half-link phase shift) to satisfy the backward-derivative
   * Gauss constraint exactly. eps_+ is the positive-helicity state: i k x eps_+ = +|k| eps_+.
   * On the k_z axis (where ph is undefined) ph = 0 is used.
   *
   * The vector is always evaluated with the site's own coordinates, also on the redundant (Hermitian partner) sites of
   * the R2C layout. Away from Nyquist components eps(-k) = eps(k)^*, so a field eps(k) xi(k) with Hermitian xi is
   * Hermitian. On the Nyquist planes sin(pi n/N) is not odd under n -> -n; there the c2r transform keeps the Hermitian
   * part of the two partners, each of which is transverse to its own lattice momentum, so that the backward-derivative
   * Gauss constraint still holds exactly (the helicity of these few modes is then only approximate). Evaluating on the
   * partner and conjugating instead would violate the Gauss constraint on those planes.
   *
   * Only implemented for NDim = 3.
   */
  template <typename T, size_t NDim, bool positive> class HelicityVectorHelper : public DimensionCountRecorder<NDim>
  {
    static_assert(NDim == 3, "HelicityVector is only defined for NDim = 3.");

  public:
    using ToolboxPtr = device::memory::host_ptr<MemoryToolBox<NDim>>;

    /** @param component Cartesian component, 1, 2 or 3. */
    HelicityVectorHelper(ToolboxPtr pToolBox, int component)
        : DimensionCountRecorder<NDim>(SpaceStateType::Fourier), mLayout(pToolBox->mLayouts.getFourierSpaceLayout()),
          mN(pToolBox->mNGridPointsVec[0]), mComp(component)
    {
    }

    template <typename... IDX>
      requires IsVariadicNDIndex<NDim, IDX...>
    DEVICE_FORCEINLINE_FUNCTION complex<T> eval(const IDX &...idx) const
    {
      device::IdxArray<NDim> coord;
      mLayout.putSpatialLocationFromMemoryIndexInto(coord, idx...);

      // Direction of the lattice momentum ktilde_i = sin(pi n_i / N).
      T s[NDim];
      for (size_t d = 0; d < NDim; ++d)
        s[d] = sin(Constants::pi<T> * T(coord[d]) / T(mN));

      const T rho2 = s[0] * s[0] + s[1] * s[1];
      const T kn = sqrt(rho2 + s[2] * s[2]);
      if (kn == T(0)) return complex<T>(T(0), T(0));

      const T rho = sqrt(rho2);
      const T cosTh = s[2] / kn;
      const T sinTh = rho / kn;
      const T cosPh = rho > T(0) ? s[0] / rho : T(1);
      const T sinPh = rho > T(0) ? s[1] / rho : T(0);

      const T sgn = positive ? T(1) : T(-1);
      const T invSqrt2 = T(1) / sqrt(T(2));

      T re, im;
      if (mComp == 1) {
        re = cosTh * cosPh;
        im = -sgn * sinPh;
      } else if (mComp == 2) {
        re = cosTh * sinPh;
        im = sgn * cosPh;
      } else {
        re = -sinTh;
        im = T(0);
      }

      return complex<T>(invSqrt2 * re, invSqrt2 * im);
    }

  private:
    LayoutStruct<NDim> mLayout;
    device::Idx mN;
    int mComp;
  };

  template <typename T, size_t NDim> using HelicityVectorPlus = HelicityVectorHelper<T, NDim, true>;
  template <typename T, size_t NDim> using HelicityVectorMinus = HelicityVectorHelper<T, NDim, false>;

} // namespace TempLat

#endif
