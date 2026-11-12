#ifndef AXIONU1CS_M2PHI2_H // Usual macro guard to prevent multiple inclusion
#define AXIONU1CS_M2PHI2_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Nicolas Loayza, Year: 2026

#include "CosmoInterface/cosmointerface.h"

namespace TempLat
{
  /////////
  // Model: axion phi with quadratic potential, coupled to a U(1) gauge field A through (alpha/4 Lambda) phi F Fdual,
  // plus a complex scalar chi charged under the SAME U(1), with a mass term m_chi^2 |chi|^2.
  //
  //   V(phi, chi) = 1/2 m^2 phi^2 + m_chi^2 |chi|^2
  //
  // Program variables: fStar = phi_0 (initial axion amplitude), omegaStar = m (axion mass), alpha = 0.
  //   V / (fStar^2 omegaStar^2) = 1/2 phi~^2 + (m_chi/m)^2 |chi~|^2
  /////////

  struct ModelPars : public TempLat::DefaultModelPars {
    static constexpr size_t NScalars = 1;  // the axion
    static constexpr size_t NCScalars = 1; // the charged complex scalar
    static constexpr size_t NU1Flds = 1;   // the U(1) gauge field
    static constexpr size_t NPotTerms = 2;

    using FloatType = double;

    // The complex scalar (index 0) is charged under the U(1) field (index 0)...
    typedef TempLat::CouplingsManager<NCScalars, NU1Flds, true> CsU1Couplings;
    // ... and the axion (index 0) couples to the same U(1) field through phi F Fdual.
    typedef TempLat::CouplingsManager<NScalars, NU1Flds, true> ScalarU1AxionCouplings;
  };

#define MODELNAME axionU1CS_m2phi2

  template <class R> using Model = MakeModel(R, ModelPars);

  class MODELNAME : public Model<MODELNAME>
  {
  private:
    FloatType mass;   // axion mass (GeV)
    FloatType massCS; // complex scalar mass (GeV)
    FloatType muCS2;  // (massCS / omegaStar)^2: complex scalar mass squared in program units

  public:
    static constexpr size_t NDim = Model<MODELNAME>::NDim;

    MODELNAME(ParameterParser &parser, RunParameters<FloatType> &runPar,
              device::memory::host_ptr<MemoryToolBox<NDim>> toolBox)
        : Model<MODELNAME>(parser, runPar.getLatParams(), toolBox, runPar.dt, STRINGIFY(MODELLABEL))
    {
      /////////
      // Initial homogeneous components of the fields
      /////////

      // Axion (GeV and GeV^2)
      fldS0 = parser.get<FloatType, 1>("initial_amplitudes");
      piS0 = parser.get<FloatType, 1>("initial_momenta");

      // Complex scalar norm |chi| and its derivative (GeV and GeV^2), split equally between the two components as
      // in lphi4U1. Default: no homogeneous component, only vacuum fluctuations.
      FloatType normCmplx0 = parser.get<FloatType>("cmplx_field_initial_norm", 0);
      FloatType normPiCmplx0 = parser.get<FloatType>("cmplx_momentum_initial_norm", 0);
      fldCS0(0_c) = Complexify(normCmplx0 / sqrt(FloatType(2)), normCmplx0 / sqrt(FloatType(2)));
      piCS0(0_c) = Complexify(normPiCmplx0 / sqrt(FloatType(2)), normPiCmplx0 / sqrt(FloatType(2)));

      /////////
      // Parameters of the model
      /////////
      // --> The U(1) gauge coupling and charge of chi are set in the parameter file (gU1s, CSU1_charges), the
      //     axion-gauge coupling through alphaLambda_AxionU1.
      mass = parser.get<FloatType>("Mass");
      massCS = parser.get<FloatType>("MassCS");

      /////////
      // Rescaling for program variables
      /////////
      alpha = 0;
      fStar = fldS0(0_c);
      omegaStar = mass;
      muCS2 = pow<2>(massCS / omegaStar);

      setInitialPotentialAndMassesFromPotential();

      // The complex scalar mass is set by hand:
      //  - with the convention used for complex scalars (L = |d chi|^2 - V, chi = c0 + i c1, EoM c0'' = ... - 1/2
      //    dV/dc0) the mass of each component is m_chi, i.e. mu^2 = (m_chi/omegaStar)^2 in program units, which is
      //    what the vacuum fluctuations (omega_k^2 = k^2 + masses2CS) must use;
      //  - the generic computation from the second derivative of the potential with respect to the norm would give
      //    2 mu^2, and is ill-defined (0/0) when the homogeneous component of chi vanishes.
      masses2CS(0_c) = Complexify(muCS2, muCS2);
    }

    /////////
    // Program potential
    /////////

    auto potentialTerms(Tag<0>) // Term 0: axion mass term
    {
      return FloatType(0.5) * pow<2>(fldS(0_c));
    }

    auto potentialTerms(Tag<1>) // Term 1: complex scalar mass term, (m_chi/m)^2 |chi|^2
    {
      return muCS2 * pow<2>(norm(fldCS(0_c)));
    }

    /////////
    // First derivatives
    /////////

    auto potDeriv(Tag<0>) // with respect to the axion
    {
      return fldS(0_c);
    }

    auto potDerivNormCS(Tag<0>) // with respect to the norm of the complex scalar
    {
      return 2 * muCS2 * norm(fldCS(0_c));
    }

    /////////
    // Second derivatives
    /////////

    auto potDeriv2(Tag<0>) // with respect to the axion
    {
      return OneType();
    }

    auto potDeriv2NormCS(Tag<0>) // with respect to the norm of the complex scalar
    {
      return 2 * muCS2;
    }
  };
} // namespace TempLat

#endif // AXIONU1CS_M2PHI2_H
