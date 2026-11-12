#ifndef COSMOINTERFACE_INITIALIZERS_EXTERNALPOWERSPECTRUMINITIALIZER_H
#define COSMOINTERFACE_INITIALIZERS_EXTERNALPOWERSPECTRUMINITIALIZER_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Daniel G. Figueroa, Adrien Florio, Francisco Torrenti, Ander Urio and Nicolas Loayza
// Year: 2025

#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/random/randomgaussianfield.h"
#include "TempLat/lattice/algebra/coordinates/wavenumber.h"
#include "TempLat/lattice/algebra/operators/operators.h"
#include "TempLat/lattice/algebra/helpers/getngrid.h"
#include "TempLat/lattice/algebra/complexalgebra/asfourier.h"
#include "TempLat/lattice/algebra/constants/symbols.h"
#include "TempLat/lattice/algebra/adapter/momentuminterpolator.h"

#include "TempLat/lattice/algebra/coordinates/momentummultiplicity.h"
#include "CosmoInterface/definitions/phaseBunchDavies.h"
#include "CosmoInterface/definitions/helicitybasis.h"

#include <fstream>
#include <sstream>

namespace TempLat
{

  MakeException(PSTypeINotSupportedForNDIMDifferentFrom3);
  MakeException(ExternalGaugePSFileError);

  /** @brief  A class that contains functions related with initial random fluctuations for scalar fields given by an
   * external power spectrum.
   *
   *  The functions of this class deal with the initial random fluctuations of the scalar singlets.
   *  As described in the documentation, they consist in a sum of left-moving and right-moving waves in momentum space,
   *  whose amplitudes follow a Gaussian random distribution with variance given by and external power spectrum
   * function.
   */
  template <typename T> class ExternalPowerSpectrumInitializer
  {
  public:
    /* Put public methods here. These should change very little over time. */
    template <class Model>
    ExternalPowerSpectrumInitializer(Model &model, T pLSide, std::string pSeed) : lSide(pLSide), baseSeed(pSeed)
    {
    }

    template <class Model>
    void readSpectrumFile(Model &model, const std::string &str, std::vector<T> &kIN, std::vector<T> &PSamp,
                          std::vector<T> &PSampConj) const
    {
      int lineCounter = 0;
      std::ifstream t;
      t.open(str);
      std::string tmp;

      // This just to count number of lines of the input file
      if (t.good()) {
        while (getline(t, tmp))
          lineCounter++;
      }
      t.close();

      kIN = std::vector<T>(lineCounter);
      PSamp = std::vector<T>(lineCounter);
      PSampConj = std::vector<T>(lineCounter);

      t.open(str);
      if (t.good()) {
        for (int j = 0; j < lineCounter; j++) {
          t >> kIN[j];
          kIN[j] /= model.omegaStar; // Now is in omegaStar units, i.e. it's kappa = k/omegaStar
          t >> PSamp[j];
          PSamp[j] *= model.omegaStar; // Now is in 1/omegaStar units, i.e. \tilde{P} = omegaStar * P(k), e.g. P(k) =
                                       // 1/(2\omega_k)
          t >> PSampConj[j];
          PSampConj[j] /= model.omegaStar; // Now is in 1/omegaStar units, i.e. \tilde{P} = omegaStar * P(k), e.g. P(k)
                                           // = 1/(2\omega_k)
        }
      }
      t.close();
    }

    template <class Model, size_t NDim>
    auto getInputFluctuationsNormTypeI(const MomentumInterpolator<T, NDim> &PSinterp, Model &model, Field<T, NDim> f,
                                       T kCutOff) const
    {
      if (NDim != 3)
        throw(PSTypeINotSupportedForNDIMDifferentFrom3("Initialization with PSType = 1 is only supported for NDim=3."));

      FourierSite<NDim> ntilde(f.getToolBox());
      auto k = ntilde.norm() * f.getKIR();
      auto Hcut = heaviside(kCutOff - k);

      auto invmult = MomentumMultiplicity<T, NDim>(f.getToolBox());

      // @label:externalps_typeI_norm
      return Hcut * (model.omegaStar / model.fStar) * (T(1) / pow(model.aI, T(1.5))) *
             pow(lSide / pow<2>(f.getDx()), T(1.5)) *
             safeSqrt(T(4) * Constants::pi<T> * pow(ntilde.norm(), T(2)) * invmult * PSinterp / T(2));
      // @endlabel
    }

    template <class Model, size_t NDim>
    auto getInputFluctuationsNormTypeII(const MomentumInterpolator<T, NDim> &PSinterp, Model &model, Field<T, NDim> f,
                                        T kCutOff) const
    {
      FourierSite<NDim> ntilde(f.getToolBox());
      auto k = ntilde.norm() * f.getKIR();
      auto Hcut = heaviside(kCutOff - k);

      // @label:externalps_typeII_norm
      return Hcut * (model.omegaStar / model.fStar) * (T(1) / pow(model.aI, T(1.5))) *
             pow(lSide / pow<2>(f.getDx()), T(1.5)) * safeSqrt(PSinterp) / sqrt(T(2));
      // @endlabel
    }

    template <class Model, size_t NDim>
    auto getInputFluctuationsTypeI(const MomentumInterpolator<T, NDim> &PSinterp, Model &model, Field<T, NDim> f,
                                   const std::string &mySeed, T kCutOff) const
    {
      auto norm = getInputFluctuationsNormTypeI<Model, NDim>(PSinterp, model, f, kCutOff);
      return norm * RandomGaussianField<T, NDim>(baseSeed + mySeed + f.toString(), f.getToolBox());
    }

    template <class Model, size_t NDim>
    auto getInputFluctuationsTypeII(const MomentumInterpolator<T, NDim> &PSinterp, Model &model, Field<T, NDim> f,
                                    const std::string &mySeed, T kCutOff) const
    {
      // @label:externalps_random_gaussian_field
      auto norm = getInputFluctuationsNormTypeII<Model, NDim>(PSinterp, model, f, kCutOff);
      return norm * RandomGaussianField<T, NDim>(baseSeed + mySeed + f.toString(), f.getToolBox());
      // @endlabel
    }

    template <class Model, size_t NDim>
    void conjugateGaussianInputFluctuations(Model &model, Field<T, NDim> f, Field<T, NDim> p, const std::string &str,
                                            T kCutOff, int PSType) const
    {
      std::vector<T> kIN, PSfld, PSmom;
      readSpectrumFile(model, str, kIN, PSfld, PSmom);

      // @label:externalps_interpolators
      MomentumInterpolator<T, NDim> PSinterpFld(kIN, PSfld, f.getToolBox(), f.getKIR());
      MomentumInterpolator<T, NDim> PSinterpMom(kIN, PSmom, f.getToolBox(), f.getKIR());
      // @endlabel

      if (PSType == 1) {
        auto fLeft = getInputFluctuationsTypeI<Model, NDim>(PSinterpFld, model, f, "Random left", kCutOff);
        auto fRight = getInputFluctuationsTypeI<Model, NDim>(PSinterpFld, model, f, "Random right", kCutOff);

        f.inFourierSpace() = (fLeft + fRight) / sqrt(T(2));
        f.inFourierSpace().setZeroMode(0);

        auto pLeft = getInputFluctuationsTypeI<Model, NDim>(PSinterpMom, model, p, "Random left", kCutOff);
        auto pRight = getInputFluctuationsTypeI<Model, NDim>(PSinterpMom, model, p, "Random right", kCutOff);

        p.inFourierSpace() = pow(model.aI, model.alpha) * (pLeft + pRight) / sqrt(T(2));
        p.inFourierSpace().setZeroMode(0);
      } else {
        // @label:externalps_sample_field_modes
        auto fLeft = getInputFluctuationsTypeII<Model, NDim>(PSinterpFld, model, f, "Random left", kCutOff);
        auto fRight = getInputFluctuationsTypeII<Model, NDim>(PSinterpFld, model, f, "Random right", kCutOff);

        f.inFourierSpace() = (fLeft + fRight) / sqrt(T(2));
        f.inFourierSpace().setZeroMode(0);
        // @endlabel

        // @label:externalps_sample_momentum_modes
        auto pLeft = getInputFluctuationsTypeII<Model, NDim>(PSinterpMom, model, p, "Random left", kCutOff);
        auto pRight = getInputFluctuationsTypeII<Model, NDim>(PSinterpMom, model, p, "Random right", kCutOff);

        p.inFourierSpace() = pow(model.aI, model.alpha) * (pLeft + pRight) / sqrt(T(2));
        p.inFourierSpace().setZeroMode(0);
        // @endlabel
      }
    }

    template <class Model, size_t NDim>
    auto getInputFluctuationsNormTypeI(Model &model, Field<T, NDim> f, T kCutOff) const
    {
      if (NDim != 3)
        throw(PSTypeINotSupportedForNDIMDifferentFrom3("Initialization with PSType = 1 is only supported for NDim=3."));

      FourierSite<NDim> ntilde(f.getToolBox());
      auto k = ntilde.norm() * f.getKIR();
      auto Hcut = heaviside(kCutOff - k);

      auto invmult = MomentumMultiplicity<T, NDim>(f.getToolBox());

      return Hcut * safeSqrt(T(4) * Constants::pi<T> * pow(ntilde.norm(), T(2)) * invmult) *
             (T(1) / pow(model.aI, T(1.5))) * pow(lSide / pow<2>(f.getDx()), T(1.5)) * sqrt(T(0.5));
    }

    template <class Model, size_t NDim>
    auto getInputFluctuationsTypeI(Model &model, Field<T, NDim> f, const std::string &mySeed, T kCutOff) const
    {
      auto norm = getInputFluctuationsNormTypeI<Model, NDim>(model, f, kCutOff);
      return norm * RandomGaussianField<T, NDim>(baseSeed + mySeed + f.toString(), f.getToolBox());
    }

    template <class Model, size_t NDim>
    auto ConjugateGaussianFluctuationsVarianceOne(Model &model, Field<T, NDim> f, Field<T, NDim> p, T kCutOff) const
    {
      auto fLeft = getInputFluctuationsTypeI<Model, NDim>(model, f, "Random left", kCutOff);
      auto fRight = getInputFluctuationsTypeI<Model, NDim>(model, f, "Random right", kCutOff);

      f.inFourierSpace() = (fLeft + fRight) / sqrt(T(2));
      f.inFourierSpace().setZeroMode(0);

      auto pLeft = getInputFluctuationsTypeI<Model, NDim>(model, p, "Random left", kCutOff);
      auto pRight = getInputFluctuationsTypeI<Model, NDim>(model, p, "Random right", kCutOff);

      p.inFourierSpace() = (pLeft + pRight) / sqrt(T(2));
      p.inFourierSpace().setZeroMode(0);
    }

    template <class Model, class VF>
    void BunchDavisTransverseU1(Model &model, VF f, VF p, VF tmpF, VF tmpP, T aDot, T kCutOff) const
    {

      size_t N = GetNGrid::get(model); // Reads N (number of points per dimension)
      FourierSite<Model::NDim> ntilde(model.getToolBox());
      auto knorm = ntilde.norm() * model.fldU1(0_c)(1_c).getKIR();

      auto keffm = MakeVector(i, 1, Model::NDim, sin(Constants::pi<T> / N * ntilde(i)));
      auto keffm2 = Total(i, 1, Model::NDim, pow<2>(keffm(i)));

      // definition of the phase shifts to properly have gauge fields living on the links when doing fourier
      // transformations
      auto extraphaseminus = MakeVector(i, 1, Model::NDim, complexPhase(-Constants::pi<T> / N * ntilde(i)));
      auto extraphase = MakeVector(i, 1, Model::NDim, complexPhase(Constants::pi<T> / N * ntilde(i)));

      /* Initialization of Gauge Field */

      // We set the RGF fluctuations with variance = 1 for the gauge and electric field
      ForLoop(i, 1, Model::NDim, ConjugateGaussianFluctuationsVarianceOne(model, f(i), p(i), kCutOff););

      ForLoop(k, 1, Model::NDim, f(k).inFourierSpace() = asFourier(extraphaseminus(k)) * f(k).inFourierSpace();
              f(k).inFourierSpace().setZeroMode(0);
              p(k).inFourierSpace() = asFourier(extraphaseminus(k)) * p(k).inFourierSpace();
              p(k).inFourierSpace().setZeroMode(0););
      // Saves the transverse projected vector field in a auxiliar field
      ForLoop(k, 1, Model::NDim,
              tmpF(k).inFourierSpace() = Total(
                  j, 1, Model::NDim, ((k == j ? T(1) : T(0)) - keffm(k) * keffm(j) / keffm2) * f(j).inFourierSpace());
              tmpP(k).inFourierSpace() = Total(
                  j, 1, Model::NDim, ((k == j ? T(1) : T(0)) - keffm(k) * keffm(j) / keffm2) * p(j).inFourierSpace()););

      auto phaseA = BDPhasePi2A<T, Model::NDim>(model.getToolBox(), model.kIR, model.aDotI);
      auto phaseE = BDPhasePi2E<T, Model::NDim>(model.getToolBox(), model.kIR, model.aDotI);

      /* Initialization of Gauge Field */
      // Sets the purely transverse initial vector field with its Bunch-Davis amplitude and appropiate phase
      //
      ForLoop(k, 1, Model::NDim, f(k).inFourierSpace() = sqrt(T(1) / (T(2) * knorm)) * phaseA * tmpF(k).inFourierSpace();
              f(k).inFourierSpace() = asFourier(extraphase(k)) * f(k).inFourierSpace();
              f(k).inFourierSpace().setZeroMode(0);

              // f(k).inFourierSpace() = heaviside(kCutOff - knorm)*f(k).inFourierSpace();
      );

      /* Initialization of Electric Field */
      //
      ForLoop(k, 1, Model::NDim, p(k).inFourierSpace() = sqrt(knorm / T(2)) * phaseE * tmpP(k).inFourierSpace();
              p(k).inFourierSpace() = asFourier(extraphase(k)) * p(k).inFourierSpace();
              p(k).inFourierSpace().setZeroMode(0);

              // p(k).inFourierSpace() = heaviside(kCutOff - knorm)*p(k).inFourierSpace();
      );
    }

    /** @brief Reads a chiral gauge-field spectrum file with five columns per line:
     *
     *     k    P_A+(k)    P_E+(k)    P_A-(k)    P_E-(k)
     *
     * Empty lines and lines starting with '#' are skipped; k must be strictly increasing.
     * If physicalUnits is false (default), the columns are taken to be already in program units: k in units of
     * omegaStar and the spectra normalised such that the Bunch-Davies vacuum reads P_A = 1/(2k), P_E = k/2.
     * If physicalUnits is true, the columns are in GeV units and are converted as for the scalar external spectra:
     * k -> k/omegaStar, P_A -> omegaStar P_A, P_E -> P_E/omegaStar.
     */
    template <class Model>
    void readSpectrumFileChiralU1(Model &model, const std::string &str, bool physicalUnits, std::vector<T> &kIN,
                                  std::vector<T> &PAplus, std::vector<T> &PEplus, std::vector<T> &PAminus,
                                  std::vector<T> &PEminus) const
    {
      std::ifstream t(str);
      if (!t.good())
        throw(ExternalGaugePSFileError("Could not open the external gauge-field power spectrum file '" + str + "'."));

      kIN.clear();
      PAplus.clear();
      PEplus.clear();
      PAminus.clear();
      PEminus.clear();

      std::string line;
      size_t lineNumber = 0;
      while (std::getline(t, line)) {
        ++lineNumber;
        const auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#') continue;

        std::istringstream ss(line);
        T k, pap, pep, pam, pem;
        if (!(ss >> k >> pap >> pep >> pam >> pem))
          throw(ExternalGaugePSFileError("Line " + std::to_string(lineNumber) + " of '" + str +
                                         "' does not contain the 5 columns k, P_A+, P_E+, P_A-, P_E-."));
        if (!kIN.empty() && k <= kIN.back())
          throw(ExternalGaugePSFileError("The momenta in '" + str + "' must be strictly increasing (line " +
                                         std::to_string(lineNumber) + ")."));

        if (physicalUnits) {
          k /= model.omegaStar;
          pap *= model.omegaStar;
          pam *= model.omegaStar;
          pep /= model.omegaStar;
          pem /= model.omegaStar;
        }
        kIN.push_back(k);
        PAplus.push_back(pap);
        PEplus.push_back(pep);
        PAminus.push_back(pam);
        PEminus.push_back(pem);
      }
      if (kIN.size() < 3)
        throw(ExternalGaugePSFileError("The external gauge-field power spectrum file '" + str +
                                       "' must contain at least 3 data lines."));
    }

    /** @brief Initialises a U(1) gauge field and its electric field from external spectra of the two helicities.
     *
     * Port of the axion-inflation initialisation of the former implementation (axionu1initializer.h):
     *
     *   A_i(k) = e^{i pi n_i/N} [ sqrt(P_A+) e^{ik/aH}              eps+_i(k) xi_+(k)
     *                           + sqrt(P_A-) e^{ik/aH}              eps-_i(k) xi_-(k) ]
     *   E_i(k) = e^{i pi n_i/N} [ sqrt(P_E+) e^{i(k/aH - pi/2 s(k))} eps+_i(k) xi_+(k)
     *                           + sqrt(P_E-) (-i) e^{ik/aH}         eps-_i(k) xi_-(k) ]
     *
     *  - xi_+/- are two independent Gaussian random fields, normalised per |k| shell exactly like in
     *    BunchDavisTransverseU1 (so that a Bunch-Davies input file reproduces that initialisation statistically).
     *    The SAME xi_+/- are used for A and E: each mode of A and E comes from a single random amplitude, as in the
     *    former implementation, and only the spectral amplitude and the phase differ between them.
     *  - eps_+/- are the helicity vectors transverse to the lattice momentum sin(pi n/N) (see helicitybasis.h),
     *    and e^{i pi n_i/N} places the components on the links, so the Gauss constraint holds exactly.
     *  - s(k) interpolates between amplified (k << kExc) and Bunch-Davies (k >> kExc) positive-helicity modes,
     *    see ExcitedBDPhaseE. The negative helicity always obeys the Bunch-Davies relation.
     *  - Modes with k > kCutOff are set to zero.
     *
     * tmp is a scratch vector field (only its first two components are used, to store xi_+ and xi_-).
     */
    template <class Model, class VF>
    void ExternalPSChiralU1(Model &model, VF f, VF p, VF tmp, const std::string &file, bool physicalUnits, T aDot,
                            T kExc, T vExc, T kCutOff) const
    {
      static_assert(Model::NDim == 3, "ExternalPSChiralU1 initial conditions are only implemented for NDim = 3.");
      constexpr size_t NDim = Model::NDim;

      std::vector<T> kIN, PAplus, PEplus, PAminus, PEminus;
      readSpectrumFileChiralU1(model, file, physicalUnits, kIN, PAplus, PEplus, PAminus, PEminus);

      auto toolBox = model.getToolBox();
      const T kIR = f(1_c).getKIR();

      MomentumInterpolator<T, NDim> PSAplus(kIN, PAplus, toolBox, kIR);
      MomentumInterpolator<T, NDim> PSEplus(kIN, PEplus, toolBox, kIR);
      MomentumInterpolator<T, NDim> PSAminus(kIN, PAminus, toolBox, kIR);
      MomentumInterpolator<T, NDim> PSEminus(kIN, PEminus, toolBox, kIR);

      // 1. Two independent Gaussian random fields, one per helicity, shared by A and E.
      //    The seed is tagged with the gauge field name (not the scratch field), so that the realisation does not
      //    depend on which scratch memory is used.
      const std::string tag = "ChiralU1 " + f(1_c).toString();
      auto norm = getInputFluctuationsNormTypeI<Model, NDim>(model, tmp(1_c), kCutOff);
      tmp(1_c).inFourierSpace() = norm * RandomGaussianField<T, NDim>(baseSeed + tag + " plus", toolBox);
      tmp(1_c).inFourierSpace().setZeroMode(0);
      tmp(2_c).inFourierSpace() = norm * RandomGaussianField<T, NDim>(baseSeed + tag + " minus", toolBox);
      tmp(2_c).inFourierSpace().setZeroMode(0);

      // 2. Phases: half-link shift, Bunch-Davies phases and the excited-mode phase of E+.
      const size_t N = GetNGrid::get(model);
      FourierSite<NDim> ntilde(toolBox);
      auto extraphase = MakeVector(i, 1, NDim, complexPhase(Constants::pi<T> / N * ntilde(i)));

      auto phaseA = BDPhasePi2A<T, NDim>(toolBox, kIR, aDot);
      auto phaseEminus = BDPhasePi2E<T, NDim>(toolBox, kIR, aDot);
      auto phaseEplus = ExcitedBDPhaseE<T, NDim>(toolBox, kIR, aDot, kExc, vExc);

      // 3. Project on the helicity basis, apply the spectra and shift to the links.
      using EpsPlus = HelicityVectorPlus<T, NDim>;
      using EpsMinus = HelicityVectorMinus<T, NDim>;
      ForLoop(i, 1, NDim, auto epsPlus = EpsPlus(toolBox, int(i)); auto epsMinus = EpsMinus(toolBox, int(i));

              f(i).inFourierSpace() = safeSqrt(PSAplus) * phaseA * epsPlus * tmp(1_c).inFourierSpace() +
                                      safeSqrt(PSAminus) * phaseA * epsMinus * tmp(2_c).inFourierSpace();
              f(i).inFourierSpace() = asFourier(extraphase(i)) * f(i).inFourierSpace();
              f(i).inFourierSpace().setZeroMode(0);

              p(i).inFourierSpace() = safeSqrt(PSEplus) * phaseEplus * epsPlus * tmp(1_c).inFourierSpace() +
                                      safeSqrt(PSEminus) * phaseEminus * epsMinus * tmp(2_c).inFourierSpace();
              p(i).inFourierSpace() = asFourier(extraphase(i)) * p(i).inFourierSpace();
              p(i).inFourierSpace().setZeroMode(0););
    }

    std::string getBaseSeed() const { return baseSeed; }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */

    T lSide;
    std::string baseSeed;
  };
} // namespace TempLat

#endif
