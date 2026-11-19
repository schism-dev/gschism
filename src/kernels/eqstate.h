
// Equation of state (eqstate). Returns full density [kg/m^3].
//    T: temperature [C], S: salinity [PSU], z: depth of the point [m] (only used if ieos_pres=1)
// Branches are resolved at compile time from ieos_type/ieos_pres.
// T and S are clamped to [tempmin,tempmax] and [saltmin,saltmax] (UNESCO validity range), as in Fortran.

constexpr double tempmin = -2.0;  // [C]
constexpr double tempmax = 40.0;  // [C]
constexpr double saltmin =  0.0;  // [PSU]
constexpr double saltmax = 42.0;  // [PSU]

__inline__ __device__ double eqstate(double T, double S, const double z) {
   assert(T >= -98.0 && S >= -98.0 && "eqstate: impossible dry (T or S is a -99 sentinel)");

   T = fmax(tempmin, fmin(T, tempmax));
   S = fmax(saltmin, fmin(S, saltmax));

   if constexpr(ieos_type == 1) {
      return eos_b + eos_a*T;
   }
   else if constexpr(ieos_type == 0) {
      // UNESCO 1980 (Millero & Poisson 1981), one-atmosphere density
      const double sqrtS = sqrt(S);
      const double rhow = 999.842594 + T*(6.793952e-2 + T*(-9.095290e-3 + T*(1.001685e-4
                        + T*(-1.120083e-6 + T*6.536332e-9))));
      const double A    = 8.24493e-1 + T*(-4.0899e-3 + T*(7.6438e-5 + T*(-8.2467e-7 + T*5.3875e-9)));
      const double B    = -5.72466e-3 + T*(1.0227e-4 + T*(-1.6546e-6));
      const double rhoS = rhow + S*(A + B*sqrtS + 4.8314e-4*S);

      assert(rhoS >= 980.0 && "eqstate: weird density");

      if constexpr(ieos_pres == 0) {
         return rhoS;
      }
      else {
         // Secant bulk modulus K(S,T,p), p in bars
         const double p  = 1.0e-5*rho0*grav*fabs(z);
         const double Kw = 19652.21 + T*(148.4206 + T*(-2.327105 + T*(1.360477e-2 + T*(-5.155288e-5))));
         const double K0 = Kw + S*(54.6746 + T*(-0.603459 + T*(1.09987e-2 + T*(-6.1670e-5))))
                              + S*sqrtS*(7.944e-2 + T*(1.6483e-2 + T*(-5.3009e-4)));
         const double Ap = 3.239908 + T*(1.43713e-3 + T*(1.16092e-4 + T*(-5.77905e-7)))
                         + S*(2.2838e-3 + T*(-1.0981e-5 + T*(-1.6078e-6))) + 1.91075e-4*S*sqrtS;
         const double Bp = 8.50935e-5 + T*(-6.12293e-6 + T*5.2787e-8)
                         + S*(-9.9348e-7 + T*(2.0816e-8 + T*9.1697e-10));
         const double K  = K0 + p*(Ap + p*Bp);
         return rhoS/(1.0 - p/K);
      }
   }
   else {
      static_assert(ieos_type == 0 || ieos_type == 1, "eqstate: ieos_type not implemented");
   }

}