# sedov_taylor.py
#
# 30-oct-20
# Author: F. Gent (fred.gent.ncl@gmail.com).
#
import numpy as np
import os
from pencil import read


def sedov_taylor(*args, **kwargs):
    """
    Compute analytic radial time evolution of SN blast waves for
    comparison with numerical results

    *t_sedov*:
      Time_series object read from the simulation sn_series.dat

    *par*:
      Param object containing the simulation parameters

    *time*:
      list of time in code units

    *nt*:
      Integer size of analytic arrays

    *endt*
      Real end time in code units for the time series

    *dims*:
      Dimension of the simulation default 3D

    *rho0*:
      Ambient ISM density
    """
    st_tmp = SedovTaylor()
    st_tmp.get_st(*args, **kwargs)
    return st_tmp


class SedovTaylor(object):
    """
    SedovTaylor -- holds blast wave parameters and radial profiles.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.t = []
        self.keys = []

    def get_st(
        self,
        t_sedov=0,
        par=list(),
        time=list(),
        nt=5000,
        startt=0.0,
        endt=0.005,
        dims=3,
        quiet=True,
        rho0=1.6728e-24,
        M0=10,
        lsnowplough=True,
        lcioffi=True,
    ):
        """
        Compute analytic radial time evolution of SN blast waves for
        comparison with numerical results

        *t_sedov*:
          Time_series object read from the simulation sn_series.dat

        *par*:
          Param object containing the simulation parameters

        *time*:
          list of time in code units

        *nt*:
          Integer size of analytic arrays

        *endt*
          Real end time in code units for the time series

        *dims*:
          Dimension of the simulation default 3D

        *rho0*:
          Ambient ISM density

        *lsnowplough:
          Include original snowplough profile

        *lcioffi:
          Include Cioffi et al profile
        """

        if isinstance(par, list):
            unit_length = 3.086e21  # 1pc in cm
            unit_time = 3.1557e16  # Gyr in s
            unit_velocity = unit_length / unit_time  # ~km/s in cm/s
            unit_density = 1.6728e-24  # g/cm3
            unit_energy_density = unit_density * unit_velocity ** 2  # erg/cm3
            unit_energy = unit_energy_density * unit_length ** 3  # erg
            E0 = 1e51 / unit_energy  # erg
            M0 = 10
        else:
            unit_length = par.unit_length
            unit_time = par.unit_time
            unit_velocity = par.unit_velocity
            unit_density = par.unit_density
            unit_energy_density = par.unit_energy_density
            unit_energy = par.unit_energy
            E0 = par.ampl_sn
            if par.lsn_mass:
                M0 = par.mass_sn
            else:
                M0 = 10
        if len(time) > 0:
            time = np.array(time)
        else:
            time = np.linspace(startt, endt, nt) + t_sedov
        xi = 2.026  # McKee/Ostriker 1988
        rho0 /= unit_density
        setattr(self, "unit_density", unit_density)
        setattr(self, "unit_length", unit_length)
        setattr(self, "unit_time", unit_time)
        setattr(self, "unit_velocity", unit_velocity)
        setattr(self, "unit_energy", unit_energy)
        setattr(self, "E0", E0)
        setattr(self, "M0", M0)
        setattr(self, "rho0", rho0)
        setattr(self, "xi", xi)
        setattr(self, "t_sedov", t_sedov)
        setattr(self, "dims", dims)
        setattr(self, "t", time)
        for key in self.__dict__.keys():
            print(key, self.__getattribute__(key))

        # Sedov-Taylor
        rst = (self.xi * self.E0 / self.rho0) ** (1.0 / (self.dims + 2.0)) * self.t ** (
            2.0 / (self.dims + 2.0)
        )
        rdot = (
            2.0
            / (self.dims + 2.0)
            * (self.xi * self.E0 / self.rho0) ** (1.0 / (self.dims + 2.0))
            * self.t ** (2.0 / (self.dims + 2.0) - 1)
        )
        setattr(self, "stR", rst)
        setattr(self, "stRdot", rdot)
        m_p = 1.67262158e-24
        n0 = self.rho0 * self.unit_density / m_p
        setattr(self, "n0", n0)

        def get_snpl(self, quiet):
            """
            Compute analytic radial time evolution of SN blast waves for
            comparison of snowplough solution with numerical results
            """
            # Woltier transition to momentum conservation
            vtrans = 2.0e2 * self.n0 ** (2.0 / 17.0) * self.E0 ** (1.0 / 17.0)
            ttrans = (
                2.0
                / (2.0 + self.dims)
                * (self.xi * self.E0 / self.rho0) ** (1.0 / (self.dims + 2.0))
                * vtrans ** (-1)
            ) ** (1.0 / (1 - 2.0 / (self.dims + 2.0)))
            rtrans = (self.xi * self.E0 / self.rho0) ** (
                1.0 / (self.dims + 2.0)
            ) * ttrans ** (2.0 / (self.dims + 2.0))
            if not quiet:
                print("Woltier transition to momentum conservation")
                print("ttrans {}, rtrans {}, vtrans {}".format(ttrans, rtrans, vtrans))
            setattr(self, "tWolt", ttrans)
            setattr(self, "rWolt", rtrans)
            setattr(self, "vWolt", vtrans)
            # snowplough
            rsnpl = self.stR.copy()
            vsnpl = self.stRdot.copy()
            isnpl = np.where(self.t > self.tWolt)[0]
            rsnpl[isnpl] = self.rWolt * (
                8.0 / (self.dims + 2.0) * self.t[isnpl] / self.tWolt
                - 3.0 / (self.dims + 2.0)
            ) ** (1.0 / 4.0)
            vsnpl[isnpl] = self.vWolt * (
                8.0 / (self.dims + 2.0) * self.t[isnpl] / self.tWolt
                - 3.0 / (self.dims + 2.0)
            ) ** (-3.0 / 4.0)
            setattr(self, "snplR", rsnpl)
            setattr(self, "snplRdot", vsnpl)

        def get_cioffi(self, quiet):
            """
            Compute analytic radial time evolution of SN blast waves for
            comparison with numerical results
            """
            vej = (400 * self.E0 / self.M0) ** 0.5
            setattr(self, "vejecta", vej)
            # pressure driven snowplough transition
            tpds = (
                3.61e-5 * self.E0 ** (3.0 / 14.0) / self.n0 ** (4.0 / 7.0)
            ) / 2.718281828
            rpds = (self.xi * self.E0 / self.rho0) ** (
                1.0 / (self.dims + 2.0)
            ) * tpds ** (2.0 / (self.dims + 2.0))
            vpds = (2.0 / (2 + self.dims)) * (self.xi * self.E0 / self.rho0) ** (
                2.0 / (2 + self.dims) - 1
            )
            if not quiet:
                print("Pressure-driven snowplough transition")
                print("tpds {}, rpds {}, vpds {}".format(tpds, rpds, vpds))
            setattr(self, "tpds", tpds)
            setattr(self, "rpds", rpds)
            setattr(self, "vpds", vpds)
            # momentum conserving snowplough transition
            tmcs = (
                61
                * self.vejecta ** 3
                / self.n0 ** (3.0 / 7.0)
                / self.E0 ** (3.0 / 14.0)
                * self.tpds
            )
            rmcs = self.rpds * (4.0 / 3.0 * tmcs / self.tpds - 1.0 / 3.0) ** (
                3.0 / 10.0
            )
            if not quiet:
                print("Momentum conserving snowplough")
                print("tmcs {}, rmcs {}".format(tmcs, rmcs))
            setattr(self, "tmcs", tmcs)
            setattr(self, "rmcs", rmcs)
            # Cioffi et al
            rcioffi = self.stR.copy()
            vcioffi = self.stRdot.copy()
            icioffi = np.where(self.t > self.tpds)[0]
            rcioffi[icioffi] = self.rpds * (
                4.0 / 3.0 * self.t[icioffi] / self.tpds - 1.0 / 3.0
            ) ** (3.0 / 10.0)
            vcioffi[icioffi] = (
                0.3
                * self.rpds
                * 4.0
                / 3.0
                / self.tpds
                * (4.0 / 3.0 * self.t[icioffi] / self.tpds - 1.0 / 3.0) ** (-7.0 / 10.0)
            )
            jcioffi = np.where(self.t > self.tmcs)[0]
            rcioffi[jcioffi] = (
                self.rpds
                * (
                    4.66
                    * (self.t[jcioffi] - self.tmcs)
                    / self.tpds
                    * (1.0 - 0.779 / (self.tmcs / self.tpds) ** 0.17)
                    + (self.rmcs / self.rpds) ** 4
                )
                ** 0.25
            )
            vcioffi[jcioffi] = (
                0.25
                * self.rpds
                * 4.66
                / self.tpds
                * (1.0 - 0.779 / (self.tmcs / self.tpds) ** 0.17)
                * (
                    4.66
                    * (self.t[jcioffi] - self.tmcs)
                    / self.tpds
                    * (1.0 - 0.779 / (self.tmcs / self.tpds) ** 0.17)
                    + (self.rmcs / self.rpds) ** 4
                )
                ** (-0.75)
            )
            setattr(self, "CioffiR", rcioffi)
            setattr(self, "CioffiRdot", vcioffi)

        if lsnowplough:
            get_snpl(self, quiet)
        if lcioffi:
            get_cioffi(self, quiet)
        return self
