"""Check of the gravitational-wave induced MHD waves of Duez, Liu, Shapiro &
Stephens (2005, PRD 72, 024028, Sect. IV D 1, "Paper I") against the analytic
solution of Duez et al. (2005, PRD 72, 024029, "Paper II", Table I, case 2).

A linear standing gravitational wave h+,x = h0+,x sin(kz) cos(kt) (TT gauge,
geodesic slicing) in a homogeneous fluid at rest with a uniform field B0
excites Alfven, slow and fast magnetosonic waves plus a particular solution
at the frequency of the gravitational wave. The paper's z axis is x here,
(x,y,z)_paper = (y,z,x)_code, and B is Gaussian in the paper,
B_code = B_paper/sqrt(4 pi).

The linear solution depends on the equation of state only through the
enthalpy density rho0 h = e+p and the sound speed cs^2 = Gamma P0/(rho0 h),
so the Gamma=4/3 fluid of the paper is the p = cs0^2 e fluid of this code
with e0 = rho0 h/(1+cs0^2); delta P = cs0^2 delta e in both.

The perturbations are projected on their exact spatial dependence
(sin kz for delta P and delta B, cos kz for delta v) and evaluated at z=1/8,
the point of the paper's figures. Reported are the maximal deviations from
the analytic solution over the run, relative to the maximal amplitude of
the analytic solution, for
  dP     pressure perturbation (paper Figs. 13, 14)
  dvA    velocity along the Alfven eigenvector u_A = k x v_A (Fig. 15)
  dvm2   velocity along the slow eigenvector u_m2 (Fig. 16)
  dB     magnetic field perturbation (all components)
  h+     gravitational wave, (gamma_yy-gamma_zz)/2 (paper: gamma_xx)
The analytic solution neglects terms O(h0^2) and the self-gravity of the
fluid, O(E t^2); the paper finds a slow secular drift of delta P from these
(their Eq. 80, Fig. 14), which is also present here.

usage: python check.py [--plot] [rundir ...]     (default: .)
"""
import sys, re, numpy as np


def param(run, name, default=None, files=('start.in', 'run.in'), index=0):
    """Value of a namelist parameter (index-th number of a list) in start.in or run.in."""
    for fn in files:
        txt = '\n'.join(l for l in open(run + '/' + fn) if not l.lstrip().startswith('!'))
        m = re.search(r'\b' + name + r"\s*=\s*'?([-0-9.eE+, ]+)", txt)
        if m: return float(m.group(1).split(',')[index])
    return default


def nxgrid(run):
    return int(re.search(r'nxgrid=(\d+)', open(run + '/src/cparam.local').read()).group(1))


def xgrid(run):
    """Cell-centred grid of a periodic Pencil Code run."""
    x0, L, n = param(run, 'xyz0'), param(run, 'Lxyz'), nxgrid(run)
    return x0 + (np.arange(n) + .5) * L / n


def yzaverages(run):
    """Binary data/yzaverages.dat (lwrite_avg1d_binary=T): times and {name: (nt,nx)}."""
    names = [l.split('(')[0].strip() for l in open(run + '/yzaver.in') if l.strip()]
    nx = nxgrid(run)
    raw = open(run + '/data/yzaverages.dat', 'rb').read()
    ts, bl, pos = [], [], 0
    while pos < len(raw):
        for tg in (ts, bl):
            n = np.frombuffer(raw[pos:pos + 4], 'i4')[0]; pos += 4
            tg.append(np.frombuffer(raw[pos:pos + n], 'f8').copy()); pos += n + 4
    t = np.array([x[0] for x in ts])
    A = np.array(bl).reshape(len(t), len(names), nx)
    return t, {nm: A[:, i] for i, nm in enumerate(names)}


def to_paper(v):
    """Code components (x,y,z) -> paper components (x,y,z) = code (y,z,x)."""
    return np.array([v[1], v[2], v[0]])


class Analytic:
    """Paper II, Sect. IV B (B0z /= 0), in the paper's frame and Gaussian units.
    The mode amplitudes A_A, A_m1, A_m2 follow from d_t dv(0) = k S~/E (their
    Eq. 59), which is equivalent to their Eqs. (103)-(113)."""

    def __init__(self, run):
        w, e0 = param(run, 'cs0')**2, param(run, 'e0')
        self.k = 2 * np.pi * param(run, 'nkgr') / param(run, 'Lxyz')
        self.hp, self.hx = param(run, 'hplus0'), param(run, 'hcross0')
        self.B0 = to_paper([param(run, 'B_ext', index=i) for i in range(3)]) * np.sqrt(4 * np.pi)
        self.rhoh = (1 + w) * e0                            # rho0 h = e + p
        self.GP0 = w * self.rhoh                            # Gamma P0
        self.E = self.rhoh + self.B0 @ self.B0 / (4 * np.pi)   # Eq. (22)
        self.vA = self.B0 / np.sqrt(4 * np.pi * self.E)     # Eq. (56)
        k, vA, Bx, By, Bz = self.k, self.vA, *self.B0
        kv = np.array([0, 0, k])
        cm2 = vA @ vA + w * (1 - vA @ vA)                   # Eq. (63)
        self.cm2 = cm2
        St = np.array([(Bx * self.hp + By * self.hx) * Bz,  # Eqs. (51)-(53)
                       (Bx * self.hx - By * self.hp) * Bz,
                       -(Bx**2 - By**2) * self.hp - 2 * Bx * By * self.hx]) / (4 * np.pi)
        # particular solution dv = vp cos(kz) sin(kt), Eq. (62)
        kA = kv @ vA
        M = (k**2 - kA**2) * np.eye(3) + np.outer(vA, kA * kv - k**2 * vA) \
            + np.outer(kv, kA * vA - cm2 * kv)
        self.vp = np.linalg.solve(M, k**2 * St / self.E)
        # eigenmodes, Eqs. (68)-(71)
        b = k**2 * cm2 + w * kA**2
        disc = np.sqrt(b**2 - 4 * k**2 * w * kA**2)
        om2 = np.array([(b + disc) / 2, (b - disc) / 2])    # fast, slow
        self.om = np.concatenate([[abs(kA)], np.sqrt(om2)])
        uA = np.cross(kv, vA)
        um = [vA + o2 * (1 - vA @ vA) / ((o2 - k**2) * kA) * kv for o2 in om2]
        self.U = np.array([uA, *um])                        # rows: Alfven, fast, slow
        self.A = np.linalg.solve((self.U * self.om[:, None]).T, k * St / self.E - k * self.vp)

    def amplitudes(self, t):
        """Amplitudes of dv (times cos kz), dP and dB (times sin kz) at times t."""
        t = np.asarray(t)[:, None]
        k, om, A, U, B0 = self.k, self.om, self.A, self.U, self.B0
        dv = np.sin(k * t) * self.vp + (A * np.sin(om * t)) @ U
        dP = self.GP0 * k * ((1 - np.cos(k * t))[:, 0] * self.vp[2] / k
                             + ((A * U[:, 2] / om) * (1 - np.cos(om * t))).sum(1))
        modes = np.array([u[2] * B0 - B0[2] * u for u in U])          # (u_z B0 - B0z u)
        dB = k * ((1 - np.cos(k * t)) / k * (self.vp[2] * B0 - B0[2] * self.vp)
                  + (A / om * (1 - np.cos(om * t))) @ modes)
        return dv, dP, dB


def project(f, x, k, fun):
    """Amplitude of the fun(kx) Fourier mode of profiles f (nt, nx)."""
    return 2 * (f * fun(k * x)).mean(axis=1)


def analyse(run, plot=False):
    an = Analytic(run)
    t, av = yzaverages(run)
    x, k = xgrid(run), an.k
    w, e0 = param(run, 'cs0')**2, param(run, 'e0')
    Bc = [param(run, 'B_ext', index=i) for i in range(3)]
    dv = to_paper([project(av['u' + c + 'mx'], x, k, np.cos) for c in 'xyz'])
    dB = to_paper([project(av['b' + c + 'mx'] - Bc[i], x, k, np.sin)
                   for i, c in enumerate('xyz')]) * np.sqrt(4 * np.pi)
    dP = w * project(av['rhomx'] - e0, x, k, np.sin)
    hp = .5 * project(av['gyymx'] - av['gzzmx'], x, k, np.sin)
    hx = project(av['gyzmx'], x, k, np.sin)
    dva, dPa, dBa = an.amplitudes(t)
    c, s = np.cos(k / 8), np.sin(k / 8)                   # z = 1/8
    uA, um2 = an.U[0] / np.linalg.norm(an.U[0]), an.U[2] / np.linalg.norm(an.U[2])
    q = {'dP': (s * dP, s * dPa),
         'dvA': (c * uA @ dv, c * dva @ uA),
         'dvm2': (c * um2 @ dv, c * dva @ um2),
         'dB': (s * dB.T, s * dBa),
         'h+': (s * hp, s * an.hp * np.cos(k * t)),
         'hx': (s * hx, s * an.hx * np.cos(k * t))}
    err = {n: np.max(abs(a - b)) / np.max(abs(b)) for n, (a, b) in q.items()}
    if plot:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        h0, P0 = np.hypot(an.hp, an.hx), an.GP0 / (4 / 3)
        fig, ax = plt.subplots(3, 1, figsize=(6, 8), sharex=True)
        for a, (n, norm, lab) in zip(ax, [('dP', h0 * P0, r'$\delta P/(h_0 P_0)$'),
                                          ('dvA', h0, r'$\delta v_A/h_0$'),
                                          ('dvm2', h0, r'$\delta v_{m2}/h_0$')]):
            a.plot(t, q[n][1] / norm, 'k-', lw=2, label='analytic (Paper II)')
            a.plot(t, q[n][0] / norm, 'C1:', lw=1.5, label='Pencil Code')
            a.set_ylabel(lab)
        ax[0].legend(); ax[-1].set_xlabel(r'$t/T_{\rm GW}$')
        ax[0].set_title(f'Duez et al. (2005) GW-induced MHD waves, z=1/8, N={nxgrid(run)}')
        fig.tight_layout(); fig.savefig(run + '/gw_mhd_waves.png', dpi=120)
        print('  wrote', run + '/gw_mhd_waves.png')
    return t, err, an


args = [a for a in sys.argv[1:] if a != '--plot']
res = []
for run in args or ['.']:
    t, err, an = analyse(run, '--plot' in sys.argv)
    n = nxgrid(run)
    if not res:
        print(f"E={an.E:.4e}, v_A=({', '.join(f'{v:.4f}' for v in an.vA)}), c_m^2={an.cm2:.4f}, "
              f"omega/k: Alfven {an.om[0]/an.k:.4f}, fast {an.om[1]/an.k:.4f}, slow {an.om[2]/an.k:.4f}")
    print(f"{run}: N={n}, t_end={t[-1]:.2f} GW periods; max error / max |analytic| at z=1/8:")
    print('   ' + ', '.join(f'{q} {e:.2e}' for q, e in err.items()))
    res.append((n, err))
if len(res) > 1:
    for q in res[0][1]:
        o = [np.log(a[1][q] / b[1][q]) / np.log(b[0] / a[0]) for a, b in zip(res, res[1:])]
        print(f'convergence order {q:5s}:', ', '.join(f'{v:.2f}' for v in o))
