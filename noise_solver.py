"""
Small-signal NOISE solver for the AFE, built on the same validated MNA used by
ac_solver/ac_mna.

Method (matches what Spectre does):
  1. DC solve -> operating point of every device (reuse ac_solver.ac_solve).
  2. Build the SAME 6-node small-signal Y matrix as the AC analysis. For noise
     the inputs vip/vin carry no signal, so the M7/M8 gates are AC ground; the
     Y matrix is therefore identical to the AC Y (vip/vin only ever sat in the
     AC RHS).
  3. Each transistor's channel produces a drain-source current-noise PSD
     S_id(f) = S_thermal + S_flicker(f)   [A^2/Hz]   (pmos_tft_model.get_noise_psd)
     Inject a unit current between that device's drain and source nodes, solve,
     and read the transimpedance to the differential output:
         Z_k(f) = (vop - von) / i_inject
     The device's output-noise contribution PSD is |Z_k|^2 * S_id_k.
  4. Sum over devices (uncorrelated) -> total differential output noise PSD.
  5. Refer to the amplifier input through the validated gain |H_amp(f)| to get
     the input-referred noise (IRN), the spec quantity.

Ground-truth check (Cadence Spectre, afe_gt/tb_noise.raw/noiseAnal.noise):
  total output 0.05-100 Hz = 2010 uVrms ; IRN = 209.6 uVrms ;
  M12=M13=47%, M14=M15=1.7%, M7=M8=1.1%, M9=M10=0.3%.
"""
import numpy as np
from pmos_tft_model import PMOS_TFT
from ac_mna import _stamp_mos, _stamp_adm
from ac_solver import ac_solve

# solved-node indices (same as ac_mna / ac_solver)
VOP, VON, VFBP, VFBN, NET20, NET2 = 0, 1, 2, 3, 4, 5
NN = 6
GND = ("v", 0.0)
VDDg = ("v", 0.0)
CL = 5e-12


def _n(i):
    return ("n", i)


def device_psd(W, L, Vs, Vd, Vg, freqs):
    """Drain-current noise PSD A^2/Hz over freqs: S_th + S_fl_1Hz/f."""
    t = PMOS_TFT(W=W, L=L)
    try:
        S_th, S_fl_1 = t.get_noise_psd(Vs, Vd, Vg, frequency=1.0)
    except Exception:
        return np.zeros_like(freqs), 0.0, 0.0
    return S_th + S_fl_1 / freqs, S_th, S_fl_1


def noise_analysis(sizes, bias, freqs):
    # ── 1. DC + small-signal params + gain (reuse the validated AC solver) ──
    ac = ac_solve(sizes, bias, freqs)
    if ac is None:
        return None
    dc = ac["dc_op"]
    ss = ac["ss"]
    Hmag = ac["gains"]                      # |vop-von|/vin_diff at each freq
    n2, vop, vfb, n20 = dc["net2"], dc["VOP"], dc["vfb"], dc["n20"]
    VCM, VDD, VB, VC = bias["VCM"], bias["VDD"], bias["VB"], bias["VC"]

    # per-device DC bias (Vs, Vd, Vg) — same mapping as ac_solver
    bpts = {
        "M6":  (VDD, n2,  VB),
        "M7":  (n2,  vop, VCM),
        "M8":  (n2,  vop, VCM),
        "M9":  (vop, 0.0, vfb),
        "M10": (vop, 0.0, vfb),
        "M11": (VDD, n20, VC),
        "M12": (n20, vfb, vop),
        "M13": (n20, vfb, vop),
        "M14": (vfb, 0.0, 0.0),
        "M15": (vfb, 0.0, 0.0),
    }
    # device (drain, gate, source) for Y stamping — inputs are AC ground
    devs = [
        ("M6",  _n(NET2),  GND,        VDDg),
        ("M7",  _n(VOP),   GND,        _n(NET2)),
        ("M8",  _n(VON),   GND,        _n(NET2)),
        ("M9",  GND,       _n(VFBP),   _n(VOP)),
        ("M10", GND,       _n(VFBN),   _n(VON)),
        ("M11", _n(NET20), GND,        VDDg),
        ("M12", _n(VFBN),  _n(VOP),    _n(NET20)),
        ("M13", _n(VFBP),  _n(VON),    _n(NET20)),
        ("M14", GND,       GND,        _n(VFBN)),
        ("M15", GND,       GND,        _n(VFBP)),
    ]
    # drain/source terminal of each device for current-noise injection
    inj = {name: (d, s) for name, d, g, s in devs}

    # per-device noise PSD
    psd = {}
    psd_split = {}
    for name in bpts_order(sizes):
        W, L = sizes[name]
        Vs, Vd, Vg = bpts[name]
        S, S_th, S_fl1 = device_psd(W, L, Vs, Vd, Vg, freqs)
        psd[name] = S
        psd_split[name] = (S_th, S_fl1)

    # ── 2/3. per-frequency: build Y, get transimpedance per device ──
    out_psd = np.zeros(len(freqs))                       # total output V^2/Hz
    dev_psd = {name: np.zeros(len(freqs)) for name in bpts}  # per-device output V^2/Hz

    for fi, f in enumerate(freqs):
        jw = 2j * np.pi * f
        Y = np.zeros((NN, NN), dtype=complex)
        RHS = np.zeros(NN, dtype=complex)  # unused for Y build
        for name, d, g, s in devs:
            p = ss[name]
            _stamp_mos(Y, RHS, d, g, s, p["gm"], p["gds"], p["Cgs"], p["Cgd"], jw)
        _stamp_adm(Y, RHS, _n(VOP), GND, jw * CL)
        _stamp_adm(Y, RHS, _n(VON), GND, jw * CL)

        Yinv = np.linalg.inv(Y)
        # transfer from injecting unit current at node j to (vop - von):
        #   t[j] = Yinv[VOP, j] - Yinv[VON, j]
        tvec = Yinv[VOP, :] - Yinv[VON, :]

        for name in bpts:
            d, s = inj[name]
            Z = 0.0 + 0.0j
            if d[0] == "n":
                Z += tvec[d[1]]
            if s[0] == "n":
                Z -= tvec[s[1]]
            contrib = (abs(Z) ** 2) * psd[name][fi]
            dev_psd[name][fi] = contrib
            out_psd[fi] += contrib

    # ── 4/5. integrate + input-refer ──
    return {
        "freqs": freqs,
        "out_psd": out_psd,          # differential output noise PSD V^2/Hz
        "dev_psd": dev_psd,          # per-device output PSD
        "Hmag": Hmag,                # |amplifier gain|
        "irn_psd": out_psd / Hmag ** 2,
        "psd_split": psd_split,
        "dc": dc,
    }


def bpts_order(sizes):
    return [m for m in ["M6","M7","M8","M9","M10","M11","M12","M13","M14","M15"] if m in sizes]


def band_rms(freqs, psd, f_lo, f_hi):
    mask = (freqs >= f_lo) & (freqs <= f_hi)
    return float(np.sqrt(np.trapezoid(psd[mask], freqs[mask])))


# ── run + compare to Cadence ground truth ────────────────────────────
if __name__ == "__main__":
    sizes = {
        "M6": (3000, 150), "M7": (25000, 150), "M8": (25000, 150),
        "M9": (12000, 500), "M10": (12000, 500),
        "M11": (300, 100), "M12": (500, 80), "M13": (500, 80),
        "M14": (2000, 500), "M15": (2000, 500),
    }
    bias = {"VDD": 40.0, "VCM": 32.0, "VB": 20.0, "VC": 26.0}
    freqs = np.logspace(-2, 4, 121)   # match Cadence grid (0.01..10k, dec=20)

    r = noise_analysis(sizes, bias, freqs)
    F_LO, F_HI = 0.05, 100.0

    Vout = band_rms(freqs, r["out_psd"], F_LO, F_HI)
    IRN = band_rms(freqs, r["irn_psd"], F_LO, F_HI)
    print(f"Python total OUTPUT noise {F_LO}-{F_HI} Hz = {Vout*1e6:.1f} uVrms  (Cadence 2010)")
    print(f"Python IRN               {F_LO}-{F_HI} Hz = {IRN*1e6:.1f} uVrms  (Cadence 209.6, spec<=44.5)")
    print(f"midband gain |H| = {r['Hmag'].max():.3f}")

    print("\nPer-device contribution (0.05-100 Hz):")
    var = {nm: band_rms(freqs, p, F_LO, F_HI)**2 for nm, p in r["dev_psd"].items()}
    tot = sum(var.values())
    cad = {"M12":47.0,"M13":47.0,"M14":1.7,"M15":1.7,"M7":1.1,"M8":1.1,"M9":0.3,"M10":0.3,"M6":0,"M11":0}
    for nm in sorted(var, key=lambda k:-var[k]):
        print(f"  {nm:<4} {np.sqrt(var[nm])*1e6:8.1f} uVrms  {var[nm]/tot*100:5.1f}%   (Cadence {cad.get(nm,0):.1f}%)")
