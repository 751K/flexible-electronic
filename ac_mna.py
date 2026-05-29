"""
Clean full-circuit small-signal MNA for the AFE.

No half-circuit approximation: the tail node net2 is a solved node, so the
differential virtual-ground behaviour falls out of the solve instead of being
assumed. One uniform transistor stamp is applied to all 10 devices.

Terminal encoding: each transistor terminal is either
    ("n", i)  -> solved node index i
    ("v", x)  -> known AC voltage x (AC grounds use x=0.0; driven inputs use x=Vin)
"""
import numpy as np
from pmos_tft_model import PMOS_TFT

# Solved nodes
VOP, VON, VFBP, VFBN, NET20, NET2 = 0, 1, 2, 3, 4, 5
NNODES = 6
GND = ("v", 0.0)        # AC ground
VDD = ("v", 0.0)        # ideal supply -> AC ground


def ss_params(W, L, Vs, Vd, Vg):
    """Small-signal params at a DC operating point.

    gm/gds are the *terminal* transconductance/output-conductance, extracted by
    finite-differencing the full terminal current get_Idc (which solves the
    OTFT's internal contact-network nodes s1/d1). The channel gm from
    _eval_channel is NOT used: the contact resistance degenerates it, and
    Spectre's AC analysis sees the terminal value. Using channel gm under-
    estimates the gain by ~0.8 dB and pushes the -3dB BW high (70 vs 52 Hz).
    Caps Cgss/Cgdd are the terminal caps and match Cadence directly.
    """
    t = PMOS_TFT(W=W, L=L)
    h = 1e-3
    Id = lambda vs, vd, vg: t.get_Idc(vs, vd, vg)
    gm  = (Id(Vs, Vd, Vg + h) - Id(Vs, Vd, Vg - h)) / (2 * h)
    gds = (Id(Vs, Vd + h, Vg) - Id(Vs, Vd - h, Vg)) / (2 * h)
    Cgss, Cgdd = t.get_capacitances(Vs, Vd, Vg)
    return gm, gds, Cgss, Cgdd


def _stamp_adm(Y, RHS, P, Q, y):
    """Two-terminal admittance y between terminals P and Q."""
    if P[0] == "n":
        Y[P[1], P[1]] += y
        if Q[0] == "n":
            Y[P[1], Q[1]] -= y
        else:
            RHS[P[1]] += y * Q[1]
    if Q[0] == "n":
        Y[Q[1], Q[1]] += y
        if P[0] == "n":
            Y[Q[1], P[1]] -= y
        else:
            RHS[Q[1]] += y * P[1]


def _stamp_vccs(Y, RHS, d, g, s, gm):
    """Transconductance VCCS: drain current i_d = gm*(Vg-Vs) flows from drain
    to source through the device. Canonical MNA stamp:

        Y[d,g] += gm   Y[d,s] -= gm
        Y[s,g] -= gm   Y[s,s] += gm

    Terminals that are known AC voltages move to the RHS as -c*Vknown.
    """
    def addrow(node, term, c):
        # adds c*term to the KCL equation of `node`
        if node[0] != "n":
            return
        if term[0] == "n":
            Y[node[1], term[1]] += c
        else:                       # known AC voltage -> move to RHS
            RHS[node[1]] -= c * term[1]
    addrow(d, g, +gm); addrow(d, s, -gm)
    addrow(s, g, -gm); addrow(s, s, +gm)


def _stamp_mos(Y, RHS, d, g, s, gm, gds, Cgs, Cgd, jw):
    # gds between drain and source
    _stamp_adm(Y, RHS, d, s, gds)
    # Cgs between gate and source, Cgd between gate and drain
    _stamp_adm(Y, RHS, g, s, jw * Cgs)
    _stamp_adm(Y, RHS, g, d, jw * Cgd)
    # transconductance
    _stamp_vccs(Y, RHS, d, g, s, gm)


def ac_solve_full(sizes, dc, bias, freqs, vin_diff=1.0):
    """
    sizes: {name:(W,L)}
    dc: solved DC node voltages {net2, VOP, vfb, net20}
    bias: {VDD, VCM, VB, VC}
    Returns gains array, dc gain dB, -3dB BW.
    """
    VCM = bias["VCM"]; VDDv = bias["VDD"]; VB = bias["VB"]; VC = bias["VC"]
    n2 = dc["net2"]; vop = dc["VOP"]; vfb = dc["vfb"]; n20 = dc["net20"]

    # small-signal params at the DC operating point of each device
    P = {}
    P["M6"]  = ss_params(*sizes["M6"],  VDDv, n2,  VB)
    P["M7"]  = ss_params(*sizes["M7"],  n2,   vop, VCM)
    P["M8"]  = ss_params(*sizes["M8"],  n2,   vop, VCM)
    P["M9"]  = ss_params(*sizes["M9"],  vop,  0.0, vfb)
    P["M10"] = ss_params(*sizes["M10"], vop,  0.0, vfb)
    P["M11"] = ss_params(*sizes["M11"], VDDv, n20, VC)
    P["M12"] = ss_params(*sizes["M12"], n20,  vfb, vop)
    P["M13"] = ss_params(*sizes["M13"], n20,  vfb, vop)
    P["M14"] = ss_params(*sizes["M14"], vfb,  0.0, 0.0)
    P["M15"] = ss_params(*sizes["M15"], vfb,  0.0, 0.0)

    CL = 5e-12
    half = vin_diff / 2.0
    VIP = ("v", +half); VIN = ("v", -half)

    gains = []
    for f in freqs:
        jw = 2j * np.pi * f
        Y = np.zeros((NNODES, NNODES), dtype=complex)
        RHS = np.zeros(NNODES, dtype=complex)

        n = lambda i: ("n", i)
        # device: (drain, gate, source)
        devs = [
            ("M6",  n(NET2),  VB_node(VB), VDD),
            ("M7",  n(VOP),   VIP,         n(NET2)),
            ("M8",  n(VON),   VIN,         n(NET2)),
            ("M9",  GND,      n(VFBP),     n(VOP)),
            ("M10", GND,      n(VFBN),     n(VON)),
            ("M11", n(NET20), VC_node(VC), VDD),
            ("M12", n(VFBN),  n(VOP),      n(NET20)),
            ("M13", n(VFBP),  n(VON),      n(NET20)),
            ("M14", GND,      GND,         n(VFBN)),
            ("M15", GND,      GND,         n(VFBP)),
        ]
        for name, d, g, s in devs:
            gm, gds, Cgs, Cgd = P[name]
            _stamp_mos(Y, RHS, d, g, s, gm, gds, Cgs, Cgd, jw)

        # load caps on the two outputs
        _stamp_adm(Y, RHS, n(VOP), GND, jw * CL)
        _stamp_adm(Y, RHS, n(VON), GND, jw * CL)

        V = np.linalg.solve(Y, RHS)
        gains.append(abs(V[VOP] - V[VON]) / vin_diff)

    gains = np.array(gains)
    Av0 = gains[0]
    Av0_dB = 20 * np.log10(max(Av0, 1e-12))
    peak = gains.max()
    a3 = peak / np.sqrt(2)
    # upper -3dB: first freq above the peak index where gain drops below a3
    ipk = int(np.argmax(gains))
    bw = freqs[-1]
    for i in range(ipk, len(gains)):
        if gains[i] < a3:
            bw = freqs[i]; break
    return {"gains": gains, "freqs": freqs, "Av0_dB": Av0_dB,
            "peak_dB": 20*np.log10(peak), "bw_Hz": bw, "params": P}


# gate of M6/M11 is a DC bias source -> AC ground
def VB_node(_): return ("v", 0.0)
def VC_node(_): return ("v", 0.0)


if __name__ == "__main__":
    sizes = {
        "M6": (3000, 150), "M7": (25000, 150), "M8": (25000, 150),
        "M9": (12000, 500), "M10": (12000, 500),
        "M11": (300, 100), "M12": (500, 80), "M13": (500, 80),
        "M14": (2000, 500), "M15": (2000, 500),
    }
    dc = {"net2": 37.9958, "VOP": 19.599, "vfb": 10.2871, "net20": 28.5842}
    bias = {"VDD": 40.0, "VCM": 32.0, "VB": 20.0, "VC": 26.0}
    freqs = np.logspace(-2, 4, 400)
    r = ac_solve_full(sizes, dc, bias, freqs)
    print(f"DC(0.01Hz) gain = {r['Av0_dB']:.2f} dB")
    print(f"peak gain      = {r['peak_dB']:.2f} dB   (Cadence 19.96 dB)")
    print(f"upper -3dB BW  = {r['bw_Hz']:.1f} Hz      (Cadence 52.3 Hz)")
