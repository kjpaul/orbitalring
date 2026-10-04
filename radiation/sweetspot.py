"""Dose at a point below the cable: direction-by-direction sum over the sky.
Local frame: x east, y north, z up. The cable is a square bar 8.20 m on a side running east-west, bottom face a height h
above the point. The room is a tube along the ring: half-width a north-south, floor zf below the point, with a floor slab
and side walls of given areal density. Earth blocks cosmic rays within the limb; trapped protons arrive from all directions
allowed by their pitch angle (their gyration circles are tens of km across, they do not come from the ground)."""
import numpy as np
W = 8.20; RHO_CABLE = 1.70; C_TO_WATER = 0.89     # g/cm3; carbon to water equivalence from PSTAR ranges (7.72/8.63 at 100 MeV)
def sphere(n=60000):
    i = np.arange(n)+0.5; z = 1-2*i/n; ph = np.pi*(1+5**0.5)*i; r = np.sqrt(1-z*z)
    return np.stack([r*np.cos(ph), r*np.sin(ph), z], -1)
U = sphere()
def thickness(h, a=2.0, zf=1.0, t_wall=5.0, t_floor=5.0, t_ceil=0.0, cable=True, y0=0.0):
    """water-equivalent g/cm2 along each direction U; also returns a mask of directions that pass through the cable"""
    ux, uy, uz = U.T; t = np.zeros(len(U)); hit = np.zeros(len(U), bool)
    up = uz > 1e-9
    with np.errstate(divide="ignore", invalid="ignore"):
        y_at = y0 + h*uy/uz                                   # where the ray crosses the cable's bottom plane
        hit = up & (np.abs(y_at) < W/2) & cable
        y_top = y_at + W*uy/uz
        path = np.where(np.abs(y_top) <= W/2, W/uz, (W/2*np.sign(uy)-y_at)/uy)   # metres inside the bar
        t_c = np.where(hit, path*100*RHO_CABLE*C_TO_WATER + t_ceil/np.maximum(uz, 1e-9), 0.0)
        down = uz < -1e-9
        floor = down & (np.abs(y0 + zf*uy/(-uz)) < a)
        t_f = np.where(floor, t_floor/np.maximum(-uz, 1e-9), 0.0)
        wall = ~hit & ~floor
        t_w = np.where(wall, t_wall/np.maximum(np.abs(uy), 1e-6), 0.0)
    return t_c + t_f + t_w, hit
def earth_blocked(alt_km):
    return U[:,2] < -np.sqrt(1-(6378.137/(6378.137+alt_km))**2)
def stormer(alt_km, decl_deg=0.0, K=57.0):
    """Stormer cutoff rigidity (GV) at the magnetic equator for each arrival direction. K = 59.6 x (present dipole / 8.06e22)"""
    r = (6371.2+alt_km)/6371.2; d = np.radians(decl_deg)
    east = np.array([np.cos(d), -np.sin(d), 0.0])             # magnetic east (perpendicular to horizontal B)
    cg = U @ east
    return K/(r*r*(1+np.sqrt(1-cg))**2)
