"""Depth-dose response of a thick water slab to mono-energetic primaries (Geant4 via geant4_pybind).

A pencil beam enters a slab that is wide enough to contain the whole cascade. Energy deposited is
scored in depth bins over the full lateral extent, which by reciprocity equals the dose at a point
at that depth under a broad parallel beam, i.e. the dose at the centre of a solid sphere of that
radius under isotropic flux. Scoring is by particle class and kinetic-energy band (Geant4 command
scorers), so a quality factor Q(LET) can be applied afterwards.

usage: python g4_depthdose.py <physics list> <particle> <outdir> <nev:E_MeV_total> [more nev:E ...]
"""
import sys, os
from geant4_pybind import *

DEPTH_CM = float(os.environ.get("DEPTH_CM", 600.0))      # water, 1 g/cm3
NBIN = int(os.environ.get("NBIN", 300))
HALF_XY_M = 15.0

# kinetic-energy bands (MeV, total kinetic energy of the particle) for LET-dependent species
PBANDS = [0, 0.3, 1, 3, 10, 30, 100, 300, 1e9]
ABANDS = [0, 2, 8, 30, 100, 400, 1e9]

SCORERS = []   # (name, particle list or None, (Emin,Emax) or None)
SCORERS.append(("tot", None, None))
SCORERS.append(("em", ["e-", "e+", "gamma"], None))
SCORERS.append(("mupi", ["mu-", "mu+", "pi-", "pi+", "kaon-", "kaon+"], None))
for i in range(len(PBANDS)-1): SCORERS.append((f"p{i}", ["proton"], (PBANDS[i], PBANDS[i+1])))
SCORERS.append(("dt", ["deuteron", "triton"], None))
for i in range(len(ABANDS)-1): SCORERS.append((f"a{i}", ["alpha", "He3"], (ABANDS[i], ABANDS[i+1])))
SCORERS.append(("ion", ["GenericIon"], None))
KEEP = []

class Det(G4VUserDetectorConstruction):
    def Construct(self):
        nist = G4NistManager.Instance()
        vac = nist.FindOrBuildMaterial("G4_Galactic"); wat = nist.FindOrBuildMaterial("G4_WATER")
        w = G4Box("W", (HALF_XY_M+1)*m, (HALF_XY_M+1)*m, (DEPTH_CM/2+50)*cm)
        lw = G4LogicalVolume(w, vac, "W")
        pw = G4PVPlacement(None, G4ThreeVector(), lw, "W", None, False, 0)
        s = G4Box("S", HALF_XY_M*m, HALF_XY_M*m, DEPTH_CM/2*cm)
        ls = G4LogicalVolume(s, wat, "S")
        G4PVPlacement(None, G4ThreeVector(), ls, "S", lw, False, 0)
        l = G4Box("L", HALF_XY_M*m, HALF_XY_M*m, DEPTH_CM/2/NBIN*cm)
        self.ll = G4LogicalVolume(l, wat, "L")
        G4PVReplica("L", self.ll, ls, kZAxis, NBIN, DEPTH_CM/NBIN*cm)
        KEEP.extend([w, lw, s, ls, l, self.ll])
        return pw
    def ConstructSDandField(self):
        mfd = G4MultiFunctionalDetector("mfd")
        G4SDManager.GetSDMpointer().AddNewDetector(mfd)
        for name, parts, band in SCORERS:
            ps = G4PSEnergyDeposit(name)
            if parts:
                if band:
                    f = G4SDParticleWithEnergyFilter("f"+name, band[0]*MeV, band[1]*MeV)
                else:
                    f = G4SDParticleFilter("f"+name)
                for p in parts: f.add(p)
                ps.SetFilter(f); KEEP.append(f)
            mfd.RegisterPrimitive(ps); KEEP.append(ps)
        if particle.startswith("ion:"):
            Z, A = map(int, particle.split(":")[1:])
            ps = G4PSEnergyDeposit("prim"); f = G4SDParticleFilter("fprim"); f.addIon(Z, A)
            ps.SetFilter(f); mfd.RegisterPrimitive(ps); KEEP.extend([ps, f]); SCORERS.append(("prim", None, None))
        KEEP.append(mfd)
        self.SetSensitiveDetector(self.ll, mfd)

import numpy as np
ACC = {}
class Ev(G4UserEventAction):
    def __init__(self):
        super().__init__(); self.ids = None
    def EndOfEventAction(self, ev):
        sdm = G4SDManager.GetSDMpointer()
        if self.ids is None:
            self.ids = [(n, sdm.GetCollectionID("mfd/"+n)) for n, _, _ in SCORERS]
        hce = ev.GetHCofThisEvent()
        for n, cid in self.ids:
            hm = hce.GetHC(cid)
            a = ACC[n]
            for k, v in hm: a[k] += v / MeV

class Gen(G4VUserPrimaryGeneratorAction):
    def __init__(self):
        super().__init__(); self.gun = G4ParticleGun(1)
        self.gun.SetParticlePosition(G4ThreeVector(0, 0, -(DEPTH_CM/2+1)*cm))
        self.gun.SetParticleMomentumDirection(G4ThreeVector(0, 0, 1))
    def GeneratePrimaries(self, ev): self.gun.GeneratePrimaryVertex(ev)

class Act(G4VUserActionInitialization):
    def Build(self): self.SetUserAction(Gen()); self.SetUserAction(Ev())

plname, particle, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
jobs = [(int(a.split(":")[0]), float(a.split(":")[1])) for a in sys.argv[4:]]
os.makedirs(outdir, exist_ok=True)
rm = G4RunManagerFactory.CreateRunManager(G4RunManagerType.Serial)
rm.SetUserInitialization(Det())
import geant4_pybind as _g
pl = getattr(_g, plname)()
pl.SetVerboseLevel(0)
rm.SetUserInitialization(pl); rm.SetUserInitialization(Act())
ui = G4UImanager.GetUIpointer()
ui.ApplyCommand("/run/setCut 1 mm")
if plname == "FTFP_BERT": ui.ApplyCommand("/run/setCutForAGivenParticle e- 5 cm")
rm.Initialize()
c = ui.ApplyCommand
gun_particle = particle
for nev, E in jobs:
    for n,_,_ in SCORERS: ACC[n] = np.zeros(NBIN)
    if particle.startswith("ion:"):
        Z, A = map(int, particle.split(":")[1:]); c(f"/gun/particle ion"); c(f"/gun/ion {Z} {A} {Z}")
    else:
        c(f"/gun/particle {particle}")
    c(f"/gun/energy {E} MeV")  # total kinetic energy
    rm.BeamOn(nev)
    tag = f"{outdir}/{particle.replace(':','_')}_{E:.2f}MeV_{nev}"
    np.savez(tag+".npz", E=E, nev=nev, names=[n for n,_,_ in SCORERS], data=np.array([ACC[n] for n,_,_ in SCORERS]))
    print("DONE", tag, f"sum tot per primary {ACC['tot'].sum()/nev:.1f} MeV", flush=True)
