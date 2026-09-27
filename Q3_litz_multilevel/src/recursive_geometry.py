"""Recursive normal-plane stranding; geometry only, never a field solver.

Factors and pitches are listed from individual wires outwards. A child axis is
C_child(t) = C_parent(t) + r [cos(phi) e1_parent + sin(phi) e2_parent].
Each child is offset in its parent's *normal plane*, not a global xy plane.
The radial direction defines a material frame transported to the next level.
Pitch is the advance of the global axial parameter t per relative 360-degree
turn in that parent material frame; local arc-length pitch is derived.
"""
from dataclasses import dataclass
from functools import reduce
from math import lcm
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.spatial import cKDTree


def layout(n):
    """Minimum unit-circle-center separation = 1; circular bundle envelopes."""
    if n == 1:
        return np.zeros((1, 2))
    if n == 7:
        angles = np.arange(6) * np.pi / 3
        return np.vstack(([0., 0.], np.c_[np.cos(angles), np.sin(angles)]))
    if 2 <= n <= 6:
        angles = np.arange(n) * 2 * np.pi / n
        return np.c_[np.cos(angles), np.sin(angles)] / (2 * np.sin(np.pi / n))
    # Compact hexagonal arrangement. Whole shells are retained when possible;
    # residual points are deterministic. No symmetry is assumed for these sets.
    k = int(np.ceil(np.sqrt(n)))
    pts = np.array([(q + r/2, np.sqrt(3)*r/2) for q in range(-k, k+1) for r in range(-k, k+1)])
    angle = np.arctan2(pts[:, 1], pts[:, 0])
    ids = np.lexsort((angle, np.round((pts**2).sum(1), 10)))[:n]
    return pts[ids] - pts[ids].mean(0)


@dataclass
class RecursiveCable:
    factors: tuple
    pitches_mm: tuple
    copper_area_mm2: float = 6.
    gap_um: float = 20.
    samples_per_fastest_turn: int = 128

    def __post_init__(self):
        if len(self.factors) != len(self.pitches_mm) or min(self.factors) < 1:
            raise ValueError("One nonzero signed pitch per grouping level is required")
        self.n = int(np.prod(self.factors))
        self.a = np.sqrt(self.copper_area_mm2 * 1e-6 / (np.pi*self.n))
        gap = self.gap_um * 1e-6
        self.offsets = []
        env = self.a
        for factor in self.factors:
            xy = layout(factor) * (2*env + gap)
            self.offsets.append(xy)
            env += np.linalg.norm(xy, axis=1).max()
        self.envelope_radius = env
        pu = [round(abs(p)*1000) for p in self.pitches_mm]
        if min(pu) <= 0 or any(abs(round(abs(p)*1000)-abs(p)*1000) > 1e-7 for p in self.pitches_mm):
            raise ValueError("Pitch must be nonzero and specified to integer micrometres")
        self.period = reduce(lcm, pu) * 1e-6
        nt = max(129, int(np.ceil(self.period / (min(pu)*1e-6) * self.samples_per_fastest_turn)) + 1)
        if nt*self.n > 8_000_000:
            raise ValueError("Full-period sampling budget exceeded; use a verified shorter screw-period cell")
        self.t = np.linspace(0, self.period, nt)
        c = np.zeros((1, nt, 3)); c[0, :, 2] = self.t
        e1 = np.zeros_like(c); e1[:, :, 0] = 1
        e2 = np.zeros_like(c); e2[:, :, 1] = 1
        identity = [()]
        for level in range(len(self.factors)-1, -1, -1):
            xy = self.offsets[level]
            phi = (2*np.pi*self.t[None, :] / (self.pitches_mm[level]*1e-3)
                   + np.arctan2(xy[:, 1], xy[:, 0])[:, None])
            u = np.cos(phi)[None, :, :, None]*e1[:, None] + np.sin(phi)[None, :, :, None]*e2[:, None]
            child = c[:, None] + np.linalg.norm(xy, axis=1)[None, :, None, None]*u
            c = child.reshape(-1, nt, 3)
            u = u.reshape(-1, nt, 3)
            tangent = self._spline(c)(self.t, 1) + np.array([0., 0., 1.])
            tangent /= np.linalg.norm(tangent, axis=2, keepdims=True)
            e1 = u - np.sum(u*tangent, axis=2, keepdims=True)*tangent
            e1 /= np.linalg.norm(e1, axis=2, keepdims=True)
            e2 = np.cross(tangent, e1)
            identity = [parent + (j,) for parent in identity for j in range(len(xy))]
        self.xyz = c
        self.ids_outer_to_inner = identity
        self.spline = self._spline(c)

    def _spline(self, coordinates):
        v = coordinates - np.array([0., 0., 1.])[None, None, :]*self.t[None, :, None]
        v[:, -1] = v[:, 0]
        return CubicSpline(self.t, v, axis=1, bc_type="periodic", extrapolate="periodic")

    def evaluate(self, t, derivative=0):
        out = self.spline(t, derivative)
        if derivative == 0:
            out[..., 2] += t
        elif derivative == 1:
            out[..., 2] += 1
        return out

    def at_physical_z(self, z):
        """Monotone inverse of the axis coordinate, with Newton refinement."""
        z = np.atleast_1d(z)
        t = np.broadcast_to(z, (self.n, len(z))).copy()
        # Evaluate each strand's spline separately to avoid N-by-N expansion.
        result = np.empty((self.n, len(z), 3))
        tangents = np.empty_like(result)
        for j in range(self.n):
            cj = CubicSpline(self.t, self.xyz[j] - self.t[:, None]*[0., 0., 1.], bc_type="periodic", extrapolate="periodic")
            for _ in range(8):
                v = cj(t[j]); dz = cj(t[j], 1)[:, 2] + 1
                if np.min(dz) <= 0:
                    raise ValueError("Nonmonotone axial path; graph representation invalid")
                step = (v[:, 2]+t[j]-z)/dz
                t[j] -= step
                if np.max(abs(step)) < 1e-13:
                    break
            result[j] = cj(t[j]); result[j, :, 2] += t[j]
            tangents[j] = cj(t[j], 1) + [0., 0., 1.]
        tangents /= np.linalg.norm(tangents, axis=2, keepdims=True)
        return result, tangents, t

    def screw_mapping(self, cell_length, rotation):
        """Test actual planar cut positions and tangents, including identity map."""
        pos, tan, _ = self.at_physical_z([0., cell_length])
        co, si = np.cos(rotation), np.sin(rotation)
        rot = np.array([[co, -si, 0], [si, co, 0], [0, 0, 1]])
        end = pos[:, 1].copy(); end[:, 2] -= cell_length
        back = end @ rot
        dist, ids = cKDTree(pos[:, 0]).query(back)
        t_err = np.linalg.norm(tan[:, 1] @ rot - tan[ids, 0], axis=1)
        return dict(cell_length_m=float(cell_length), rotation_deg=float(np.rad2deg(rotation)),
                    strand_destination_to_source=ids.tolist(), permutation=bool(len(set(ids)) == self.n),
                    max_position_error_m=float(dist.max()), max_tangent_error=float(t_err.max()),
                    geometry_mapping_pass=bool(len(set(ids)) == self.n and dist.max()<1e-8 and t_err.max()<1e-5),
                    electrical_mapping_status="Not yet implemented or validated; geometry check alone is insufficient")

    def diagnostics(self, cross_sections=241):
        d1 = self.evaluate(self.t, 1)
        d2 = self.evaluate(self.t, 2)
        speed = np.linalg.norm(d1, axis=2)
        if np.min(d1[:, :, 2]) <= 0:
            raise ValueError("A strand reverses its axial direction")
        growth = np.trapz(speed, self.t, axis=1)/self.period
        theta = np.arctan2(np.linalg.norm(d1[:, :, :2], axis=2), d1[:, :, 2])
        curvature = np.linalg.norm(np.cross(d1, d2), axis=2)/speed**3
        xyz, _, _ = self.at_physical_z(np.linspace(0, self.period, cross_sections))
        mindist = min(cKDTree(xyz[:, k, :2]).query(xyz[:, k, :2], k=2)[0][:, 1].min()
                      for k in range(cross_sections))
        # A conservative projected tube radius at any global z plane.
        projected_radius = self.a/np.cos(theta.max())
        radial = np.linalg.norm(xyz[:, :, :2], axis=2)
        self.stats = dict(total_strands=self.n, grouping_inner_to_outer=list(self.factors),
          recursive_levels=len(self.factors), signed_relative_pitches_mm=list(self.pitches_mm),
          diameter_um=2*self.a*1e6, gap_design_um=self.gap_um,
          copper_area_normal_mm2=self.copper_area_mm2, full_identity_repeat_mm=self.period*1e3,
          outer_diameter_bound_mm=2*(radial.max()+projected_radius)*1e3,
          growth_min=float(growth.min()), growth_mean=float(growth.mean()), growth_max=float(growth.max()),
          total_copper_volume_per_axial_m_mm3=float(self.copper_area_mm2*1000*growth.mean()),
          Rdc_length_parallel_ohm_per_m=float(1/(5.8e7*np.pi*self.a**2*np.sum(1/growth))),
          max_tangent_angle_deg=float(np.rad2deg(theta.max())), min_bend_radius_mm=float(1/curvature.max()*1e3),
          min_sampled_axis_spacing_um=float(mindist*1e6),
          conservative_sampled_projected_clearance_um=float((mindist-2*projected_radius)*1e6),
          radial_excursion_min_mm=float(np.ptp(radial,axis=1).min()*1e3),
          radial_excursion_max_mm=float(np.ptp(radial,axis=1).max()*1e3),
          all_strands_reach_center=bool(np.all(radial.min(axis=1)<self.a)),
          note="Centerlines and geometric diagnostics only. Sampled clearance is not a rigorous continuous collision certificate. No electromagnetic optimum is claimed.")
        return self.stats


if __name__ == "__main__":
    root=Path(__file__).resolve().parents[1]
    folder=root/'geometry';folder.mkdir(exist_ok=True)
    for label,fac,pitches in [('candidate_343',(7,7,7),(10.,30.,90.)),('candidate_125_gentle',(5,5,5),(20.,60.,120.)),('pilot_9',(3,3),(4.,12.))]:
        c=RecursiveCable(fac,pitches,copper_area_mm2=6*9/343 if label=='pilot_9' else 6.)
        stat=c.diagnostics()
        # These cell choices are candidates for subsequent electrical validation.
        ell={343:5e-3,125:12e-3,9:4e-3/3}[c.n]
        stat['screw_mapping']=c.screw_mapping(ell,2*np.pi*ell/(pitches[-1]*1e-3))
        (folder/(label+'.json')).write_text(json.dumps(stat,indent=2),encoding='utf-8')
        np.savez_compressed(folder/(label+'_centerlines.npz'), t=c.t, xyz=c.xyz, radius=c.a, ids=c.ids_outer_to_inner)
        print(label,{k:v for k,v in stat.items() if k!='screw_mapping'},flush=True)
