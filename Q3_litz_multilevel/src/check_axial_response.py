"""Independent Joule-loss integral versus the closed-form axial response."""
from pathlib import Path
import json
import numpy as np
from scipy.special import jv
from scipy.integrate import quad
from multipole_rom import axial_polarizability,MU,SIGMA

ROOT=Path(__file__).resolve().parents[1]

def check():
    rows=[]
    for diameter_um in [149.24,345.494149]:
        for frequency in [50000.,200000.,1000000.]:
            a=diameter_um*1e-6/2;w=2*np.pi*frequency
            k=(1-1j)*np.sqrt(w*MU*SIGMA/2)
            # H_applied has unit peak amplitude. J_phi=-dH_z/dr.
            loss_integral=np.pi/SIGMA*quad(lambda r:r*abs(k*jv(1,k*r)/jv(0,k*a))**2,0,a,epsabs=1e-24)[0]
            loss_alpha=.5*w*MU*(-axial_polarizability(a,frequency).imag)
            error=abs(loss_alpha/loss_integral-1)
            if error>1e-9:raise RuntimeError('Axial cylinder power identity failed')
            rows.append(dict(diameter_um=diameter_um,frequency_Hz=frequency,
                loss_for_unit_peak_axial_H_W_per_m=loss_integral,
                polarizability_loss_W_per_m=loss_alpha,relative_error=error))
    result=dict(status='analytic_identity_pass',fit_coefficients=None,rows=rows,
                scope='Local uniform axial magnetic field on an isolated infinite cylinder; not validation of the complete curved-wire ROM')
    (ROOT/'data/axial_response_verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result

if __name__=='__main__':print(json.dumps(check(),indent=2))
