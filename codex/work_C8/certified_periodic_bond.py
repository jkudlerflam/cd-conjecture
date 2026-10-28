"""The C8 flat bond with a uniform, nonperturbative selector-defect gap.

Only the forbidden-bond coefficient changes. The legal restriction is
identical to PeriodicBond(flat_vacancies=True). See selector_repair.md.
"""
import numpy as np

from periodic_clean_bond import PeriodicBond


class CertifiedPeriodicBond(PeriodicBond):
    def __init__(self, q=4, epsilon_max=1.):
        if not 0 < epsilon_max <= 1:
            raise ValueError('epsilon_max must be in (0,1]')
        super().__init__(q=q, U=1., J0=1., flat_vacancies=True)
        self.epsilon_max = epsilon_max
        self.original_selector_penalty = self.Ksel
        self.full_perturbation_bound = (1+2*self.m+4*(self.k+2)
                                       +409600*np.pi**4*self.b**4*self.m**3
                                       +320*np.pi**2*self.b**2*self.m)
        self.base_without_penalty_bound = abs(self.marker_mu)+self.t+49
        self.repair_norm_bound = (self.base_without_penalty_bound
                                  +epsilon_max*self.full_perturbation_bound)
        self.Ksel = max(self.Ksel, 8*self.repair_norm_bound)
        self.selector_defect_gap = self.Ksel-4*self.repair_norm_bound


if __name__ == '__main__':
    import json
    H = CertifiedPeriodicBond()
    print(json.dumps(dict(q=H.q, dimension=H.d,
                          original_penalty=H.original_selector_penalty,
                          strengthened_penalty=H.Ksel,
                          norm_without_penalty_bound=H.repair_norm_bound,
                          defect_gap_bound=H.selector_defect_gap), indent=2))
