import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import dataclass, field
import typing

from .mechanistic_params import PARAMS, INITIAL_HORMONES, PARA_ODE

@dataclass
class Follicle:
    size: float
    fsh_sensitivity: float
    birth_time: float
    destiny: int = -1  # -1: active, -2: atretic, 1: ovulated, 3: waiting, 4: ready to ovulate
    time_decrease: float = 0.0

class MenstrualMechanisticModel:
    def __init__(self, seed: int = 42):
        self.rng = np.random.default_rng(seed)
        self.params = PARAMS.copy()
        self.para_ode = PARA_ODE.copy()
        self.hormones = INITIAL_HORMONES.copy()
        self.follicles: list[Follicle] = []
        self.t = 0.0
        self.last_ovulation_t = -14.0
        
        # Stress factor (0 to 1, where 1 is total suppression of GnRH)
        self.stress_input = 0.0

    def get_active_follicle_indices(self):
        return [i for i, f in enumerate(self.follicles) if f.destiny in [-1, 3, 4]]

    def compute_targets(self, t, foll_sizes):
        if len(foll_sizes) > 0:
            sf = np.pi * np.sum((foll_sizes**5 / (foll_sizes**5 + 15**5)) * (foll_sizes**2))
        else:
            sf = 0.0
            
        e2_target = 13.0 * (1.0 + 0.02 * sf) + 150.0 * np.exp(-0.06 * (t - (self.last_ovulation_t + 7.0))**2)
        p4_target = 15.0 * np.exp(-0.06 * (t - (self.last_ovulation_t + 7.0))**2)
        
        return e2_target, p4_target

    def ode_system(self, t, y, active_indices):
        num_foll = len(active_indices)
        foll_sizes = y[:num_foll]
        
        dy = np.zeros_like(y)
        
        I_E2 = -53
        I_P4 = -51
        I_LH = -49
        I_FSH = -47
        I_RP_LH = -50
        I_RP_FSH = -48
        I_FSHFOLL = -46
        I_RLH = -19
        I_RFSH = -15
        I_RFSH_DES = -13
        I_FSHR = -14
        I_GNRH = -5
        I_REC_GA = -4
        I_REC_GI = -3
        I_GREC_A = -2
        I_GREC_I = -1
        
        y_e2 = max(0, y[I_E2])
        y_p4 = max(0, y[I_P4])
        y_lh_val = max(0, y[I_LH])
        y_fsh_val = max(0, y[I_FSH])
        
        e2_target, p4_target = self.compute_targets(t, foll_sizes)
        dy[I_E2] = (e2_target - y_e2) * 50.0 
        dy[I_P4] = (p4_target - y_p4) * 50.0
        
        stress_suppression = 1.0 - self.stress_input
        yGfreq = (16.0 / (1.0 + (y_p4 / self.params[203])**self.params[204])) * \
                 (1.0 + y_e2**self.params[206] / (self.params[205]**self.params[206] + y_e2**self.params[206]))
        yGfreq *= max(0, stress_suppression)
        
        # numerically stable Hill functions
        def hill_pos(s, k, n):
            if s <= 0: return 0.0
            ratio = s / k
            if ratio > 10.0: return 1.0 # saturation
            return ratio**n / (1.0 + ratio**n)

        def hill_neg(s, k, n):
            if s <= 0: return 1.0
            ratio = s / k
            if ratio > 2.0: return 0.0 # sharp cutoff for high n
            return 1.0 / (1.0 + ratio**n)

        yGmass = hill_pos(y_e2, self.params[208], self.params[209]) + \
                 hill_neg(y_e2, self.params[210], self.params[211])
        
        hp_e2 = hill_pos(y_e2, self.params[3], self.params[6])
        f_LH_prod = (self.params[1] + self.params[2] * hp_e2) / (1.0 + (y_p4 / self.params[4])**self.params[7])
        f_LH_rel = self.params[16] + self.params[5] * hill_pos(y[I_GREC_A], self.params[8], self.params[9])
        
        dy[I_RP_LH] = f_LH_prod - f_LH_rel * y[I_RP_LH]
        dy[I_LH] = (1.0 / self.params[12]) * f_LH_rel * y[I_RP_LH] - self.params[230] * y_lh_val * y[I_RLH] - self.params[231] * y_lh_val
        
        hm_freq = 1.0 / (1.0 + (yGfreq / self.params[11])**self.params[13])
        f_FSH_prod = (self.params[21] / (1.0 + (y_p4 / 5.0)**2)) * hm_freq
        f_FSH_rel = self.params[17] + self.params[28] * hill_pos(y[I_GREC_A], self.params[18], self.params[20])
        
        dy[I_RP_FSH] = f_FSH_prod - f_FSH_rel * y[I_RP_FSH]
        dy[I_FSH] = (1.0 / self.params[12]) * f_FSH_rel * y[I_RP_FSH] - (self.params[243] + self.params[241]) * y_fsh_val
        
        dy[I_FSHFOLL] = self.params[243] * y_fsh_val * self.params[246] - self.params[240] * y[I_FSHFOLL] * y[I_RFSH]
        dy[I_RFSH] = self.params[242] * y[I_RFSH_DES] - self.params[240] * y[I_FSHFOLL] * y[I_RFSH]
        dy[I_FSHR] = self.params[240] * y[I_FSHFOLL] * y[I_RFSH] - self.params[244] * y[I_FSHR]
        dy[I_RFSH_DES] = self.params[244] * y[I_FSHR] - self.params[242] * y[I_RFSH_DES]
        
        dy[I_GNRH] = self.params[301] * yGmass * yGfreq - self.params[302] * y[I_GNRH] * y[I_REC_GA] + \
                     self.params[303] * y[I_GREC_A] - self.params[300] * y[I_GNRH]
        dy[I_REC_GA] = self.params[303] * y[I_GREC_A] - self.params[302] * y[I_GNRH] * y[I_REC_GA] - \
                       self.params[306] * y[I_REC_GA] + self.params[307] * y[I_REC_GI]
        dy[I_REC_GI] = self.params[311] + self.params[306] * y[I_REC_GA] - self.params[307] * y[I_REC_GI] + \
                       self.params[305] * y[I_GREC_I] - self.params[308] * y[I_REC_GI]
        dy[I_GREC_A] = self.params[302] * y[I_GNRH] * y[I_REC_GA] - self.params[303] * y[I_GREC_A] - \
                       self.params[309] * y[I_GREC_A] + self.params[310] * y[I_GREC_I]
        dy[I_GREC_I] = self.params[309] * y[I_GREC_A] - self.params[310] * y[I_GREC_I] - \
                       self.params[305] * y[I_GREC_I] - self.params[304] * y[I_GREC_I]
        
        sum_v = np.sum(foll_sizes**self.para_ode['v'])
        fsh_rez_comp = y[I_FSHR]
        
        for idx in range(num_foll):
            foll_obj = self.follicles[active_indices[idx]]
            f_fsh = (fsh_rez_comp)**1 / (fsh_rez_comp**1 + foll_obj.fsh_sensitivity**1)
            gamma_p4 = self.para_ode['gamma'] * ((1.0 / (1.0 + (y_p4 / 3.0)**5)) + (fsh_rez_comp**10 / (0.59**10 + fsh_rez_comp**10)))
            kappa_fsh = self.para_ode['k'] * (0.85**25 / (0.85**25 + fsh_rez_comp**25))
            
            x_growth = f_fsh * (self.para_ode['xi'] - y[idx]) * y[idx] * \
                       (gamma_p4 - (kappa_fsh * (sum_v - (self.para_ode['mu'] * y[idx]**self.para_ode['v']))))
            
            if foll_obj.destiny == -2:
                dy[idx] = -0.05 * y[idx] * (t - foll_obj.time_decrease)
            else:
                dy[idx] = x_growth
                
        return dy

    def step(self, dt=0.25):
        # 1. Recruit
        y_fsh_val = self.hormones[-47]
        fsh_imp = y_fsh_val**3 / (y_fsh_val**3 + 12.35**3)
        lambda_eff = (20.0/14.0) + 4.0 * (20.0/14.0) * fsh_imp
        
        n_new = self.rng.poisson(lambda_eff * dt)
        for _ in range(n_new):
            sens = self.rng.normal(0.6, 0.55 * 0.6)
            self.follicles.append(Follicle(size=4.0, fsh_sensitivity=max(0.1, sens), birth_time=self.t))
            
        # 2. Setup Integration
        active_indices = self.get_active_follicle_indices()
        num_foll = len(active_indices)
        
        y0 = np.zeros(num_foll + 54)
        for i, idx in enumerate(active_indices):
            y0[i] = self.follicles[idx].size
        y0[num_foll:] = self.hormones
        
        # Event function for ovulation
        def ovulation_event(t, y):
            num_f = len(active_indices)
            lh_val = y[num_f - 49]
            if lh_val < 25.0: return 1.0
            # Check if any follicle > folmax
            f_sizes = y[:num_f]
            if len(f_sizes) == 0: return 1.0
            max_size = np.max(f_sizes)
            return self.para_ode['folmax'] - max_size
            
        ovulation_event.terminal = True
        ovulation_event.direction = -1

        sol = solve_ivp(
            fun=lambda t, y: self.ode_system(t, y, active_indices),
            t_span=[self.t, self.t + dt],
            y0=y0,
            events=ovulation_event,
            method='BDF',
            rtol=1e-3, atol=1e-6
        )
        
        # 3. Update state
        y_final = sol.y[:, -1]
        for i, idx in enumerate(active_indices):
            self.follicles[idx].size = y_final[i]
        self.hormones = y_final[num_foll:]
        self.t = float(sol.t[-1]) if hasattr(sol.t, '__len__') else float(sol.t)
        
        ovulated = False
        if sol.status == 1: # Event triggered
            # Mark the largest follicle as ovulated
            sizes = y_final[:num_foll]
            if len(sizes) > 0:
                best_i = np.argmax(sizes)
                self.follicles[active_indices[best_i]].destiny = 1
                self.last_ovulation_t = self.t
                ovulated = True
        
        # 4. Atresia and other cleanup
        for idx in self.get_active_follicle_indices():
            f = self.follicles[idx]
            if f.size < 2.0:
                f.destiny = -2
                f.time_decrease = self.t
                
        return ovulated

def run_mechanistic_trial(days=100, stress_timeline=None):
    model = MenstrualMechanisticModel()
    ovulation_times = []
    
    while model.t < days:
        day_idx = min(int(model.t), days - 1)
        if stress_timeline is not None:
            model.stress_input = stress_timeline[day_idx]
        
        if model.step(0.25):
            ovulation_times.append(model.t)
                
    return ovulation_times
