"""High-precision nonlinear controls of THERMAL_HISTORY.md.

Standard library only. Checks four trajectories, including a two-slot
memory model and a thermal derivative too large for L+g<rho.
These controls do not establish any microscopic RG derivative bound.
"""

from decimal import Decimal, localcontext
from pathlib import Path
import json


def dec(x):
    return Decimal(str(x))


def run_case(depth, memory=False, large_thermal_derivative=False):
    with localcontext() as context:
        context.prec = 240
        two, epsilon, A = dec(2), dec("0.02"), dec("1.1")
        p, rho = 4, dec(8)
        temperature = two**(-depth)
        xline, yline, xthermal, ythermal = [], [], [], []
        g_current = dec(20) if large_thermal_derivative else dec("0.006")
        g_old = dec("0.003") if memory else dec(0)
        g_value = dec("0.03") + dec("0.1")*(g_current+g_old)
        a_current = dec("0.52") if memory else dec("0.6")
        a_old = dec("0.05") if memory else dec(0)
        hessian = dec("0.1") if memory else dec("0.2")
        L_eff = a_current+a_old/rho
        L_p = a_current+a_old/dec(16)
        G_delta = g_current+g_old/rho
        G_line = g_current+two*g_old
        B = g_value/(dec(16)-L_p)
        F = A*(B*hessian+G_line)
        D = F/(rho-L_eff)*(G_delta/(rho*(1-dec(1)/16))).exp()

        scalar_c_old = dec("0.1") if memory else dec(0)
        scalar_e_old = dec("0.04") if memory else dec(0)
        C_line = dec("1.1")+two*scalar_c_old
        C_delta = dec("1.1")+scalar_c_old/rho
        E_line = dec("0.42")+two*scalar_e_old
        E_delta = dec("0.42")+scalar_e_old/rho
        scalar_bound = D*(C_delta+E_delta)+A*B*dec("0.6")+A*E_line

        max_delta_ratio, max_scalar_ratio = dec(0), dec(0)
        max_thermal_ratio, max_line_ratio = dec(0), dec(0)
        x0, y0, xt, yt = epsilon, dec(0), epsilon, dec(0)

        for step in range(depth+1):
            scale = two**(-step)
            thermal_power = (temperature/scale)**p
            older = lambda values: values[-1] if values else dec(0)
            xo, yo, xto, yto = map(older, (xline, yline, xthermal, ythermal))
            defect_scale = epsilon*scale
            marked_scale = defect_scale*thermal_power
            delta = (xt-x0)-(yt-y0)
            ratio = abs(delta)/marked_scale
            line_ratio = abs(x0-y0)/defect_scale
            thermal_ratio = abs(yt-y0)/thermal_power
            assert ratio <= D*(1+dec("1e-60"))
            assert line_ratio <= A
            assert thermal_ratio <= B
            assert max(map(abs, (x0, y0, xt, yt, xo, yo, xto, yto))) <= dec("0.1")

            def phi(z, old):
                return z+dec("0.3")*z*z+scalar_c_old*old

            def psi(z, old):
                return thermal_power*(dec("0.4")*(z+dec("0.2")*z*z)+scalar_e_old*old)

            r0 = phi(x0, xo)-phi(y0, yo)
            rt = phi(xt, xto)+psi(xt, xto)-phi(yt, yto)-psi(yt, yto)
            scalar_ratio = abs(rt-r0)/marked_scale
            assert abs(r0)/defect_scale <= A*C_line
            assert scalar_ratio <= scalar_bound*(1+dec("1e-60"))
            max_delta_ratio = max(max_delta_ratio, ratio)
            max_scalar_ratio = max(max_scalar_ratio, scalar_ratio)
            max_line_ratio = max(max_line_ratio, line_ratio)
            max_thermal_ratio = max(max_thermal_ratio, thermal_ratio)

            if step == depth:
                break

            def line(z, old):
                if memory:
                    linear = dec("0.5")*z if step == 0 else dec("0.4")*z+dec("0.05")*old
                    return linear+dec("0.05")*z*z
                return dec("0.5")*z+dec("0.1")*z*z

            def thermal(z, old):
                return line(z, old)+thermal_power*(dec("0.03")+g_current*z+g_old*old)

            next_values = (line(x0, xo), line(y0, yo), thermal(xt, xto), thermal(yt, yto))
            xline.append(x0)
            yline.append(y0)
            xthermal.append(xt)
            ythermal.append(yt)
            x0, y0, xt, yt = next_values

        return dict(depth=depth, memory=memory,
                    large_thermal_derivative=large_thermal_derivative,
                    L_eff=float(L_eff), rho=float(rho),
                    thermal_derivative_coefficient=float(g_current+g_old),
                    four_trajectory_ratio_max=float(max_delta_ratio),
                    four_trajectory_bound=float(D),
                    scalar_thermal_ratio_max=float(max_scalar_ratio),
                    scalar_thermal_bound=float(scalar_bound),
                    line_defect_ratio_max=float(max_line_ratio),
                    reference_thermal_ratio_max=float(max_thermal_ratio),
                    reference_thermal_bound=float(B))


def main():
    cases = [run_case(depth, memory, large)
             for depth in (1, 5, 10, 20, 40, 80, 120)
             for memory, large in ((False, False), (True, False), (False, True))]
    result = dict(status="all high-precision recursion checks passed", cases=len(cases),
                  precision_digits=240,
                  scope="Abstract nonlinear maps; no microscopic RG hypothesis tested",
                  results=cases)
    Path(__file__).with_name("history_checks.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k: v for k, v in result.items() if k != "results"}, indent=2))


if __name__ == "__main__":
    main()
