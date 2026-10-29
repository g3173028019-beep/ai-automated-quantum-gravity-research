"""Finite checks of exact Abel block formulas; no sampling or scientific packages."""
import json
import math
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "reproduced-output"
OUT.mkdir(exist_ok=True)

def exp_minus_exp(log_z):
    # exp(-z), avoiding huge intermediate z.
    return 0.0 if log_z > math.log(746.0) else math.exp(-math.exp(log_z))

def block_probability(log_epsilon, k):
    log_start = (2**k) * math.log(2.0)
    log_end = (2**(k + 1)) * math.log(2.0)
    return exp_minus_exp(log_epsilon + log_start) - exp_minus_exp(log_epsilon + log_end)

rows = []
for k in range(1, 9):
    log_epsilon = -0.5 * ((2**k) + (2**(k + 1))) * math.log(2.0)
    epsilon = math.exp(log_epsilon)
    mass = block_probability(log_epsilon, k)
    u = math.exp(-0.5 * (2**k) * math.log(2.0))
    analytic_error = -math.expm1(-u) + math.exp(-1/u)
    value = sum((j % 2) * block_probability(log_epsilon, j) for j in range(16))
    omitted_tail = exp_minus_exp(log_epsilon + (2**16) * math.log(2.0))
    assert abs(value - (k % 2)) <= analytic_error + omitted_tail + 2e-14
    assert omitted_tail == 0.0
    normalized_denominator = epsilon / (-math.expm1(-epsilon)) * (1 + value)
    rows.append(dict(k=k, epsilon=epsilon, block_value=k % 2,
                     block_probability=mass, error_bound=analytic_error,
                     abel_witness=value,
                     epsilon_Z_for_coefficients_1_plus_x=normalized_denominator))

parity = [dict(epsilon=e, even_probability=1/(1+math.exp(-e)),
               shift_functional_norm=-2*math.expm1(-e))
          for e in [1.0, 0.1, 0.01, 0.001]]

# The fixed-word/microlocal ray witnesses failure of total-ideal membership.
fixed_word_ray = [dict(gamma=0, k=k, total_complexity=k, lift=1)
                  for k in [1, 10, 100, 1000]]
result = dict(status="passed finite exact-formula checks", block_rows=rows,
              parity_rows=parity, word_ideal_counterexample=fixed_word_ray,
              proof_boundary="Finite rows illustrate analytic bounds; no finite test proves a weak-star limit or physical admissibility. Block mass and witness round to endpoints in late rows; the separately stable analytic error bound stays positive.")
(OUT / "calculation-output.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
