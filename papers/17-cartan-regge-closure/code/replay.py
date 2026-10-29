"""Replay existing producers read-only, persisting main-claim evidence here."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
import verify_first_order_angle_exact as exact
import sourced_regge_dynamics as dyn
import verify_cartan_regge_ranks as ranks
import sympy as sp

inputs=json.loads((Path(__file__).resolve().parents[1]/'inputs'/'calibration.json').read_text(encoding='utf-8'))
expected_points=[sp.Matrix([sp.Rational(x) for x in row]) for row in inputs['coordinates']]
assert expected_points==exact.POINTS
source=float(sp.Rational(inputs['source']))
rank_output={}
for side in ('initial','final'):
    matching,boundary=ranks.exact_matrices(side)
    combined=matching.col_join(boundary)
    mr=ranks.DomainMatrix.from_Matrix(matching).rank()
    cr=ranks.DomainMatrix.from_Matrix(combined).rank()
    rank_output[side]=dict(matching_rank=mr,combined_rank=cr,nullity=combined.cols-cr)
assert rank_output['initial']==dict(matching_rank=25,combined_rank=39,nullity=1)
assert rank_output['final']==dict(matching_rank=6,combined_rank=20,nullity=0)
blocks=[]
for i in range(2,6):
    multiplier,det,contribution=exact.exact_simplex_block(i)
    assert det!=0
    blocks.append(dict(omitted=i,multiplier=str(multiplier),determinant=str(det),schur=str(contribution)))
h=sp.factor(sum(exact.exact_simplex_block(i)[2] for i in range(2,6)))
c,D,v,target=exact.source_target()
assert h==sp.Rational(-224,16875) and h/target==sp.Rational(-1,8)
lengths=dyn.squared_lengths(dyn.canonical_radon_configuration())
_,upper=dyn.valid_q_interval(lengths)
roots=dyn.find_sourced_saddles(lengths,source,q_bounds=(2.25,upper))
assert len(roots)==2
clean=[]
for root in roots:
    assert abs(root['equation_residual'])<1e-10
    assert abs(root['cartan_internal_curvature']-root['regge_internal_curvature'])<3e-9
    clean.append({k:root[k] for k in ('q_01','hessian','equation_residual','regge_internal_curvature','cartan_internal_curvature')})
    clean[-1]['deficits']={str(k):v for k,v in root['deficits'].items()}
out=dict(metric_ranks=rank_output,blocks=blocks,hessian=str(h),D=str(D),source_target=str(target),ratio=str(h/target),
         source=source,saddles=clean,fold=dyn.positive_source_fold(lengths),
         coverage='Existing exact Hessian producer replay and numerical nonlinear main-claim witnesses; not a fresh independent first-order derivation or a global fold proof.')
expected=json.loads((Path(__file__).resolve().parents[1]/'outputs'/'established_verification.json').read_text(encoding='utf-8'))
assert out['blocks']==expected['blocks']
assert out['hessian']==expected['hessian'] and out['ratio']==expected['ratio']
for actual,saved in zip(out['saddles'],expected['saddles']):
    assert abs(actual['q_01']-saved['q_01'])<2e-10
    assert abs(actual['hessian']-saved['hessian'])<1e-7
assert abs(out['fold']['source_critical']-expected['fold']['source_critical'])<1e-10
(Path(__file__).resolve().parents[1]/'outputs'/'verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))
