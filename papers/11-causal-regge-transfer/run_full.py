"""Optional complete-action replay; deliberately separate from the cheap default."""
from pathlib import Path
import argparse,json,sys,time
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'provenance/frozen_producer'))
import numpy as np
from scripts.multilevel_causal_spin2 import multilevel_spin_two_certificate
def encode(x):
    if isinstance(x,np.ndarray):return encode(x.tolist())
    if isinstance(x,complex):return {'real':x.real,'imag':x.imag}
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [encode(v) for v in x]
    if isinstance(x,np.generic):return encode(x.item())
    return x
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--level',type=int,choices=[1,5],required=True)
args=parser.parse_args()
start=time.perf_counter(); result=multilevel_spin_two_certificate(args.level,1.2e-4)
result['elapsed_seconds']=time.perf_counter()-start
(HERE/'results').mkdir(exist_ok=True)
(HERE/f'results/full-level-{args.level}.json').write_text(json.dumps(encode(result),indent=2))
print(f'Completed optional level {args.level} in {result["elapsed_seconds"]:.2f} seconds.')
