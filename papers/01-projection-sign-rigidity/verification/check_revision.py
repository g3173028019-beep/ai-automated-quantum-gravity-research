"""Exact revision check. Python >=3.10; standard library only.
Run from any working directory: python verification/check_revision.py
Writes verification/revision-checks.json. No network or historical project needed.
"""
from pathlib import Path
from fractions import Fraction
from itertools import combinations
import json
N=6
edges=list(combinations(range(N),2))
free=[e for e in edges if 0 not in e]
triangles=list(combinations(range(N),3))
def field(mask):
    signs={e:1 for e in edges}
    for a,e in enumerate(free):signs[e]=1 if mask>>a&1 else -1
    def get(i,j):return signs[tuple(sorted((i,j)))]
    h={t:get(t[0],t[1])*get(t[1],t[2])*get(t[0],t[2]) for t in triangles}
    rows=tuple(sum(h[tuple(sorted((i,j,k)))] for k in range(N) if k not in (i,j)) for i,j in edges)
    return tuple(h[t] for t in triangles),rows
K=[[Fraction(1,6) for _ in range(N)] for _ in range(N)]
assert all(sum(K[i][k]*K[k][j] for k in range(N))==K[i][j] for i in range(N) for j in range(N))
h1,r1=field(58);h2,r2=field(86)
assert h1!=h2 and r1==r2
expected=(0,0,0,-2,-2,0,0,2,2,0,2,2,-2,-2,0)
assert r1==expected and r1!=(4,)*15
count=sum(field(mask)[1]==(4,)*15 for mask in range(1<<len(free)))
assert count==1
out={"passed":True,"arithmetic":"exact integers and fractions","projection":"ones(6)/6","rank":1,"full_support":True,"gauge_fixed_assignments":1024,"target_fiber_count":count,"masks":[58,86],"common_Ah_times_36":r1,"h1":h1,"h2":h2,"triangles":triangles,"all_edges":edges,"free_edges":free,"scope":"Ancillary full-support fixed-fiber characterization; rank one is outside the generic interior-rank theorem. The all-size generic theorem rests on its analytic proof, not this finite check."}
p=Path(__file__).resolve().with_name('revision-checks.json');p.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:out[k] for k in ['passed','gauge_fixed_assignments','target_fiber_count','rank']}))
