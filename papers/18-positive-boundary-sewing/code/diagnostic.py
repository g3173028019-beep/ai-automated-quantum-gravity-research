"""Independent small certificates for scalar versus amplified positivity."""
import numpy as np

X=np.array([[0,1],[1,0]],complex)
units=[]
for i in range(2):
    row=[]
    for j in range(2):
        e=np.zeros((2,2),complex);e[i,j]=1;row.append(e)
    units.append(row)
choi=sum(np.kron(units[i][j],units[i][j].T) for i in range(2) for j in range(2))
assert np.allclose(np.linalg.eigvalsh(choi),[-1,1,1,1])
print('transpose_choi_spectrum=-1,1,1,1')
print('transpose_ancillary_control=PASS')
plus=np.array([[.5,.5],[.5,.5]])
minus=np.array([[.5,-.5],[-.5,.5]])
assert np.allclose(np.diag(plus),np.diag(minus))
assert np.isclose(np.trace((plus-minus)@X),2)
print('observer_restriction_witness=PASS')
comp=lambda a:a[:1,:1]
defect=np.linalg.norm(comp(units[0][1]@units[1][0])-comp(units[0][1])@comp(units[1][0]))
assert np.isclose(defect,1)
print('full_matrix_compression_defect=1')
print('compression_invalid_control=PASS')
d=3;z=np.exp(2j*np.pi/d);inv2=pow(2,-1,d)
shift=np.roll(np.eye(d),1,axis=0);clock=np.diag(z**np.arange(d))
res=max(np.linalg.norm(z**((a*b*inv2)%d)*np.linalg.matrix_power(shift,a)@np.linalg.matrix_power(clock,b)
    -z**((-b*a*inv2)%d)*np.linalg.matrix_power(clock,b)@np.linalg.matrix_power(shift,a))
    for a in range(d) for b in range(d))
assert res<1e-12
print('weyl_dictionary_residual='+str(res))
print('weyl_dictionary=PASS')
print('overall=PASS')
