from pathlib import Path
import subprocess,sys,json
root=Path(__file__).resolve().parent
results=[]
for file in ['doubled_instrument.py', 'weyl_charge_cut.py', 'diagnostic.py']:
    cp=subprocess.run([sys.executable,str(root/'code'/file)],capture_output=True,text=True)
    assert cp.returncode==0 and '=FAIL' not in cp.stdout and '=PASS' in cp.stdout, (file,cp.stderr,cp.stdout)
    (root/'outputs').mkdir(exist_ok=True)
    (root/'outputs'/(Path(file).stem+'.txt')).write_text(cp.stdout,encoding='utf-8')
    results.append({'file':file,'exit_code':cp.returncode,'pass_rows':cp.stdout.count('=PASS')})
(root/'outputs'/'receipt.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf-8')
print(json.dumps(results))
