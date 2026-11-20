"""Admit an actual fixed-geometry P255 state, without declaring convergence."""
import json
from screen_io import require,manifest


def load_state(run,role,state,*,h_A,repeat):
    directory=run/'calculation'/f'{role}-{state}'
    path=directory/'state.json';record=json.loads(path.read_bytes())
    require(record['schema']=='pdd-static-quartet-state-v1','Wrong static state schema')
    require(record['role']==role and record['state']==state,'Wrong static state role')
    require(record['charge']==(1 if state=='probe' else 0),'Wrong static-state charge')
    require(manifest(directory/'state.gpw')==record['gpw'],'Changed actual static GPW')
    spec=record['specification']
    require(spec['h_A']==h_A and spec['start_R_A']==.1,'Wrong approved h/radius control')
    require(spec['repeat']==repeat,'Wrong approved static control cell')
    return record,manifest(path)
