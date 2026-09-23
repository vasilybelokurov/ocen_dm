"""Original failure-isolation script, retained as the debugging record.

Run from the project root with PYTHONPATH=src in the project environment.
Uses the preserved smoke-test checkpoint under tests/fixtures/. For a shorter
replay of the saved orbit/potential, use bin/diagnose_compact_frequency.py.
"""
import json
from pathlib import Path
import numpy as np
from ocen_dm.kinematics.df_fit import model_config_from_dict, build_df_model
from ocen_dm.kinematics.compact_recovery import refined_config
from ocen_dm.kinematics.positive_df import agama_pc
ag = agama_pc()
original = ag.actions
original_mapper = ag.ActionMapper
out = Path('results/diagnostics/compact_frequency_failure_20260923')
out.mkdir(parents=True, exist_ok=True)
class Mapper:
    def __init__(self, pot):
        self.pot = pot
        self.mapper = original_mapper(pot)
    def __call__(self, aa):
        result = self.mapper(aa)
        bad = np.any(~np.isfinite(result), axis=1)
        if np.any(bad):
            self.pot.export(str(out/'potential.ini'))
            print('BAD MAPPER INPUT', np.asarray(aa)[bad].tolist(), flush=True)
            (out/'actions.json').write_text(json.dumps(np.asarray(aa)[bad].tolist()))
        return result
ag.ActionMapper = Mapper
def probe(potential, xv, *a, **kw):
    try:
        return original(potential, xv, *a, **kw)
    except RuntimeError as exc:
        potential.export(str(out/'potential.ini'))
        print('FAILED VECTOR', np.asarray(xv).shape, str(exc), flush=True)
        for index, row in enumerate(np.atleast_2d(xv)):
            try:
                original(potential, row, *a, **kw)
            except RuntimeError:
                print('FAILED ROW', index, row.tolist(), flush=True)
                (out/'orbit.json').write_text(json.dumps(dict(index=index, xv=row.tolist(), error=str(exc))))
                break
        raise
ag.actions = probe
payload = json.loads(Path('tests/fixtures/compact_frequency/recovery_checkpoint.json').read_text())
c = refined_config(model_config_from_dict(payload['config']))
build_df_model(c)
