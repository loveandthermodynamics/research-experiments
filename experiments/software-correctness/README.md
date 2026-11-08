# Does feedback change what Hestia learns?

For another participant's answer to matter, it has to change the state from which an agent later acts. Hestia proposes that learning about others also changes its account of itself, with costs and feelings shaping conduct. These checks ask whether the reviewed software carries out some necessary steps of that proposal.

The October 22, 2026 probes exercise real Hestia services with synthetic scores and agents. All three are **characterization tests**: they assert the observed defective behavior. A pass therefore reproduces a defect, rather than showing that the mechanism works as intended. The November 7 replay reproduced the three saved outputs; it did not rerun the separate historical 255-test suite.

## Findings and the claims they constrain

The saved findings are:

1. After inferred feedback has resolved a prediction at `0.9`, explicit negative feedback changes the message's recorded evaluation to `0.2`, while the resolved prediction and learned state retain the earlier outcome. This checks the view and database; no browser display was tested.
2. Eight epochs containing only `mean_fitness = 100 / epoch` produce a LIFE-L6 compliance label despite containing no group-size or resource-cost observations. That label is not evidence of negentropic scaling.
3. When agent B handles the next turn, agent A's pending peer prediction remains unresolved at the tested resolver gate.

The first and third results show failures in the tested learning path. The second shows that the scaling acceptance label is not a valid measurement of the proposed resource–group-size relation. These are substantive problems to repair before using this implementation to test developmental claims. They do not measure human benefit, emergent collective identity or reciprocal development, and do not constitute a comparison of the full reflexive loop against another learning architecture.

The wider ambition concerns what new capacities agents and people could develop together, while remaining answerable for care and costs. See [what the negative results test](../../docs/interpreting-results.md), including why the historical fictional adversaries are conceptual cases rather than results of these software probes.

[results.json](results.json) preserves the three synthetic outputs and the separate historical 255-test suite summary. [manifest.json](manifest.json) identifies the reviewed working-tree bytes as well as the repository's then-current commit. A commit ID alone does not identify a dirty working tree. The initial blank-key test failure and later successful placeholder-key suite are separate historical attempts; the raw logs remain outside this public snapshot.

## Method: what “offline” means

“Offline” describes the execution conditions, not the subject, maturity or availability of the project. The probes use an isolated in-memory test database and supplied feedback. Network and DNS access are blocked, so no model provider or human participant generates the answers. They check specified service and state transitions in the matching runtime, not a live conversation or population experiment.

## Reproduction requires the Hestia runtime

The two scripts are not a standalone implementation. They import Django, DRF and private Hestia modules. Arrange access to the matching Hestia source and dependencies separately, and compare its file hashes with the manifest before interpreting a replay.

From an authorized checkout's `api/` directory, with the archive located elsewhere:

```sh
export HESTIA_AUDIT_DIR=/path/to/research-experiments/experiments/software-correctness
PYTHONPATH="$HESTIA_AUDIT_DIR" uv run --frozen --offline python "$HESTIA_AUDIT_DIR/run_offline_tests.py" audit_probes
```

The runner forces in-memory SQLite, channels and Celery configuration, supplies an explicitly noncredential constructor placeholder, and denies socket connections and DNS. It does not need a provider key, production database or live model call. Keep any new local stack traces private; report the `AUDIT_PROBE` JSON rows and the exact source basis.

The proposed next scientific study compares the full reflexive loop, an ablation and a competent independent-self baseline under matched feedback and resources. These characterization probes do not run that comparison. See the [social-learning digest](../../docs/topics/social-learning-moral-agency.md) for its role in the wider program.

Read the repository's [reuse notice](../../NOTICE.md) before treating public visibility as a software license.
