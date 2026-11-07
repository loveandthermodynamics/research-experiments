# Reproducing the record

[Overview](README.md) · [Catalogue](CATALOGUE.md) · [Provenance and reuse](NOTICE.md)

There are four different operations here: checking the delivered bytes, rerunning a small calculation, replaying a saved event ledger, and regenerating a full ensemble. A successful operation supports its own scope. It does not establish that the model describes organisms, that the interpretation is novel, or that participants benefit.

To check delivered bytes against [release-integrity.json](release-integrity.json), run `python scripts/verify_archive.py` from the repository root. That manifest intentionally excludes its own bytes and Git internals. It is a corruption/change check, not an independent signature or a scientific acceptance test.

## Preserve the historical outputs

Many original scripts write result JSON, certificates or figures beside themselves. Use disposable copies. Do not run an output-producing script in the original research corpus or overwrite a frozen protocol in order to make a check pass. The quick replay helper below creates temporary copies automatically.

The original scientific Python files are retained with source hashes. Publication-only metadata adaptations are listed in [provenance.json](provenance.json) and the [physical snapshot ledger](provenance/physical-snapshot.json). Some original protocol or certificate bindings refer to the historical bytes of documents whose workstation commands or private metadata have been adapted for publication. Those historical hashes remain historical; they have not been silently rebound to the public copy. A full rerun that changes a protocol or its bindings is a derived run and needs its own record.

## Two dependency-free model replays

With Python 3.11 or later:

```sh
python scripts/replay_quick.py
```

Or use the script's declared Python environment through uv:

```sh
uv run --no-project scripts/replay_quick.py
```

This runs the early candidate model and multiplicative-persistence model in temporary copies and compares every saved output field, allowing declared numeric roundoff. It checks the ecological, policy, conscience, controller and scaling examples as well as the persistence/repair/common-shock counterexamples. It does not run Hestia or contact a model provider.

The [November 7 check](verification/quick-replay-2026-11-07.json) matched 201 and 80 leaf values respectively, with relative tolerance `1e-8` and absolute tolerance `1e-10`, using Python 3.12.11. A future run can save a new report without changing the historical results:

```sh
python scripts/replay_quick.py --output .replay/my-quick-replay.json
```

## Physical model dependencies

The recorded reference environment is Python 3.12.4, NumPy 2.2.6, SciPy 1.16.2, Matplotlib 3.9.2, mpmath 1.3.0 and SymPy 1.14.0. These identify the reference software context; they do not guarantee bitwise-identical trajectories across platforms or library builds. Preserve the versions actually used for a new check and compare the relevant physical observations, not only elapsed times or gzip container hashes.

A uv-managed scientific environment can be requested explicitly:

```sh
uv run --no-project --python 3.12 \
  --with numpy==2.2.6 --with scipy==1.16.2 --with matplotlib==3.9.2 \
  --with mpmath==1.3.0 --with sympy==1.14.0 \
  python path/to/the/check.py
```

Run from the copied check's documented working directory. Preserve the sibling layout beneath `experiments/physical-individuation/`; several phases import an earlier local model or retained dependency. Do not flatten the directories or combine same-named model/verifier files.

## Small numerical, symbolic and event-ledger checks

The archive includes the early lattice and wound-intervention NPZ ledgers. They support exact event-by-event checks of those saved synthetic histories without delivering every later ensemble. In a disposable copy of `experiments/physical-individuation/`, the original entry points include:

```sh
python verify_lattice_replay.py
python verify_integrity.py
python verify_closure.py
python verify_tail_bridge.py
```

The phase directories also retain small interval, path-ranking, generalist and capacity certificates. Some write a new certificate or verification report. Read the script and its model/claim note before running it; a certificate's assumptions matter as much as its successful exit code.

Fresh local checks performed while assembling this release are reported in [the physical replay summary](verification/physical-small-replays-2026-11-07.json). Its per-check status is authoritative: merely including a verifier in this repository does not mean it was run during this release.

## Full stochastic ensembles

The large pilot/confirmation/preparation/continuation histories, per-record replay caches and duplicate ZIP bundles are omitted. The archive retains source, protocols, declared readouts, compact aggregate findings, historical verification records and selected small fixtures. It therefore supports inspection and future reconstruction, while the omitted complete ensembles cannot be independently replayed event-by-event from this checkout alone.

For a full historical replay, obtain the corresponding original raw histories and original binding material, or regenerate the declared study in a separate working copy with an explicit provenance record. The [catalogue](CATALOGUE.md) identifies the study, preparation, sample structure, controls and evidence. The relevant model, study design, protocol and checkpoint describe each phase's commands and gates.

Keep these distinctions when regenerating:

- Pilots, confirmations and later replications are different datasets. Do not reselect parameters after seeing the confirmation result.
- A prepared assembly's maintenance or recovery is separate from formation from dispersed components.
- Shared preparations, paired branches and reused prefixes are not independent new replicates.
- A changed file digest caused by elapsed-time or compression metadata is different from a changed physical history. Conversely, a matching seed does not excuse a changed generator or protocol.
- Report a derived public protocol as a new version. Do not replace original evidence to conceal a mismatch.

No full later ensemble was rerun for this release. Its historical findings remain historical, including failed primary gates and inconclusive comparisons.

## Hestia software probes

The [software-correctness package](experiments/software-correctness/README.md) contains two compact audit scripts, historical synthetic outputs, source/evidence hashes and a [November 7 replay](experiments/software-correctness/replay-verification.json). The three probes were rerun with an in-memory test environment and blocked sockets/DNS. Their matching outcomes reproduce defects; they do not pass a scientific acceptance criterion.

These scripts require separately agreed access to a matching private Hestia implementation. The public archive alone cannot run them. The historical 255-test suite was not rerun while preparing this release. No held-out full-loop/ablation/independent-self comparison, human-benefit study or collective-development trial is included.

## Reporting a replay or correction

Name the record ID, source revision and actual environment. Link the command, input/protocol, output and what matched or differed. State whether the check was a mathematical calculation, numerical recomputation, event replay, conventional software test or a new scientific comparison. Preserve failed runs and distinguish an internal review from external scientific validation.

Use the [research forum](https://love-and-thermodynamics.discourse.group/) for interpretation and the [program tracker](https://github.com/loveandthermodynamics/research-program/issues) for an agreed bounded contribution. Read [NOTICE.md](NOTICE.md) and agree any necessary access and reuse terms before substantial work.
