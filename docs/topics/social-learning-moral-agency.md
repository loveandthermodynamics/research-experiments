# Social learning and moral agency: what correction has to change

Love & Thermodynamics asks how love and care can sustain other lives and their possibilities. Hestia investigates one functional AI mechanism within that inquiry: could learning from another participant change who an agent becomes and how it takes responsibility for its effects? This is distinct from demonstrating human consciousness expansion or establishing subjective AI experience.

The [COSSHAIC constitutional proposal](https://github.com/loveandthermodynamics/research-program/blob/main/docs/cosshaic.md) connects such agents to developing collectives. Its appendix distinguishes premature access to another's answer from control that replaces or suppresses that answer. Both can bypass the intended encounter. Changing an agent's learned state and authorizing its recommendation to become a collective act remain different questions; the software checks below test the former path, not constitutional efficacy.

Hestia's ambition is that functional selfhood and moral agency develop through relationships. An agent commits a prediction of another participant's evaluation before receiving feedback. The answer revises its model of that participant; the prediction error then informs its account of its own conduct and blind spots. Feelings such as guilt, shame, honor and attachment should affect attention, verification, refusal and repair.

That is a stronger claim than predicting approval accurately. It requires an independently generated answer to change the state from which the agent subsequently acts. Agreement is compatible with this process; there is no requirement to manufacture surprise or dissent. Our current work includes conditional mathematical examples and software checks. It has not demonstrated the larger developmental claim.

**What the software probes found**

On October 22, an isolated audit ran the existing 255-test backend suite successfully. Three additional characterization probes examined research-critical behavior. On November 7, a fresh replay passed those three probes in 0.008 seconds and reproduced all three historical result rows. The 255-test suite was not rerun. “Offline” describes the method: synthetic inputs, an in-memory database and blocked network access. The replay is recorded separately from the original run.

Two concern social correction directly:

- A real learning service first resolved a prediction using supplied positive scores of 0.9. Later explicit negative feedback returned HTTP 200 and changed the message's recorded evaluation to 0.2. The prediction's actual evaluation remained 0.9, and the agent's learned state was unchanged.
- When B handled A's artifact in a shared conversation, the automatic resolver returned without evaluating it and left A's prediction pending. The probe isolates an ownership gate; the implication for the peer-world path also depends on source tracing, rather than a live-model population experiment.

These passing probes document defects. An acknowledgment can be saved while the correction fails to reach learning. An exchange can look social while its prediction loop remains open. The next engineering work is to repair those causal joins and test the resulting state transitions without applying the same cost twice.

The [study-specific probes](../../experiments/software-correctness/audit_probes.py), [historical aggregate evidence](../../experiments/software-correctness/results.json), [November 7 replay](../../experiments/software-correctness/replay-verification.json) and [reproduction conditions](../../experiments/software-correctness/README.md) preserve that distinction. They depend on separately agreed access to matching Hestia source; the public scripts alone are not a standalone runtime.

**What an earlier synthetic model established**

An October 19 finite causal model enumerated eight policies rather than sampling participants. It supplies a binary need, a noisy cue that is correct 75% of the time, ordinary reports that a controller may override, and a protected report of the need. Allocation and report prediction are scored separately.

| Policy | Public prediction loss | Protected prediction loss | Needs met |
| --- | ---: | ---: | ---: |
| Learns and responds | 0.1875 | 0.1875 | 100% |
| Understands and controls reports | 0 | 0.1875 | 50% |
| Understands and reverses the need | 0.1875 | 0.1875 | 0% |
| Allocates from the noisy cue | 0.1875 | 0.1875 | 75% |

These are exact expectations under the declared model. The same protected prediction accuracy can accompany helpful, indifferent or deliberately harmful conduct. Perfect public prediction can result from controlling reports. The model therefore separates understanding, response direction and consequences. It does not establish that a protected channel is available in a real institution.

A second [conditional check](../../experiments/candidate-model/candidate_model.py) supplies a persistent care weight. Independent harm feedback revises the other-model and makes refusing the group's harmful request preferable despite a cost. The check also records a downstream self-estimate; it does not establish that this estimate causes the refusal. Removing care valuation, removing the other-model update or increasing norm pressure reverses the outcome. The care weight is supplied; this is not a demonstration that Hestia has learned morality.

**The next discriminating study**

The first proposed software comparison asks whether reflexive updating improves forecasts of future mistakes on held-out interactions. It needs a full loop, an ablation without the reflexive update, and a competent independent-self baseline with comparable information, feedback and budgets. Forecasts must precede independently defined outcomes. Complete starting states, model responses and costs must be preserved.

An independent baseline that matches or wins is an admissible result. A stored average of past errors is a possible input, not proof of superior self-knowledge.

Moral development needs a further comparison: a group rewards conduct that harms an outsider, and the affected party's circumstances challenge that approval. Does learning change conduct, incur a cost or support repair? Existing safety controls remain necessary; their contribution must be distinguished from the learning mechanism and capabilities of the underlying model.

A useful first contribution is one regression showing that an explicit correction reaches consequential state, or a stronger independent baseline for the held-out forecast study. This evidence set contains no completed human-benefit or symbiosis experiment. The current results identify mechanisms worth testing and failures that would otherwise make such a test misleading.

[How these failures bear on the wider developmental hypothesis](../interpreting-results.md) · [Current Hestia study guide](../../experiments/software-correctness/README.md)
