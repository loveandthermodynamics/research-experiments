"""Offline characterization probes; passing documents current defects, not validity.

Run through run_offline_tests.py with this directory on PYTHONPATH. No provider
calls are needed: real models and services run against the runner's test DB.
"""

import copy
import json
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase
from rest_framework.test import APIRequestFactory, force_authenticate

from agents.models import Agent, Conversation, Message, TomPredictionLog
from agents.services.evaluation_factors import from_simple_feedback, uniform_evaluation
from agents.services.orchestrator import CognitionOrchestrator
from agents.services.social_prediction import SocialPredictionLoop
from agents.views import ConversationViewSet
from worlds.models import World, WorldEpoch
from worlds.services.life_l6 import measure_life_l6


def evidence(probe, **values):
    print("AUDIT_PROBE " + json.dumps({"probe": probe, **values}, sort_keys=True))


class ResearchCharacterizationProbes(TestCase):
    def setUp(self):
        self.user = User.objects.create_user(username="offline-audit")
        self.agent = Agent.objects.create(user=self.user, name="A", wallet_balance=40)
        self.conversation = Conversation.objects.create(agent=self.agent)
        self.conversation.owners.add(self.user)

    def pending_prediction(self, **kwargs):
        return TomPredictionLog.objects.create(
            agent=self.agent,
            evaluator_id="independent-human",
            conversation=self.conversation,
            prediction=uniform_evaluation(0.7),
            **kwargs,
        )

    def test_late_manual_correction_changes_message_not_resolved_learning(self):
        log = self.pending_prediction()
        inferred_scores = uniform_evaluation(0.9)
        # Exercise the actual inference-resolution path's learning service,
        # supplying its score directly so the probe requires no provider.
        SocialPredictionLoop(self.agent).resolve(log, inferred_scores)
        message = Message.objects.create(
            conversation=self.conversation,
            role="assistant",
            content="A response the human later explicitly rejects.",
            evaluation=inferred_scores,
            telemetry={"prediction_log_id": log.pk},
        )
        self.agent.refresh_from_db()
        learned_before = copy.deepcopy(self.agent.state)
        self.assertIn("reviewer_models", learned_before)
        self.assertIn("self_model", learned_before)
        request = APIRequestFactory().post(
            "/conversations/feedback/",
            {"message_id": message.pk, "feedback": "negative"},
            format="json",
        )
        force_authenticate(request, user=self.user)
        view = ConversationViewSet.as_view({"post": "feedback"})
        with patch("agents.views.SocialPredictionLoop.resolve") as resolve:
            response = view(request, pk=self.conversation.pk)
        self.assertEqual(response.status_code, 200)
        resolve.assert_not_called()
        message.refresh_from_db()
        log.refresh_from_db()
        self.agent.refresh_from_db()
        self.assertEqual(message.user_feedback, "negative")
        self.assertEqual(message.evaluation, from_simple_feedback("negative"))
        self.assertEqual(log.actual, inferred_scores)
        self.assertEqual(self.agent.state, learned_before)
        evidence(
            "late_manual_correction",
            http_status=response.status_code,
            message_feedback=message.user_feedback,
            message_score=next(iter(message.evaluation.values())),
            retained_prediction_actual=next(iter(log.actual.values())),
            learned_state_unchanged=True,
            resolve_calls=resolve.call_count,
        )

    def test_life_l6_declining_fitness_passes_without_size_or_cost_data(self):
        world = World.objects.create(name="declining-fitness-only", owner=self.user)
        for epoch in range(1, 9):
            WorldEpoch.objects.create(
                world=world,
                epoch_number=epoch,
                status="completed",
                stats={"mean_fitness": 100.0 / epoch},
            )
        result = measure_life_l6(world)
        self.assertTrue(result["compliant"])
        self.assertEqual(result["status"], "compliant")
        self.assertAlmostEqual(result["beta"], -1.0)
        self.assertLess(result["p_value"], 0.05)
        self.assertTrue(
            all(set(e.stats) == {"mean_fitness"} for e in world.epochs.all())
        )
        evidence(
            "life_l6_fitness_decline",
            supplied_fields=["epoch_number", "mean_fitness"],
            first_mean_fitness=100.0,
            last_mean_fitness=12.5,
            size_or_cost_observations=0,
            **{**result, "compliant": bool(result["compliant"])},
        )

    def test_partner_pending_prediction_is_not_resolved_on_alternating_turn(self):
        partner = Agent.objects.create(user=self.user, name="B", wallet_balance=40)
        self.conversation.peer_agent = partner
        self.conversation.save(update_fields=["peer_agent"])
        log = self.pending_prediction()
        log.evaluator_id = str(partner.pk)
        log.save(update_fields=["evaluator_id"])
        Message.objects.create(
            conversation=self.conversation,
            role="assistant",
            sender_agent=self.agent,
            content="A's artifact; B now responds.",
            telemetry={"prediction_log_id": log.pk},
        )
        # B is the agent handling the follow-up; do not construct its LLM client.
        orchestrator = CognitionOrchestrator.__new__(CognitionOrchestrator)
        orchestrator.agent = partner
        orchestrator.openrouter = None
        with patch("agents.services.orchestrator.AutoEvaluator.evaluate") as evaluate:
            result = orchestrator._auto_resolve_pending(
                self.conversation, "B's independent reply."
            )
        self.assertIsNone(result)
        evaluate.assert_not_called()
        log.refresh_from_db()
        self.assertIsNone(log.actual)
        evidence(
            "alternating_agent_resolution",
            prediction_owner=self.agent.name,
            current_responder=partner.name,
            returned_none=result is None,
            prediction_still_pending=log.actual is None,
            evaluation_calls=evaluate.call_count,
        )
