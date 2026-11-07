"""Run the existing Django suite with isolated storage and no network access."""

import os
import socket
import sys

os.environ.update(
    {
        "DATABASE_URL": "sqlite:///:memory:",
        # Constructor-only placeholder. Socket guards below prevent real calls.
        "OPENROUTER_API_KEY": "offline-audit-placeholder-not-a-real-key",
        "HESTIA_MODEL_ENV": "test",
        "CHANNEL_LAYER_BACKEND": "channels.layers.InMemoryChannelLayer",
        "REDIS_URL": "",
        "CELERY_BROKER_URL": "memory://",
        "CELERY_RESULT_BACKEND": "cache+memory://",
        "DJANGO_SETTINGS_MODULE": "hestia.settings",
    }
)


def deny_network(*args, **kwargs):
    raise RuntimeError("Network disabled for the isolated research audit")


socket.socket.connect = deny_network
socket.socket.connect_ex = deny_network
socket.create_connection = deny_network
socket.getaddrinfo = deny_network

sys.path.insert(0, os.getcwd())

from django.core.management import execute_from_command_line

execute_from_command_line(
    ["manage.py", "test", "--noinput", "--verbosity=1", *sys.argv[1:]]
)
