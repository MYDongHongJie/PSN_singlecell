#!/usr/bin/env python
"""Tests for `taskserver` package."""

import pytest

from taskserver import core
from fastapi.testclient import TestClient

client = TestClient(core.app)


@pytest.fixture
def response():
    """Sample pytest fixture.
    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_add_task():
    data = {
        "cmd": "test_cmd",
        "workroot": "test_workroot",
        "projectId": "test_projectId",
        "taskId": "test_taskId",
        "priority": -2,
        "queue": "test_big",
        "slot": 4
    }
    response = client.post('/tasks/add', json=data)
    assert response.status_code == 200


def test_remove():
    data = {'qsubId': "testid"}
    response = client.post('/tasks/remove', json=data)
    assert response.status_code == 200
