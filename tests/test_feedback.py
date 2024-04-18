from unittest.mock import Mock

from brdr.auto_referencer import AutoReferencer


def test_feedback():
    feedback_mock = Mock()
    AutoReferencer(feedback=feedback_mock)
    # assert pushinfo called with parameter "AutoReferencer initialized"
    feedback_mock.pushInfo.assert_called_with("AutoReferencer initialized")
