from unittest.mock import Mock

from brdr.aligner import Aligner


def test_feedback():
    feedback_mock = Mock()
    Aligner(feedback=feedback_mock)
    # assert pushinfo called with parameter "Aligner initialized"
    feedback_mock.pushInfo.assert_called_with("Aligner initialized")
