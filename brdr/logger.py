import logging
import sys

LOGGER = logging.getLogger("BRDR")
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)
formatter = logging.Formatter(
    fmt="%(asctime)s - %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
)
handler.setFormatter(formatter)
LOGGER.addHandler(handler)


class Logger:
    def __init__(self, feedback=None):
        self.feedback = feedback

    def feedback_debug(self, text):
        if self.feedback is not None:
            # self.feedback.pushInfo(text)
            return
        LOGGER.debug(text)

    def feedback_info(self, text):
        if self.feedback is not None:
            self.feedback.pushInfo(text)
            return
        LOGGER.info(text)

    def feedback_warning(self, text):
        if self.feedback is not None:
            self.feedback.pushInfo(text)
            return
        LOGGER.warning(text)
