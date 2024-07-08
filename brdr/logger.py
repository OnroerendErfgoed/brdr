import logging

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S"
)


class Logger:
    def __init__(self, feedback=None):
        self.feedback = feedback

    def feedback_debug(self, text):
        if self.feedback is not None:
            # self.feedback.pushInfo(text)
            return
        logging.debug(text)

    def feedback_info(self, text):
        if self.feedback is not None:
            self.feedback.pushInfo(text)
            return
        logging.info(text)

    def feedback_warning(self, text):
        if self.feedback is not None:
            self.feedback.pushInfo(text)
            return
        logging.warning(text)
