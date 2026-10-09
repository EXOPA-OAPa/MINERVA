"""Redirection of the printed output to a log file."""
import sys


class ToggleableLogger(object):
    """Replacement for sys.stdout that writes to a log file and/or to the terminal."""

    def __init__(self, filename, to_file=True, to_screen=False):
        self.terminal = sys.stdout  # original terminal
        self.log = open(filename, "w", encoding='utf-8')   # new log at every run

        # Where the output goes
        self.log_to_file = to_file
        self.log_to_screen = to_screen

    def write(self, message):
        if self.log_to_file:
            self.log.write(message)
            self.log.flush()

        if self.log_to_screen:
            self.terminal.write(message)
            self.terminal.flush()

    def flush(self):
        self.log.flush()
        self.terminal.flush()

    def set_mode(self, to_file, to_screen):
        """Change where the output goes."""
        self.log_to_file = to_file
        self.log_to_screen = to_screen

    def restore_stdout(self):
        """Restore the original terminal and close the log file."""
        sys.stdout = self.terminal
        self.log.close()
