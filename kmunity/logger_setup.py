#!/usr/bin/env python

"""...

"""

import sys
from loguru import logger
import kmunity


def colorize():
    """colorize the logger if stderr is IPython/Jupyter or a terminal (TTY)"""
    try:
        import IPython
        tty1 = bool(IPython.get_ipython())
    except ImportError:
        tty1 = False
    tty2 = sys.stderr.isatty()
    if tty1 or tty2:
        return True
    return False


LOGGERS = [0]


class Formatter:
    def __init__(self):
        self.padding = 0
        self.fmt = (
            "{level.icon} {module}:{function}{extra[padding]} | "
            "<level>{message}</level>\n{exception}"
        )

    def format(self, record):
        length = len("{module}:{function}".format(**record))
        self.padding = max(self.padding, length)
        record["extra"]["padding"] = " " * (self.padding - length)
        return self.fmt


def set_log_level(log_level="INFO"):
    """Set the log level for loguru logger bound to kmunity.
    """
    for idx in LOGGERS:
        try:
            logger.remove(idx)
        except ValueError:
            pass

    formatter = Formatter()
    if log_level in ("DEBUG", "TRACE"):
        idx = logger.add(
            sink=sys.stderr,
            level=log_level,
            colorize=colorize(),
            format=formatter.format,
            filter=lambda x: x['extra'].get("name") == "kmunity",
        )
    else:
        idx = logger.add(
            sink=sys.stderr,
            level=log_level,
            colorize=colorize(),
            format=formatter.format,
            filter=lambda x: x['extra'].get("name") == "kmunity",
        )
    LOGGERS.append(idx)
    logger.enable("kmunity")
    # logger.bind(name="kmunity").debug(
    #     f"phyloshape v.{kmunity.__version__} logging enabled"
    # )


if __name__ == "__main__":

    kmunity.set_log_level("DEBUG")
    logger = logger.bind(name="phyloshape")
    logger.info("HI")
    logger.warning("HIIIIII")

    # catch and raise exceptions with logger
    try:
        logger.fake()
    except AttributeError as exc:
        logger.exception("ERROR", exc)
