import time
from dataclasses import dataclass


@dataclass
class TimeInterval:
    start: int
    msg: str

    def __enter__(self):
        self.start = time.time()

    def __exit__(self, t, v, tr):
        elapsed = float(time.time() - self.start)
        if self.msg == "":
            print(f'Block took {elapsed:.1f}s')
        else:
            print(f'{self.msg} took {elapsed:.1f}s')


def timing(msg: str = ""):
    return TimeInterval(0, msg)


def is_binary(file: str):
    """ Based on https://stackoverflow.com/a/7392391 """
    textchars = {7, 8, 9, 10, 12, 13, 27} | set(range(0x20, 0x100)) - {0x7f}
    with open(file) as f:
        return any(ord(c) not in textchars for c in f.read(1024))
