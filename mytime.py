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
